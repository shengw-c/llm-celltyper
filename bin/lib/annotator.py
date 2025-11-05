#!/usr/bin/env python3
"""
Main cell type annotation module with hierarchical annotation support.
"""

import scanpy as sc
import pandas as pd
import json
import os
from typing import Dict, List, Union, Optional

from .logger import PipelineLogger
from .llm_client import CellTypeAnnotationClient
from .data_processing import prepare_subset_dataset, prepare_celltype_inputs
from .tree_utils import find_node_data, get_immediate_children, TreeValidationError
from .prompts import cluster_PROMPT, Celltyper_Instruction

logger = PipelineLogger.get_logger(__name__)


class HierarchicalAnnotation:
    """
    Manages hierarchical cell type annotations across multiple levels.
    
    This class collects and organizes annotations from all levels of the
    cell type hierarchy, allowing for comprehensive tracking and export
    of the complete annotation tree.
    
    Attributes:
        output_dir: Directory for output files
        annotations: Dictionary mapping nametag to annotation DataFrame
        hierarchy_levels: Dictionary mapping nametag to hierarchy level number
    
    Example:
        >>> collector = HierarchicalAnnotation()
        >>> collector.add_annotation("lev0_test", df_lev0, level=0)
        >>> collector.add_annotation("lev1_epithelial", df_lev1, level=1)
        >>> collector.export_final_annotations("final.tsv")
    """
    
    def __init__(self, output_dir: str = "work"):
        """
        Initialize the hierarchical annotation collector.
        
        Args:
            output_dir: Directory for output files (default: "work")
        """
        self.output_dir = output_dir
        self.annotations = {}  # {nametag: annotation_dataframe}
        self.hierarchy_levels = {}  # {nametag: level_number}
    
    def add_annotation(self, nametag: str, annotation_df: pd.DataFrame, level: int):
        """
        Add an annotation result to the collection.
        
        Args:
            nametag: Unique identifier for this annotation run
            annotation_df: DataFrame with cell annotations (index=cell_ids)
            level: Hierarchy level number (0=root, 1=first children, etc.)
        """
        self.annotations[nametag] = annotation_df
        self.hierarchy_levels[nametag] = level
        logger.debug(f"Added annotation for {nametag} at level {level}")
    
    def export_final_annotations(self, output_file: str = None) -> str:
        """
        Export all hierarchical annotations to a single TSV file in pivoted format.
        
        Creates a wide-format table where each hierarchy level becomes a column.
        The finest (deepest) annotation for each cell is also provided.
        
        Args:
            output_file: Path to output file (default: work/final_annotations.tsv)
        
        Returns:
            Path to the exported file
        
        Example:
            Output file format:
            ```
            cell_id              ann_level0       ann_level1      ann_level2     ann_finest
            AAACCTGAGCGATATA_1  Epithelial Cell  Alveolar Cell   AT2 Cell       AT2 Cell
            BBACCTGAGCGATATA_1  Immune Cell      T Cell          CD8 T Cell     CD8 T Cell
            ```
        """
        if output_file is None:
            output_file = os.path.join(self.output_dir, "final_annotations.tsv")
        
        logger.info(f"Exporting final annotations to {output_file}")
        
        if not self.annotations:
            logger.warning("No annotations to export")
            return output_file
        
        # Group annotations by level
        levels_data = {}
        max_level = 0
        for nametag, df in self.annotations.items():
            level = self.hierarchy_levels.get(nametag, 0)
            max_level = max(max_level, level)
            
            # Create a simplified DataFrame with just cell_id and annotation
            temp_df = df[['ann']].copy()
            temp_df.columns = [f'ann_level{level}']
            
            if level not in levels_data:
                levels_data[level] = []
            levels_data[level].append(temp_df)
        
        # Merge all DataFrames at each level
        level_dfs = {}
        for level in sorted(levels_data.keys()):
            # Concatenate all annotations at this level, keeping the last one for duplicate cells
            combined = pd.concat(levels_data[level], axis=0)
            # Remove duplicates, keeping the last annotation
            level_dfs[level] = combined[~combined.index.duplicated(keep='last')]
        
        # Start with the first level
        if 0 in level_dfs:
            result_df = level_dfs[0]
        else:
            # If no level 0, start with the lowest available level
            min_level = min(level_dfs.keys())
            result_df = level_dfs[min_level]
        
        # Progressively join all levels
        for level in sorted(level_dfs.keys()):
            if level == 0 or (level == min_level if 0 not in level_dfs else False):
                continue
            result_df = result_df.join(level_dfs[level], how='outer')
        
        # Add cell_id as first column
        result_df.insert(0, 'cell_id', result_df.index)
        
        # Add ann_finest column - the deepest annotation for each cell
        ann_cols = [col for col in result_df.columns if col.startswith('ann_level')]
        result_df['ann_finest'] = result_df[ann_cols].ffill(axis=1).iloc[:, -1]

        # Save to file
        result_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Exported {len(result_df)} cells with {len(ann_cols)} annotation levels to {output_file}")
        
        return output_file
    
    def get_level_annotations(self, level: int) -> Dict[str, pd.DataFrame]:
        """
        Get all annotations at a specific hierarchy level.
        
        Args:
            level: Hierarchy level to retrieve (0=root, 1=children, etc.)
        
        Returns:
            Dictionary mapping nametag to annotation DataFrame for that level
        """
        return {
            nametag: df 
            for nametag, df in self.annotations.items() 
            if self.hierarchy_levels.get(nametag, 0) == level
        }


def annotate_cell_types(
    expected_cells: Dict[str, Dict[str, Union[bool, str]]], 
    tree_file: str, 
    general_context: str, 
    nametag: str, 
    input_file: str, 
    batch_key: Optional[List[str]] = None,
    integration: bool = False, 
    subset_cell_id_file: Optional[str] = None,
    cpus_per_task: int = 16,
    current_level: int = 0,
    hierarchical_collector: Optional[HierarchicalAnnotation] = None,
    min_cells_for_subtype: int = 1000,
    condition_num: int = 1,
    llm_model_general: str = "gemini-2.5-flash",
    llm_model_complicated: str = "gemini-2.5.pro",
    llm_max_retries: int = 3,
    llm_nreplies: int = 1,
    gsea_databases: str = "MSigDB_Hallmark_2020,KEGG_2021_Human",
    top_genes: int = 20,
    top_pathways: int = 20,
    nested: bool = True,
    output_dir: str = "."
) -> HierarchicalAnnotation:
    """
    Hierarchically annotates cell types in single-cell RNA-seq data using LLM-based annotation.
    
    This function performs automated cell type annotation at the current hierarchical level,
    then recursively annotates subtypes for any cell types that have children in the tree.
    
    Args:
        expected_cells: Dictionary of expected cell types at this level.
                       Format: {cell_name: {'has_children': bool, 'definition': str}}
        tree_file: Path to JSON file containing the complete cell type hierarchy tree.
        general_context: Biological context description (e.g., "lung tissue, healthy adult").
        nametag: Unique identifier for this annotation run, used in output filenames.
        input_file: Path to the input h5ad file containing the single-cell data.
        batch_key: List of column names in adata.obs for batch correction. None if no batches.
        integration: Whether to perform Harmony integration.
        subset_cell_id_file: Path to TSV file with cell IDs to subset. If None, uses all cells.
        cpus_per_task: Number of CPUs to use for parallel processing.
        current_level: Current hierarchy level (0 = root).
        hierarchical_collector: Collector for hierarchical annotations.
        min_cells_for_subtype: Minimum number of cells required to annotate subtypes (default: 1000).
        condition_num: Number of expected condition-driven groups or artifacts (default: 1).
        llm_model_general: Name of the LLM model for general purpose (default: "gemini-2.5-flash").
        llm_model_complicated: Name of the LLM model for complicated tasks (default: "gemini-2.5.pro").
        llm_nreplies: Number of replies to request from LLM (default: 1, used for cell type annotation).
        llm_max_retries: Maximum retries for LLM API calls (default: 3).
        gsea_databases: Comma-separated list of gene set databases (default: "MSigDB_Hallmark_2020,KEGG_2021_Human").
        top_genes: Number of top marker genes per cluster (default: 20).
        top_pathways: Number of top pathways per cluster (default: 20).
        nested: Whether to perform recursive annotation of subtypes (default: True).
               If False, only annotates current level without recursion.
        output_dir: Output directory for results (default: ".").
    
    Returns:
        HierarchicalAnnotation object containing all annotation results.
    """
    # Initialize collector if not provided
    if hierarchical_collector is None:
        hierarchical_collector = HierarchicalAnnotation(output_dir=output_dir)
    
    logger.info("=" * 80)
    logger.info(f"Starting annotation for: {nametag} (Level {current_level})")
    logger.info("=" * 80)
    
    # Set CPU usage for scanpy operations
    sc.settings.n_jobs = cpus_per_task
    
    # Create responses directory if it doesn't exist
    responses_dir = os.path.join(output_dir, "responses")
    os.makedirs(responses_dir, exist_ok=True)
    
    # Initialize LLM client with configurable parameters
    client = CellTypeAnnotationClient(
        model_name=llm_model_general,
        max_retries=llm_max_retries
    )
    
    ## Prepare the input for clustree
    logger.info("Preparing subset dataset and computing clusters")
    dataset = prepare_subset_dataset(
        nametag=nametag, 
        input_file=input_file,
        cell_id_file=subset_cell_id_file,
        batch_key=batch_key,
        integration=integration,
        working_dir=output_dir,
        cpus=cpus_per_task
    )
    
    ## Determine the good resolution
    logger.info("Selecting optimal clustering resolution")
    cls_contents = [
        cluster_PROMPT.format(
            cell_type_num=len(expected_cells),
            condition_num=condition_num,
            transition_cutoff=0.1
        ),
        dataset["cluster_adjacency_matrix"]
    ]

    cls_response = client.select_cluster_resolution(cls_contents)

    # Save cluster selection response
    cluster_response_file = os.path.join(responses_dir, f"{nametag}_cluster_selection.json")
    with open(cluster_response_file, 'w') as f:
        json.dump({
            "nametag": nametag,
            "level": current_level,
            "timestamp": pd.Timestamp.now().isoformat(),
            "expected_cell_types": list(expected_cells.keys()),
            "expected_cell_type_num": len(expected_cells),
            "cluster_selection_response": cls_response,
            "cluster_tree_file": dataset.get("cluster_tree_file", ""),
        }, f, indent=2)
    logger.info(f"Saved cluster selection response: {cluster_response_file}")
    
    ## Prepare input data for annotation
    # Load the adata object once here
    logger.info(f"Loading processed AnnData from {dataset['subset_adata']}")
    adata = sc.read_h5ad(dataset["subset_adata"])
    
    logger.info("Preparing cell type annotation inputs")
    ann_input = prepare_celltype_inputs(
        nametag=nametag, 
        adata=adata,  # Pass the AnnData object directly
        resolution=cls_response["resolution"],
        working_dir=output_dir,
        cpus=cpus_per_task,
        gsea_databases=gsea_databases,
        top_genes=top_genes,
        top_pathways=top_pathways
    )
    
    ## Run annotation
    logger.info("Running LLM-based cell type annotation")
    # Use the new simplified dictionary keys
    formated_prompt = Celltyper_Instruction.format(
        expr_context=general_context,
        candidate_cell_types=json.dumps(expected_cells),
        marker_genes_json=json.dumps(ann_input["top_genes_by_specificity"]),
        pathway_json=json.dumps(ann_input["top_pathways"]),
        cluster_adjacency_json=json.dumps(ann_input["adjacency"])
    )

    logger.debug(f"Annotation prompt: {formated_prompt}")
    ann_dict, res_df_dict = client.annotate_cell_types(
        formated_prompt, 
        nreplicate=llm_nreplies,
        custom_model=llm_model_complicated,
        custom_model_consolidation=llm_model_general
    )

    # Save cell type annotation response
    if not res_df_dict is None:
        # Save replicate responses if available
        replicate_annotation_response_file = os.path.join(responses_dir, f"{nametag}_cell_annotation_by_replicate.json")
        with open(replicate_annotation_response_file, 'w') as f:
            json.dump(res_df_dict, f, indent=2)

    annotation_response_file = os.path.join(responses_dir, f"{nametag}_cell_annotation.json")
    with open(annotation_response_file, 'w') as f:
        json.dump({
            "nametag": nametag,
            "level": current_level,
            "timestamp": pd.Timestamp.now().isoformat(),
            "general_context": general_context,
            "expected_cell_types": list(expected_cells.keys()),
            "marker_genes_used": ann_input["top_genes_by_specificity"],
            "pathways_used":  ann_input["top_pathways"],
            "resolution": cls_response["resolution"],
            "annotation_response": ann_dict,
            "umap_file": ann_input.get("umap_file", ""),
        }, f, indent=2)
    logger.info(f"Saved annotation response: {annotation_response_file}")

    ## Make annotation pandas dataframe
    logger.info("Processing annotation results")
    # Use the previously loaded AnnData object to avoid redundant I/O
    cluster_col = f'leiden_{str(cls_response["resolution"]).replace(".", "_")}'
    df = adata.obs[[cluster_col]].copy()

    # Convert ann_dict to DataFrame
    ann_df = pd.DataFrame(ann_dict)
    try:
        ann_dict = ann_df.set_index("cluster_id")["cell_type"].to_dict()
    except KeyError as e:
        logger.info("No consensus 'cell_type' found in annotation, falling back to 'cell_type_hypotheses'")
        ann_dict = ann_df.set_index("cluster_id")["cell_type_hypotheses"].to_dict()
    df["ann"] = adata.obs[cluster_col].map(ann_dict)
    
    # Add to hierarchical collector
    hierarchical_collector.add_annotation(nametag, df, current_level)
    
    # Save annotation results for this level
    output_tsv = dataset["subset_adata"].replace(".h5ad", ".tsv")
    df.to_csv(output_tsv, sep="\t")
    logger.info(f"Saved annotation results: {output_tsv}")
    
    ## Recursive annotation for subtypes (only if nested=True)
    if not nested:
        logger.info("Nested annotation disabled (--no-nest flag), splitting dataset by cell type")
        
        # Create split directory
        split_dir = os.path.join(output_dir, "split")
        os.makedirs(split_dir, exist_ok=True)
        
        # Split dataset into separate files for each annotated cell type
        split_files = []
        for cell_type in expected_cells.keys():
            # Filter cells for this type
            cells_in_type = df[df["ann"] == cell_type].index.tolist()
            
            if len(cells_in_type) == 0:
                logger.warning(f"No cells found for cell type '{cell_type}', skipping")
                continue
            
            logger.info(f"Splitting {len(cells_in_type)} cells for '{cell_type}'")
            
            # Subset the AnnData object
            adata_subset = adata[cells_in_type, :].copy()
            
            # Save to file
            cell_type_safe = cell_type.replace(' ', '_').replace('/', '_')
            split_file = os.path.join(split_dir, f"{cell_type_safe}.h5ad")
            adata_subset.write_h5ad(split_file, compression='gzip')
            split_files.append(split_file)
            logger.info(f"Saved split file: {split_file}")
        
        logger.info(f"Split complete: created {len(split_files)} files")
        logger.info(f"Completed annotation for {nametag}")
        return hierarchical_collector
    
    logger.info("Checking for subtypes to annotate")
    with open(tree_file, 'r') as f:
        cell_tree = json.load(f)
    
    for key in expected_cells.keys():
        if expected_cells[key]["has_children"]:
            logger.info(f"Cell type '{key}' has children, preparing subtype annotation")
            
            # Filter cell IDs for this specific cell type
            cells_in_current_type = df[df["ann"] == key].index.tolist()
            
            if len(cells_in_current_type) < min_cells_for_subtype:
                logger.warning(
                    f"Only {len(cells_in_current_type)} cells found for cell type '{key}' "
                    f"(minimum {min_cells_for_subtype} required), skipping subtype annotation"
                )
                continue
            
            logger.info(f"Found {len(cells_in_current_type)} cells for '{key}' subtype annotation")
            
            # Create a temporary file with cell IDs for this subtype
            sub_cell_id_file = os.path.join(output_dir, "tmp", f"{nametag}_{key}_cell_ids.tsv")
            os.makedirs(os.path.dirname(sub_cell_id_file), exist_ok=True)
            sub_cell_df = pd.DataFrame(index=cells_in_current_type)
            sub_cell_df.to_csv(sub_cell_id_file, sep="\t")
            
            # Prepare recursive call parameters
            sub_node_data = find_node_data(cell_name=key, tree_dict=cell_tree, raise_on_missing=True)
            # include parent cell type, in case no clear specific children identified
            sub_expected_cells = get_immediate_children(sub_node_data, include_parent=False, parent_name=key)
            
            if not sub_expected_cells:
                logger.warning(
                    f"Cell type '{key}' has no valid children defined in tree, skipping subtype annotation"
                )
                continue
                
            # Build context properly without duplication
            # Check if we're already in a focused context (to avoid nested "similar to")
            if "Focus on the cell type:" in general_context:
                # Already in a focused context, just update the focus
                base_context = general_context.split('Focus on the cell type:')[0].strip().rstrip('.')
                sub_general_context = f"{base_context}. Focus on the cell type: {key}."
            else:
                # First level of focus
                base_context = general_context.strip().rstrip('.')
                sub_general_context = f"{base_context}. Focus on the cell type: {key}."
            sub_nametag = f"{nametag}_{key.replace(' ', '_')}"
            sub_input_file = dataset["subset_adata"]
            
            # Recursive call for subtype annotation
            logger.info(f"Recursively annotating subtypes of '{key}'")
            annotate_cell_types(
                expected_cells=sub_expected_cells,
                tree_file=tree_file,
                general_context=sub_general_context,
                nametag=sub_nametag,
                input_file=sub_input_file,
                batch_key=batch_key,
                integration=integration,  # Preserve parent's integration setting
                subset_cell_id_file=sub_cell_id_file,
                cpus_per_task=cpus_per_task,
                current_level=current_level + 1,
                hierarchical_collector=hierarchical_collector,
                min_cells_for_subtype=min_cells_for_subtype,
                condition_num=condition_num,
                llm_model_general=llm_model_general,
                llm_model_complicated=llm_model_complicated,
                llm_max_retries=llm_max_retries,
                llm_nreplies=llm_nreplies,
                gsea_databases=gsea_databases,
                top_genes=top_genes,
                top_pathways=top_pathways,
                output_dir=output_dir
            )
    
    logger.info(f"Completed annotation for {nametag}")
    return hierarchical_collector
