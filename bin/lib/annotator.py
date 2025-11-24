import scanpy as sc
import pandas as pd
import json
import os
from typing import Dict, List, Union, Optional, Tuple

from .logger import PipelineLogger
from .llm_client import (
    llm_cell_typer,
    llm_harmonize_cell_types
)
from .data_processing import DataProcessingHandler
from .cluster_id_manager import ClusterIDManager

logger = PipelineLogger.get_logger(__name__)


def reformat_cluster_data(
    cluster_datasets: Dict[str, Dict],
    context: Union[str, List[str]]
) -> List[Dict]:
    """
    Reformats cluster datasets for LLM annotation.
    
    Args:
        cluster_datasets: Dictionary mapping cluster IDs to their marker genes and pathways.
        context: Biological context string or list of context strings (one per cluster).
    
    Returns:
        List of dictionaries with cluster, context, marker_genes, and pathways.
    
    Raises:
        ValueError: If context is a list with length not matching number of clusters.
    """
    # Transform the Object into the required Array
    cluster_data_list = []
    if isinstance(context, list) and len(context) == len(cluster_datasets):
        for (cluster_id, data), ctx in zip(cluster_datasets.items(), context):
            cluster_data_list.append({
                "cluster": cluster_id,
                "context": ctx,
                "marker_genes": data["marker_genes"],
                "pathways": data["pathways"]
            })
    elif isinstance(context, str):
        for cluster_id, data in cluster_datasets.items():
            cluster_data_list.append({
                "cluster": cluster_id,
                "context": context,
                "marker_genes": data["marker_genes"],
                "pathways": data["pathways"]
            })
    else:
        raise ValueError("Context must be either a string or a list of strings matching the number of clusters.")
    return cluster_data_list

def updateJsonFile(json_file_path: str, new_data: Union[List[dict], dict]) -> None:
    """
    Updates a JSON file with new data. If the file does not exist, it creates a new one.

    Args:
        json_file_path: Path to the JSON file to be updated.
        new_data: List[dict] or dict containing the new data to be added or updated in the JSON file.
    """
    if os.path.exists(json_file_path):
        with open(json_file_path, 'r') as f:
            existing_data = json.load(f)
        # Update existing data with new data
        if isinstance(new_data, list):
            existing_data.extend(new_data)
        else:
            existing_data.append(new_data)
        # Write back to the JSON file
        with open(json_file_path, 'w') as f:
            json.dump(existing_data, f, indent=2)
    else:
        # Create a new JSON file with the new data
        data_to_write = new_data if isinstance(new_data, list) else [new_data]
        with open(json_file_path, 'w') as f:
            json.dump(data_to_write, f, indent=2)

    logger.info(f"Updated JSON file at: {json_file_path}")


def _run_initial_annotation(
    data_handler: DataProcessingHandler,
    general_context: str,
    gsea_databases: str,
    top_genes: int,
    top_pathways: int,
    max_resolution: float,
    transition_cutoff: float,
    llm_model: str,
    data_folder: str,
    response_json: str,
    niter: int
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Runs the initial (first iteration) annotation on all cells.
    
    Args:
        data_handler: Initialized DataProcessingHandler instance.
        general_context: Biological context for annotation.
        gsea_databases: GSEA databases to use.
        top_genes: Number of top marker genes.
        top_pathways: Number of top pathways.
        max_resolution: Maximum clustering resolution.
        transition_cutoff: Stability threshold.
        llm_model: LLM model name to use.
        data_folder: Directory to save cluster data.
        response_json: Path to response JSON file.
        niter: Current iteration number.
    
    Returns:
        Tuple of (annotation_results_df, cell2cluster_df).
    """
    logger.info("Processing initial annotation (all cells)")
    data_handler.reset_current_adata()
    
    # Prepare clustering inputs
    logger.info("Preparing dataset for annotation")
    optimal_resolution_res, cluster_datasets, cell2cluster_df = data_handler.prepare_celltype_inputs(
        gsea_databases=gsea_databases,
        top_genes=top_genes,
        top_pathways=top_pathways,
        max_resolution=max_resolution,
        transition_cutoff=transition_cutoff,
        parent_cluster_id=None,  # Initial annotation
        iteration=niter
    )
    
    cluster_data_list = reformat_cluster_data(
        cluster_datasets=cluster_datasets,
        context=general_context
    )
    
    # Store cluster data
    with open(os.path.join(data_folder, f"iter{niter}_cluster_data.json"), 'w') as f:
        json.dump({
            "optimal_resolution": optimal_resolution_res,
            "cluster_data": cluster_data_list
        }, f, indent=2)
    
    # Run LLM annotation
    logger.info("Running LLM-based cell type annotation")
    ann_res = llm_cell_typer(
        cluster_data=cluster_data_list,
        custom_model=llm_model
    )
    updateJsonFile(json_file_path=response_json, new_data=ann_res)
    
    # Update cluster ID manager with annotations
    if data_handler.cluster_id_manager:
        data_handler.cluster_id_manager.update_cluster_annotations(ann_res)
    
    # Rename columns for tracking
    cell2cluster_df.columns = ["cell_id", f"iter{niter}"]
    
    return pd.DataFrame(ann_res), cell2cluster_df


def _run_hierarchical_annotation(
    sub_cluster_data_list: List[Dict],
    sub_cell2cluster_mapping: pd.DataFrame,
    allcells: pd.DataFrame,
    llm_model: str,
    response_json: str,
    cluster_id_manager=None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Runs annotation on sub-clusters from previous iteration.
    
    Args:
        sub_cluster_data_list: List of cluster data prepared in previous iteration.
        sub_cell2cluster_mapping: Mapping of cells to sub-clusters.
        allcells: Master DataFrame tracking all cells.
        llm_model: LLM model name to use.
        response_json: Path to response JSON file.
        cluster_id_manager: Optional ClusterIDManager instance for updating annotations.
    
    Returns:
        Tuple of (annotation_results_df, updated_allcells_df).
    """
    logger.info(f"Processing hierarchical annotation ({len(sub_cluster_data_list)} sub-clusters)")
    
    # Run LLM annotation
    ann_res = llm_cell_typer(
        cluster_data=sub_cluster_data_list,
        custom_model=llm_model
    )

    updateJsonFile(json_file_path=response_json, new_data=ann_res)
    
    # Update cluster ID manager with annotations
    if cluster_id_manager:
        cluster_id_manager.update_cluster_annotations(ann_res)
    
    # Merge sub-cluster results back into allcells
    allcells_updated = allcells.merge(
        sub_cell2cluster_mapping,
        on="cell_id",
        how="left"
    )
    
    return pd.DataFrame(ann_res), allcells_updated


def _identify_clusters_for_splitting(
    ann_res_df: pd.DataFrame,
    allcells: pd.DataFrame,
    current_iter_col: str,
    min_cells_for_subtype: int
) -> pd.DataFrame:
    """
    Identifies clusters that need further subdivision.
    
    Args:
        ann_res_df: DataFrame containing annotation results.
        allcells: Master DataFrame tracking all cells.
        current_iter_col: Column name for current iteration (e.g., 'iter1').
        min_cells_for_subtype: Minimum cells required for subdivision.
    
    Returns:
        DataFrame containing clusters to split with their metadata.
    """
    logger.info("Checking which clusters require further subdivision")
    
    # Get clusters marked for splitting
    clusters_to_split = ann_res_df.query("split_status == 'Yes'")
    
    logger.info(f"Found {clusters_to_split.shape[0]} clusters marked for further sub-annotation")

    ## make sure the cluster has enough cells: based on cell type hypothesis
    passed_cluster = []
    passed_cell_types = []
    for cell_type in clusters_to_split.cell_type_hypothesis.unique():
        cell_type_rows = clusters_to_split[clusters_to_split.cell_type_hypothesis == cell_type]
        cell_type_ixs = cell_type_rows['cluster'].tolist()
        sel_cells = allcells[allcells[current_iter_col].isin(cell_type_ixs)]
        if sel_cells.shape[0] >= min_cells_for_subtype:
            passed_cluster.extend(cell_type_ixs)
            passed_cell_types.append(cell_type)
    further_split = clusters_to_split[clusters_to_split['cluster'].isin(passed_cluster)]
    logger.info(f"After applying minimum cell count filter ({min_cells_for_subtype}), {further_split.shape[0]} clusters ({len(passed_cell_types)} cell types) remain for subdivision")

    if further_split.shape[0] == 0:
        logger.info("No clusters marked for further sub-annotation")
        return pd.DataFrame()
    
    logger.info(f"Found {further_split.shape[0]} clusters meeting subdivision criteria")
    
    return further_split


def _prepare_subclusters(
    further_split: pd.DataFrame,
    allcells: pd.DataFrame,
    current_iter_col: str,
    niter: int,
    data_handler: DataProcessingHandler,
    gsea_databases: str,
    top_genes: int,
    top_pathways: int,
    max_resolution: float,
    transition_cutoff: float,
    data_folder: str
) -> Tuple[List[Dict], pd.DataFrame]:
    """
    Prepares sub-cluster data for next iteration.
    
    Args:
        further_split: DataFrame of clusters to subdivide.
        allcells: Master DataFrame tracking all cells.
        current_iter_col: Column name for current iteration.
        niter: Current iteration number.
        data_handler: DataProcessingHandler instance.
        gsea_databases: GSEA databases to use.
        top_genes: Number of top marker genes.
        top_pathways: Number of top pathways.
        max_resolution: Maximum clustering resolution.
        transition_cutoff: Stability threshold.
        data_folder: Directory to save cluster data.
    
    Returns:
        Tuple of (sub_cluster_data_list, sub_cell2cluster_mapping_df).
    """
    logger.info("Preparing data for next iteration of annotation")

    sub_cluster_data_list_combined = []
    sub_cell2cluster_mapping_list = []
    
    cell_types = further_split.cell_type_hypothesis.unique()

    for i in range(len(cell_types)):
        cell_type = cell_types[i]
        logger.info(f"Preparing sub-clusters for cell type: {cell_type} ({i+1}/{len(cell_types)})")
        cell_type_rows = further_split[further_split.cell_type_hypothesis == cell_type]
        cell_type_ixs = list(set(cell_type_rows["cluster"].tolist()))
        
        # Use the first parent cluster ID for tracking
        parent_cluster_id = cell_type_ixs[0]
        
        # Get cells belonging to this cluster
        sel_cells = allcells[allcells[current_iter_col].isin(cell_type_ixs)]["cell_id"].tolist()
        logger.info(f"{len(sel_cells)} cells selected for subdivision")
        
        # Reset data handler with selected cells
        data_handler.reset_current_adata(cell_ids=sel_cells)
        
        # Prepare sub-clustering inputs with parent tracking
        logger.info(f"Running subclustering for {cell_type}")
        sub_optimal_resolution_res, sub_cluster_datasets, sub_cell2cluster_df = data_handler.prepare_celltype_inputs(
            nametag=".".join(cell_type_ixs),
            gsea_databases=gsea_databases,
            top_genes=top_genes,
            top_pathways=top_pathways,
            max_resolution=max_resolution,
            transition_cutoff=transition_cutoff,
            parent_cluster_id=parent_cluster_id,
            iteration=niter+1
        )
        
        # Check if subclustering was successful
        if not sub_cluster_datasets:
            logger.warning(f"No suitable sub-clustering found for {cell_type}, skipping")
            continue
        
        # Get the context from the parent cluster annotation
        # Use next_round_context if available, otherwise use cell_type_hypothesis
        parent_context_row = cell_type_rows.iloc[0]
        if 'next_round_context' in parent_context_row and pd.notna(parent_context_row['next_round_context']):
            parent_context = parent_context_row['next_round_context']
        else:
            parent_context = f"{cell_type} subpopulation"
        
        # Create context for sub-clusters
        cell_type_contexts = [parent_context for _ in sub_cluster_datasets]
        
        # Reformat cluster data with appropriate context
        sub_cluster_data_list = reformat_cluster_data(
            cluster_datasets=sub_cluster_datasets,
            context=cell_type_contexts
        )
        sub_cluster_data_list_combined.extend(sub_cluster_data_list)
        
        # Track cell-to-subcluster mapping
        sub_cell2cluster_df.columns = ["cell_id", f"iter{niter+1}"]
        sub_cell2cluster_mapping_list.append(sub_cell2cluster_df)
    
    # Check if we have any sub-clusters to process
    if len(sub_cluster_data_list_combined) == 0:
        logger.info("No valid sub-clusters generated")
        return [], pd.DataFrame()
    
    # Combine all sub-cluster mappings
    sub_cell2cluster_mapping = pd.concat(sub_cell2cluster_mapping_list, ignore_index=True)
    
    # Store data for next iteration
    logger.info(further_split)
    with open(os.path.join(data_folder, f"iter{niter+1}_cluster_data.json"), 'w') as f:
        json.dump({
            "parent_clusters": further_split["cluster"].tolist(),
            "cluster_data": sub_cluster_data_list_combined
        }, f, indent=2)
    
    logger.info(f"Prepared {len(sub_cluster_data_list_combined)} sub-clusters for iteration {niter+1}")
    
    return sub_cluster_data_list_combined, sub_cell2cluster_mapping

def annotate_cell_types(
    general_context: str, 
    input_file: str, 
    layer: str = "counts",
    batch_key: Optional[List[str]] = None,
    integration: bool = False, 
    cpus_per_task: int = 16,
    min_cells_for_subtype: int = 1000,
    max_resolution: float = 1.0,
    transition_cutoff: float = 0.1,
    llm_model: str = "gemini-2.5-flash",
    llm_max_retries: int = 3,
    gsea_databases: str = "MSigDB_Hallmark_2020,KEGG_2021_Human",
    top_genes: int = 20,
    top_pathways: int = 20,
    output_dir: str = "."
) -> pd.DataFrame:
    """
    Performs hierarchical cell type annotation on single-cell RNA-seq data using LLM.
    
    Iteratively clusters cells and annotates cell types, recursively subdividing
    clusters marked for further refinement until no more splits are needed.
    
    Args:
        general_context: Biological context (e.g., "lung tissue, healthy adult").
        input_file: Path to input h5ad file with single-cell data.
        layer: Layer name containing raw counts (default: "counts").
        batch_key: Column names in adata.obs for batch correction.
        integration: Whether to perform Harmony integration.
        cpus_per_task: Number of CPUs for parallel processing.
        min_cells_for_subtype: Minimum cells required for subdivision.
        max_resolution: Maximum Leiden clustering resolution.
        transition_cutoff: Stability threshold for clustering.
        llm_model: LLM model name for annotation and harmonization.
        llm_max_retries: Maximum API retry attempts.
        gsea_databases: Comma-separated gene set databases.
        top_genes: Number of top marker genes per cluster.
        top_pathways: Number of top enriched pathways per cluster.
        output_dir: Output directory for results.
    
    Returns:
        DataFrame containing hierarchical cell annotations with iteration columns,
        cell type labels, confidence scores, and harmonized labels.
    """
    logger.info("=" * 80)
    logger.info(f"Starting annotation...")
    logger.info("=" * 80)

    # Set CPU usage for scanpy operations
    sc.settings.n_jobs = cpus_per_task
    
    # Create responses directory if it doesn't exist
    responses_dir = os.path.join(output_dir, "responses")
    os.makedirs(responses_dir, exist_ok=True)

    data_folder = os.path.join(output_dir, "data")
    os.makedirs(data_folder, exist_ok=True)
    
    response_json = os.path.join(responses_dir, "cell_annotation.json")

    # Initialize LLM client with configurable parameters
    # Note: LLM annotation code is currently commented out below
    # client = CellTypeAnnotationClient(
    #     max_retries=llm_max_retries
    # )

    # Initialize cluster ID manager
    cluster_id_mapping_file = os.path.join(data_folder, "cluster_id_mappings.json")
    cluster_id_manager = ClusterIDManager(mapping_file=cluster_id_mapping_file)
    
    # Initialize the data processing handler with cluster ID manager
    data_handler = DataProcessingHandler(
        input_file=input_file,
        layer=layer,
        batch_keys=batch_key,
        integration=integration,
        working_dir=output_dir,
        cpus=cpus_per_task,
        cluster_id_manager=cluster_id_manager
    )

    ## Perform nested annotation with LLM
    niter = 0
    allcells = None  # Track all cells across iterations
    max_iterations = 10  # Safety limit to prevent infinite loops
    
    # Variables for iteration continuity
    sub_cluster_data_list_combined = []
    sub_cell2cluster_mapping = pd.DataFrame()
    
    while niter < max_iterations:
        niter += 1
        logger.info(f"=" * 80)
        logger.info(f"Starting iteration {niter}")
        logger.info(f"=" * 80)
        
        # Run annotation based on iteration type
        if niter == 1:
            # Initial annotation on all cells
            ann_res_df, cell2cluster_df = _run_initial_annotation(
                data_handler=data_handler,
                general_context=general_context,
                gsea_databases=gsea_databases,
                top_genes=top_genes,
                top_pathways=top_pathways,
                max_resolution=max_resolution,
                transition_cutoff=transition_cutoff,
                llm_model=llm_model,
                data_folder=data_folder,
                response_json=response_json,
                niter=niter
            )
            allcells = cell2cluster_df.copy()
        else:
            # Hierarchical annotation on sub-clusters
            ann_res_df, allcells = _run_hierarchical_annotation(
                sub_cluster_data_list=sub_cluster_data_list_combined,
                sub_cell2cluster_mapping=sub_cell2cluster_mapping,
                allcells=allcells,
                llm_model=llm_model,
                response_json=response_json,
                cluster_id_manager=data_handler.cluster_id_manager
            )
        
        # Identify clusters that need further subdivision
        current_iter_col = f"iter{niter}"
        further_split = _identify_clusters_for_splitting(
            ann_res_df=ann_res_df,
            allcells=allcells,
            current_iter_col=current_iter_col,
            min_cells_for_subtype=min_cells_for_subtype
        )
        
        # Check termination conditions
        if further_split.shape[0] == 0:
            logger.info("Ending hierarchical annotation - no clusters require further subdivision")
            break
        
        # Prepare sub-clusters for next iteration
        sub_cluster_data_list_combined, sub_cell2cluster_mapping = _prepare_subclusters(
            further_split=further_split,
            allcells=allcells,
            current_iter_col=current_iter_col,
            niter=niter,
            data_handler=data_handler,
            gsea_databases=gsea_databases,
            top_genes=top_genes,
            top_pathways=top_pathways,
            max_resolution=max_resolution,
            transition_cutoff=transition_cutoff,
            data_folder=data_folder
        )
        
        # Check if preparation was successful
        if len(sub_cluster_data_list_combined) == 0:
            logger.info("Ending hierarchical annotation - no valid sub-clusters generated")
            break
    
    # Final output: save the complete hierarchical annotation
    logger.info("=" * 80)
    logger.info(f"Hierarchical annotation completed after {niter} iterations")
    logger.info("=" * 80)
    
    # Load all annotation responses as DataFrame
    logger.info("Loading annotation responses to add cell type information")
    if os.path.exists(response_json):
        with open(response_json, 'r') as f:
            all_annotations = json.load(f)
        annotations_df = pd.DataFrame(all_annotations)
        logger.info(f"Loaded {len(annotations_df)} annotation records")
    else:
        logger.warning("No annotation responses found")
        all_annotations = {}
        annotations_df = pd.DataFrame()
    
    # Determine the finest (deepest) cluster ID for each cell
    iter_cols = [col for col in allcells.columns if col.startswith('iter')]
    
    def get_final_cluster(row):
        """Get the last non-null cluster ID for a cell."""
        for col in reversed(iter_cols):  # Start from the deepest iteration
            if pd.notna(row[col]):
                return row[col]
        return None  # Should never happen as iter1 should always have a value
    
    allcells['final_cluster'] = allcells.apply(get_final_cluster, axis=1)
    logger.info(f"Determined final cluster for each cell")
    
    # Merge annotation information with allcells based on final cluster
    if not annotations_df.empty:
        # Rename 'cluster' column to 'final_cluster' for merging
        annotations_df = annotations_df.rename(columns={'cluster': 'final_cluster'})
        
        # Merge to add all annotation fields (cell_type, confidence, reasoning, etc.)
        allcells = allcells.merge(
            annotations_df,
            on='final_cluster',
            how='left',
            suffixes=('', '_annotation')
        )
        logger.info(f"Merged annotation information (cell_type, confidence, reasoning, etc.)")
    else:
        logger.warning("No annotations to merge")

    ## Harmonize cell type labels
    if all_annotations:
        to_harmonize = []
        final_cls = set(allcells["final_cluster"].unique())
        for ann in all_annotations:
            if ann["cluster"] in final_cls:
                to_harmonize.append({
                    "cluster": ann["cluster"],
                    "label": ann["cell_type_hypothesis"],
                    "definition": ann["cell_type_description"]
                })
        
        logger.info("Harmonizing cell type labels across clusters")
        harmonized_response_json = os.path.join(responses_dir, "cell_harmonization.json")
        harmonization_res = llm_harmonize_cell_types(
            cell_type_annotations=to_harmonize,
            custom_model=llm_model,
            max_retries=llm_max_retries
        )
        updateJsonFile(json_file_path=harmonized_response_json, new_data=harmonization_res)
        
        # Convert harmonization results to DataFrame
        harmonization_df = pd.DataFrame(harmonization_res)
        
        # Merge harmonized labels back into allcells
        if not harmonization_df.empty:
            allcells = allcells.merge(
                harmonization_df.rename(columns={"cluster": "final_cluster"}),
                on="final_cluster",
                how="left"
            )
            logger.info("Merged harmonized cell type labels")
        else:
            logger.warning("No harmonization results to merge")
    else:
        logger.warning("No annotations available for harmonization")
    
    # Save final cell-to-cluster mapping at all levels
    final_output_path = os.path.join(data_folder, "hierarchical_annotation_complete.csv")
    allcells.to_csv(final_output_path, index=False)
    logger.info(f"Saved complete hierarchical annotation to: {final_output_path}")
    
    # Generate UMAP visualizations with final annotations
    logger.info("Generating UMAP visualizations with final annotations")
    
    # Generate UMAP with harmonized cell types if available
    if 'harmonized_cell_type' in allcells.columns and allcells['harmonized_cell_type'].notna().any():
        logger.info("Generating UMAP visualization with harmonized cell types")
        ann_dict = allcells.set_index('cell_id')['harmonized_cell_type'].to_dict()
        umap_output_file = os.path.join(output_dir, "figures", "umap_harmonized_celltypes.png")
        data_handler.generate_umap_visualization(
            ann_dict=ann_dict,
            output_file=umap_output_file
        )
        logger.info(f"UMAP with harmonized labels saved to: {umap_output_file}")
    else:
        logger.warning("No harmonized labels found, skipping harmonized UMAP visualization")
    
    # At the end of annotation, export readable cluster ID table
    logger.info("Exporting cluster ID lineage table")
    cluster_id_table_file = os.path.join(data_folder, "cluster_id_lineage.csv")
    data_handler.cluster_id_manager.export_readable_table(cluster_id_table_file)
    
    # Export final annotated h5ad file
    logger.info("Exporting final annotated h5ad file with hierarchical annotations")
    final_h5ad_path = os.path.join(output_dir, "annotated_data.h5ad")
    data_handler.export_annotated_h5ad(
        allcells=allcells,
        output_file=final_h5ad_path
    )
    logger.info(f"Final annotated h5ad saved to: {final_h5ad_path}")
    
    logger.info("=" * 80)
    logger.info(f"Hierarchical annotation completed after {niter} iterations")
    logger.info("=" * 80)
    
    return allcells
