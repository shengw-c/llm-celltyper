#!/usr/bin/env python3
"""
Data processing utilities for single-cell RNA-seq analysis.
Includes functions for normalization, dimensionality reduction, and clustering.
"""

import json
from typing import Any, Dict, List, Optional, Tuple, Union
import gseapy
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import colorcet as cc
from tqdm import tqdm
from tqdm import TqdmWarning
import os
import warnings

from .logger import PipelineLogger
from .cluster_id_manager import ClusterIDManager

warnings.filterwarnings("ignore", category=TqdmWarning)

logger = PipelineLogger.get_logger(__name__)

class DataProcessingHandler:
    """
    Handles single-cell RNA-seq data processing for cell type annotation.
    
    Provides preprocessing, dimensionality reduction, clustering, and marker
    gene extraction for LLM-based cell type annotation.
    """
    def __init__(
        self,
        input_file: str,
        layer: str = "counts",
        batch_keys: Optional[List[str]] = None,
        integration: bool = False,
        working_dir: str = ".",
        cpus: int = 16,
        gsea_dbs: Optional[str] = "MSigDB_Hallmark_2020,KEGG_2021_Human",
        cluster_id_manager: Optional[ClusterIDManager] = None
    ):
        """
        Initializes the DataProcessingHandler.

        Args:
            input_file: Path to the input h5ad file.
            layer: Name of the layer in the AnnData object containing raw counts. Defaults to "counts".
            batch_keys: Column name(s) in adata.obs for batch correction. Can be a single key or a list of keys.
            integration: If True, perform Harmony integration for batch correction.
            working_dir: Directory to save output files (figures, data).
            cpus: Number of CPUs for parallel processing in scanpy.
            gsea_dbs: Comma-separated string of gene set databases for GSEA (e.g., "MSigDB_Hallmark_2020,KEGG_2021_Human").
        """
        self.nthreads = cpus
        
        sc.settings.n_jobs = self.nthreads

        self.layer = layer
        self.batch_keys = batch_keys
        self.integration = integration
        self.gsea_dbs = gsea_dbs.split(",") if gsea_dbs else []

        ## prepare working directories
        self.working_dir = working_dir
        self.fig_folder = os.path.join(working_dir, "figures")
        self.data_folder = os.path.join(working_dir, "data")
        self.response_folder = os.path.join(working_dir, "responses")
        os.makedirs(self.fig_folder, exist_ok=True)
        os.makedirs(self.data_folder, exist_ok=True)
        os.makedirs(self.response_folder, exist_ok=True)
    
        ## load the input h5ad file
        logger.info(f"Loading input file: {input_file}")
        self.adata = sc.read_h5ad(input_file)
        
        ### make sure the given batch_keys exists
        if self.batch_keys is not None:
            if isinstance(self.batch_keys, list):
                for key in self.batch_keys:
                    assert key in self.adata.obs.columns, f"Batch key '{key}' not found in adata.obs."
            else:
                assert self.batch_keys in self.adata.obs.columns, f"Batch key '{self.batch_keys}' not found in adata.obs."
        logger.info(f"Loaded {self.adata.n_obs} cells, {self.adata.n_vars} genes")

        ## define current status
        self.current_adata = None
        self.cluster_id_manager = cluster_id_manager

    def reset_current_adata(self, cell_ids: Optional[list] = None) -> None:
        """
        Resets the current AnnData object to the original loaded data or filters to specified cells.
        
        Args:
            cell_ids: Optional list of cell IDs to filter to. If None, uses all cells.
        """
        # Filter the AnnData object for the specified cell IDs
        if cell_ids is not None:
            logger.info(f"Filtering cells using {len(cell_ids)} provided cell IDs")
            self.current_adata = self.adata[self.adata.obs.index.isin(cell_ids)].copy()
            logger.info(f"Filtered to {self.current_adata.n_obs} cells")
        else:
            logger.info(f"No cell IDs provided, using all cells from input data")
            self.current_adata = self.adata.copy()

    def _determineStableResolutions(
        self,
        cluster_adjacency_matrix_dict: Dict[str, Any],
        transition_cutoff: float
    ) -> Tuple[Dict[str, int], Dict[str, int]]:
        """
        Identifies stable and substable clustering resolutions.

        Args:
            cluster_adjacency_matrix_dict: Cluster flow matrices between resolutions.
            transition_cutoff: Minimum flow proportion for significant transitions.

        Returns:
            Tuple of (stables, substables) where:
            - stables: Resolutions with no significant splitting or merging.
            - substables: Resolutions with no significant merging (splitting allowed).
        """
        stables = {}
        substables = {}
        for key in cluster_adjacency_matrix_dict.keys():
            ## skip the first layer to avoid trivial stable resolution
            if key=='leiden_0_0_to_leiden_0_05': continue
            dat = np.array(cluster_adjacency_matrix_dict[key]["data"])

            # skip trivial single cluster
            if dat.shape[0]==1: continue
            split_resolutions = ((dat>=transition_cutoff).sum(1)>1).any()
            merge_resolutions = ((dat>=transition_cutoff).sum(0)>1).any()
            if not (split_resolutions or merge_resolutions):
                stables[key] = dat.shape[0]
            if not merge_resolutions:
                substables[key] = dat.shape[0]
        logger.info(f"""Stable resolutions found: {stables}
Substable resolutions found: {substables}""")
        return stables, substables

    def _select_cluster_resolution(
        self,
        cluster_adjacency_matrix_dict: Dict[str, Any],
        transition_cutoff: float = 0.1,
        max_resolution: float = 1.0
    ) -> Dict[str, Union[float, str]]:
        """
        Selects optimal clustering resolution based on stability criteria.

        Args:
            cluster_adjacency_matrix_dict: Cluster flow matrices.
            transition_cutoff: Flow threshold for cluster stability.
            max_resolution: Maximum resolution that was actually computed.

        Returns:
            Dictionary with 'resolution' (float) and 'justification' (str).
        """     
        logger.info(f"Requesting cluster resolution selection...")
        
        ## try to do manual determination before sending to LLM
        logger.info(f"Try manual determination of cluster resolution")
        
        ## Four tiers of criteria:
        ## 1. stable or substable clusters meeting given resolution cutoff
        ## 2. stable or substable clusters meeting relaxed resolution cutoff (2 * given cutoff)
        ## STOP if any criteria is met
        ## Raise error if none of the criteria is met

        valid_stables, valid_substables = self._determineStableResolutions(cluster_adjacency_matrix_dict, transition_cutoff)
        
        if not valid_stables and not valid_substables:
            logger.info(f"No stable or subcluster identified using given transition cutoff of {transition_cutoff}, try to relax it to {2 * transition_cutoff}...")
            valid_stables, valid_substables = self._determineStableResolutions(cluster_adjacency_matrix_dict, transition_cutoff * 2)                    

        if valid_stables:
            valid_stables = {float(key.split("_to_")[0].replace("leiden_", "").replace("_", ".")):value for key,value in valid_stables.items()}
            sel_resolution = max(valid_stables.keys())
            logger.info(f"Manual determination successful based on stable clusters (maximum). Selected resolution: {sel_resolution} with {valid_stables[sel_resolution]} clusters.")
            return {
                "resolution": sel_resolution,
                "justification": f"Manual determination successful based on stable clusters (maximum). Selected resolution: {sel_resolution} with {valid_stables[sel_resolution]} clusters."
            }
        elif valid_substables:
            valid_substables = {float(key.split("_to_")[0].replace("leiden_", "").replace("_", ".")):value for key,value in valid_substables.items()}
            sel_resolution = max(valid_substables.keys())

            logger.info(f"Manual determination successful based on substable clusters (maximum). Selected resolution: {sel_resolution} with {valid_substables[sel_resolution]} clusters.")
            return {
                "resolution": sel_resolution,
                "justification": f"Manual determination successful based on substable clusters (maximum). Selected resolution: {sel_resolution} with {valid_substables[sel_resolution]} clusters."
            }
        else:
            logger.info(f"No stable resolutions found meeting target.")
            logger.info(f"Manual determination failed. Hard resetting to resolution {max_resolution}.")
            return {
                "resolution": max_resolution,
                "justification": f"No stable resolutions found meeting target. For safety, returning maximum resolution {max_resolution}."
            }


    def _get_top_gene_and_pathway(
        self,
        sel_cls: str,
        n_genes: int = 10,
        n_pathways: int = 10,
    ) -> Dict[str, Dict[str, Union[Dict[str, Dict[str, float]], List[str]]]]:
        """
        Extracts top marker genes and enriched pathways for each cluster.

        Args:
            sel_cls: Column name in adata.obs containing cluster labels.
            n_genes: Number of top marker genes per cluster.
            n_pathways: Number of top enriched pathways per cluster.

        Returns:
            Nested dictionary: {cluster_id: {"marker_genes": {...}, "pathways": [...]}}
        """
        
        logger.info(f"Preparing inputs for annotation: {sel_cls}, top_genes={n_genes} and top_pathways={n_pathways}..")

        ## prepare object to store the final results
        final_dict = {}

        # Get the full ranked genes DataFrame from the AnnData object
        ranked_genes_df = sc.get.rank_genes_groups_df(self.current_adata, group=None)
        
        # --- Efficiently calculate expression percentage ('pct') ---
        # Get a list of all unique genes across all clusters
        all_genes = ranked_genes_df['names'].unique()
        
        # Create a DataFrame with boolean expression values (True if gene is expressed)
        logger.debug("Calculating expression percentages")
        expr_df = pd.DataFrame(
            (self.current_adata[:, all_genes].X > 0).toarray(),
            index=self.current_adata.obs[sel_cls],
            columns=all_genes
        )
        
        # For each cluster, calculate the percentage of cells expressing each gene
        clusters = ranked_genes_df["group"].drop_duplicates()
        pct_df = []
        for cls in clusters:
            in_df = pd.DataFrame(expr_df.loc[cls,].mean().T).reset_index(names="names")
            out_df = pd.DataFrame(expr_df.drop(cls).mean(0).T).reset_index(names="names")
            cls_pct_df = in_df.merge(out_df, on="names")
            cls_pct_df.columns = ["names", "pct_in", "pct_out"]
            cls_pct_df["group"] = cls
            pct_df.append(cls_pct_df)
        pct_df = pd.concat(pct_df)

        # Add the percentages of expression back to the ranked_genes_df
        ranked_genes_df = ranked_genes_df.merge(pct_df)

        # top N based on scores from rank_genes_groups
        top_genes_df = ranked_genes_df.groupby('group', observed=True).head(n_genes).copy()

        # --- Reformat into the required nested dictionary structure ---
        # Output based on rank_genes_groups scores
        cols_to_keep = ['names', 'logfoldchanges', 'pvals_adj', 'pct_in', 'pct_out']

        ## round float values for cleaner output
        top_genes_df[["logfoldchanges", "pct_in", "pct_out"]] = top_genes_df[["logfoldchanges", "pct_in", "pct_out"]].round(3).astype(str).astype(float)
        top_genes_df["pvals_adj"] = top_genes_df["pvals_adj"].apply(lambda x: f"{x:.2e}").astype(str)

        for cluster, group_df in top_genes_df.groupby('group', observed=True):
            cluster_dict = {}
            for record in group_df[cols_to_keep].to_dict(orient='records'):
                gene_symbol = record.pop('names')
                cluster_dict[gene_symbol] = record
            final_dict[str(cluster)] = {"marker_genes": cluster_dict}

        # ----GSEA based on the ranked genes----
        clusters = ranked_genes_df['group'].unique()
        
        pathway_results = []
        for cls in clusters:
            logger.debug(f"Running GSEA for cluster {cls}")
            
            # Prepare the prerank input for the current cluster
            prerank_input = ranked_genes_df[ranked_genes_df['group'] == cls][
                ['names', 'logfoldchanges']
            ]

            ## upcase gene names for consistency
            prerank_input["names"] = prerank_input["names"].str.upper()
            
            # Run GSEA Prerank
            try:
                gsea_result = gseapy.prerank(
                    prerank_input,
                    gene_sets=self.gsea_dbs,
                    threads=self.nthreads,
                    verbose=False  # Keep the output clean
                )
                
                # Extract the top N pathway terms
                top_pathways = gsea_result.res2d.query("NES>0").head(n_pathways)['Term'].tolist()
                final_dict[str(cls)]['pathways'] = top_pathways
                
            except Exception as e:
                logger.warning(f"GSEA failed for cluster {cls}: {str(e)}")
                final_dict[str(cls)]['pathways'] = []
            
            ## uncomment to skip pathway analysis
            # final_dict[str(cls)]['pathways'] = []

        return final_dict

    def _determine_optimal_resolution(
        self,
        max_resolution: float = 1.0,
        transition_cutoff: float = 0.1,
    ) -> Dict[str, Union[float, str]]:
        """
        Determines optimal clustering resolution across multiple resolutions.

        Args:
            max_resolution: Maximum resolution for Leiden clustering.
            transition_cutoff: Flow threshold for cluster stability.

        Returns:
            Dictionary with 'resolution', 'justification', and 'cluster_tree_file'.
        """
        # Normalization and preprocessing
        if self.layer is None or self.layer not in self.adata.layers:
            logger.warning(f"'{self.layer}' not found, treat adata.X as raw counts.")
        else:
            logger.info(f"Using layer '{self.layer}' as raw counts.")
            self.current_adata.X = self.current_adata.layers[self.layer].copy()
            logger.debug("Copied raw counts into adata.X")
        
        # Normalizing to median total counts
        sc.pp.normalize_total(self.current_adata)
        logger.debug("Normalization (normalize_total) complete")
        
        # Logarithmize the data
        sc.pp.log1p(self.current_adata)
        logger.debug("Log1p transformation complete")

        # Identify HVGs
        logger.info("Identifying highly variable genes")
        # If batch_keys is a list with multiple keys, create a combined batch column
        hvg_batch_key = None
        if self.batch_keys is not None:
            if isinstance(self.batch_keys, list):
                if len(self.batch_keys) == 1:
                    hvg_batch_key = self.batch_keys[0]
                elif len(self.batch_keys) > 1:
                    # Combine multiple batch keys into a single column
                    logger.info(f"Combining batch keys {self.batch_keys} into a single column for HVG detection")
                    self.current_adata.obs['_combined_batch'] = self.current_adata.obs[self.batch_keys].apply(
                        lambda row: '_'.join(row.astype(str)), axis=1
                    ).astype('category')
                    hvg_batch_key = '_combined_batch'
            else:
                hvg_batch_key = self.batch_keys
        
        sc.pp.highly_variable_genes(self.current_adata, n_top_genes=2000, batch_key=hvg_batch_key, subset=False)
        hvgs_count = int(self.current_adata.var['highly_variable'].sum())
        logger.info(f"Identified {hvgs_count} highly variable genes")

        # Dimensionality reduction and clustering
        logger.info("Computing PCA")
        sc.tl.pca(self.current_adata)

        if self.integration and self.batch_keys is not None:
            logger.info(f"Performing Harmony integration using batch key(s): {self.batch_keys}")
            sc.external.pp.harmony_integrate(self.current_adata, key=self.batch_keys)
            logger.info("Computing neighborhood graph")
            sc.pp.neighbors(self.current_adata, use_rep="X_pca_harmony")
        else:
            logger.info("Computing neighborhood graph")
            sc.pp.neighbors(self.current_adata)

        logger.info("Computing UMAP embedding")
        sc.tl.umap(self.current_adata)

        ## if the current_data is actually the full data, adding the embeddings to the full data as well
        if self.current_adata.obs_names.equals(self.adata.obs_names):
            self.adata = self.current_adata.copy()
            logger.info("Updated the full AnnData object to include the embeddings")

        # Build leiden cluster trees using different resolutions
        logger.info(f"Computing Leiden clusters at multiple resolutions up to {max_resolution}")
        resolutions = np.arange(0, max_resolution+0.05, 0.05).round(2)
        for resolution in tqdm(resolutions, desc="Computing Leiden clusters"):
            sc.tl.leiden(
                self.current_adata,
                resolution=resolution,
                flavor="igraph",
                n_iterations=2,
                key_added=f"leiden_{str(resolution).replace('.', '_')}",
            )

        ## prepare the cluster labels
        labels = [f"leiden_{str(resolution).replace('.', '_')}" for resolution in resolutions]
        
        ## build adjacency matrices to pick up a good resolution
        # Dictionary to store all the flow matrices
        flow_matrices = {}
        df = self.current_adata.obs
        # Loop through adjacent pairs of columns
        for i in range(len(labels) - 1):
            col_low = labels[i]      # e.g., 'leiden_0_55'
            col_high = labels[i+1]   # e.g., 'leiden_0_6'
            
            # Create the contingency table
            adj_matrix = pd.crosstab(df[col_low], df[col_high], normalize="index").round(2)

            # Store it in the dictionary
            matrix_name = f"{col_low}_to_{col_high}"
            flow_matrices[matrix_name] = adj_matrix
        flow_matrices_json_ready = {}
        for matrix_name, matrix_df in flow_matrices.items():
            # "split" orientation is great because it's easy to reconstruct
            flow_matrices_json_ready[matrix_name] = matrix_df.to_dict(orient="split")
        
        ## determine the optimal resolution
        logger.info("Determining optimal clustering resolution")
        resolution_selection = self._select_cluster_resolution(
            cluster_adjacency_matrix_dict=flow_matrices_json_ready,
            transition_cutoff=transition_cutoff,
            max_resolution=max_resolution
        )

        return resolution_selection

    def prepare_celltype_inputs(
        self,
        nametag: Optional[str] = None,
        gsea_databases: str = "MSigDB_Hallmark_2020,KEGG_2021_Human",
        top_genes: int = 20,
        top_pathways: int = 20,
        max_resolution: float = 1.0,
        transition_cutoff: float = 0.1,
        parent_cluster_id: Optional[str] = None,
        iteration: int = 1
    ) -> Tuple[float, Dict[str, Dict], pd.DataFrame]:
        """
        Prepares comprehensive inputs for cell type annotation.

        Args:
            nametag: Optional prefix for cluster names.
            gsea_databases: Comma-separated GSEA database names.
            top_genes: Number of top marker genes per cluster.
            top_pathways: Number of top GSEA pathways per cluster.
            max_resolution: Maximum clustering resolution.
            transition_cutoff: Stability threshold for resolution selection.
            parent_cluster_id: Parent cluster ID (for subclusters).
            iteration: Current iteration number.

        Returns:
            Tuple of (optimal_resolution_result, cluster_data, cell2cluster_df).
            Returns empty dict for cluster_data if no suitable resolution found.
        """
        if gsea_databases is None:
            gsea_databases = ['MSigDB_Hallmark_2020', 'GO_Biological_Process_2025']

        ## determine optimal resolution cluster key
        ## if first run, use twice of the max_resolution to have more clusters
        max_resolution = max_resolution * 4 if nametag is None else max_resolution
        optimal_resolution_res = self._determine_optimal_resolution(
            transition_cutoff=transition_cutoff,
            max_resolution=max_resolution
        )
        if optimal_resolution_res["resolution"] == 0:
            return {}

        ## determine the current cluster key
        sel_cls = f"leiden_{str(optimal_resolution_res['resolution']).replace('.', '_')}"
        
        # Get original numeric cluster IDs from Leiden
        original_cluster_ids = sorted([str(x) for x in self.current_adata.obs[sel_cls].unique()])
        num_clusters = len(original_cluster_ids)
        
        # Generate random IDs directly
        if self.cluster_id_manager:
            # Generate random IDs for all clusters
            random_cluster_ids = self.cluster_id_manager.generate_multiple_ids(num_clusters)
            
            # Create mapping from original numeric IDs to random IDs
            id_mapping = dict(zip(original_cluster_ids, random_cluster_ids))
            
            # Register these clusters with metadata
            self.cluster_id_manager.register_clusters_batch(
                cluster_ids=random_cluster_ids,
                iteration=iteration,
                parent_id=parent_cluster_id
            )
            
            # Apply random IDs to the AnnData object
            self.current_adata.obs[sel_cls] = self.current_adata.obs[sel_cls].astype(str).map(id_mapping).astype('category')
            
            logger.info(f"Generated {num_clusters} random cluster IDs for iteration {iteration}")
        else:
            # Fallback: use numeric IDs with optional nametag
            if nametag is not None:
                cls_map = {str(i): f"{nametag}_{i}" for i in original_cluster_ids}
                self.current_adata.obs[sel_cls] = self.current_adata.obs[sel_cls].astype(str).map(cls_map).astype('category')

        # Rank genes
        num_clusters_for_ranking = len(self.current_adata.obs[sel_cls].unique())
        
        if num_clusters_for_ranking < 2:
            logger.warning(f"Only {num_clusters_for_ranking} cluster found. Skipping differential expression analysis.")
            # Return empty cluster data for single cluster
            cluster_data = {}
            cluster_id = str(self.current_adata.obs[sel_cls].unique()[0])
            cluster_data[cluster_id] = {
                "marker_genes": {},
                "pathways": []
            }
        else:
            logger.info(f"Ranking genes for clusters using method=wilcoxon")
            sc.tl.rank_genes_groups(self.current_adata, sel_cls, method="wilcoxon", use_raw=False)
            
            # Get marker genes and pathways
            cluster_data = self._get_top_gene_and_pathway(sel_cls=sel_cls, n_genes=top_genes, n_pathways=top_pathways)
        
        logger.info(f"Extracted marker genes and pathways...")

        ## get cell 2 cluster mapping
        if "cell_id" in self.current_adata.obs.columns:
            cell2cluster_df = self.current_adata.obs.drop(columns="cell_id").reset_index(names="cell_id")[["cell_id", sel_cls]]
        else:
            cell2cluster_df = self.current_adata.obs.reset_index(names="cell_id")[["cell_id", sel_cls]]

        return optimal_resolution_res, cluster_data, cell2cluster_df
    
    def generate_umap_visualization(
        self,
        ann_dict: Dict[str, str],
        output_file: str
    ) -> None:
        """
        Generates and saves a UMAP visualization colored by a specified annotation.

        Args:
            ann_dict: A dictionary mapping cell IDs to annotation labels.
            output_file: The path to save the output PNG image.
        """
        logger.info(f"Generating UMAP with annotation labels...")
        
        # Add labels to the full AnnData object (which has UMAP coordinates)
        self.adata.obs["label"] = self.adata.obs_names.map(ann_dict).astype('category')
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Generate UMAP plot
        umap = sc.pl.umap(
            self.adata, color="label", 
            legend_loc="on data", return_fig=True,
            palette=sns.color_palette(cc.glasbey, n_colors=50)
        )
        umap.set_size_inches(10, 10)
        umap.set_dpi(100)
        umap.savefig(output_file)
        logger.info(f"UMAP saved: {output_file}")
    
    def export_annotated_h5ad(
        self,
        allcells: pd.DataFrame,
        output_file: str
    ):
        """
        Export annotated AnnData object with hierarchical annotations and UMAP embeddings.
        
        Args:
            allcells: DataFrame containing hierarchical cluster assignments and annotations.
            output_file: Path to output h5ad file.
        """
        logger.info("Preparing annotated h5ad file for export")
        
        # Reset to full dataset (not just current subset)
        self.reset_current_adata()
        
        # Ensure UMAP has been computed
        if 'X_umap' not in self.current_adata.obsm:
            logger.info("Computing UMAP for visualization")
            sc.pp.neighbors(self.current_adata, use_rep='X_pca')
            sc.tl.umap(self.current_adata)
        
        # Merge all annotations from allcells into adata.obs
        # First, set cell_id as index for merging
        allcells_indexed = allcells.set_index('cell_id')
        
        # Get all columns to add (excluding cell_id which is now the index)
        annotation_cols = [col for col in allcells.columns if col != 'cell_id']
        
        # Add each annotation column to adata.obs
        for col in annotation_cols:
            if col in allcells_indexed.columns:
                # Convert list columns to strings for HDF5 compatibility
                if col in ['key_markers_cited', 'key_pathways_cited']:
                    # Convert lists to comma-separated strings
                    self.current_adata.obs[col] = allcells_indexed.loc[self.current_adata.obs.index, col].apply(
                        lambda x: ','.join(x) if isinstance(x, list) else (x if pd.notna(x) else '')
                    )
                else:
                    self.current_adata.obs[col] = allcells_indexed.loc[self.current_adata.obs.index, col]
                logger.info(f"Added annotation column: {col}")
        
        # Add metadata about the annotation process
        self.current_adata.uns['annotation_metadata'] = {
            'annotation_date': pd.Timestamp.now().isoformat(),
            'total_iterations': len([col for col in annotation_cols if col.startswith('iter')]),
            'has_harmonized_labels': 'harmonized_cell_type' in annotation_cols,
            'annotation_columns': annotation_cols
        }
        
        # Write the annotated h5ad file
        logger.info(f"Writing annotated h5ad to: {output_file}")
        self.current_adata.write_h5ad(output_file)
        
        logger.info(f"Successfully exported annotated h5ad with {len(annotation_cols)} annotation columns")
        logger.info(f"File size: {os.path.getsize(output_file) / 1e9:.2f} GB")

