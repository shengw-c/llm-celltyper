#!/usr/bin/env python3
"""
Data processing utilities for single-cell RNA-seq analysis.
Includes functions for normalization, dimensionality reduction, and clustering.
"""

import json
from typing import Dict, List, Optional, Tuple
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import colorcet as cc
from tqdm import tqdm
from tqdm import TqdmWarning
from pyclustree import clustree
import base64
import os
import warnings

from .logger import PipelineLogger

warnings.filterwarnings("ignore", category=TqdmWarning)

logger = PipelineLogger.get_logger(__name__)


def prepare_subset_dataset(
    nametag: str, 
    input_file: str, 
    cell_id_file: Optional[str] = None, 
    batch_key: Optional[List[str]] = None, 
    integration: bool = False,
    working_dir: str = "work",
    cpus: int = 16
) -> Dict[str, str]:
    """
    Prepares a subset AnnData object based on specified cell IDs, performs normalization,
    identifies highly variable genes, computes UMAP embeddings, and generates cluster trees.
    
    This function handles the complete preprocessing pipeline from raw counts to
    clustering analysis at multiple resolutions.
    
    Args:
        nametag: A tag to identify the dataset, will be used in output file names.
                Example: "lev0_Epithelial_Cell"
        input_file: Path to the input h5ad file. Must contain a 'counts' layer
                   with raw count data.
        cell_id_file: Path to a TSV file containing cell IDs to subset (one per line,
                     with cell IDs as index). If None, uses all cells from input_file.
        batch_key: List of column names in adata.obs for batch correction during
                  HVG selection and integration. None if no batch correction needed.
        integration: Whether to perform Harmony integration after PCA. If True,
                    requires batch_key to be specified.
        working_dir: Working directory for outputs. Creates 'figures' and 'data'
                    subdirectories.
        cpus: Number of CPUs to use for scanpy parallel operations.
    
    Returns:
        Dictionary containing:
            - 'encoded_cluster_tree': Base64-encoded PNG image of the cluster tree
            - 'subset_adata': Path to the saved h5ad file containing processed subset data
            - 'cluster_tree_file': Path to cluster tree image file
    
    Example:
        >>> dataset = prepare_subset_dataset(
        ...     nametag="lev0_test",
        ...     input_file="data/test.h5ad",
        ...     batch_key=["donor_id"],
        ...     integration=True,
        ...     cpus=8
        ... )
        >>> print(dataset['subset_adata'])
        'work/data/lev0_test_subset.h5ad'
    
    Notes:
        - Normalizes to total counts per cell, then log1p transforms
        - Identifies top 2000 highly variable genes
        - Computes PCA with default parameters (50 components)
        - Optionally performs Harmony integration
        - Clusters at resolutions from 0.0 to 0.95 in 0.05 increments
        - Generates cluster tree visualization showing resolution stability
    """
    # Setup directories
    fig_folder = os.path.join(working_dir, "figures")
    data_folder = os.path.join(working_dir, "data")
    os.makedirs(fig_folder, exist_ok=True)
    os.makedirs(data_folder, exist_ok=True)
    
    # Set scanpy parallelization
    sc.settings.n_jobs = cpus
    
    # Load the input h5ad file
    logger.info(f"Loading input file: {input_file}")
    adata = sc.read_h5ad(input_file)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Filter the AnnData object for the specified cell IDs
    if cell_id_file is not None:
        logger.info(f"Filtering cells using: {cell_id_file}")
        cell_id_df = pd.read_table(cell_id_file, index_col=0)
        cell_ids = cell_id_df.index.tolist()
        adata = adata[adata.obs.index.isin(cell_ids)].copy()
        logger.info(f"Filtered to {adata.n_obs} cells")

    # Normalization and preprocessing
    logger.info("Starting normalization and preprocessing")
    adata.X = adata.layers['counts'].copy()
    logger.debug("Copied raw counts into adata.X")
    
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    logger.debug("Normalization (normalize_total) complete")
    
    # Logarithmize the data
    sc.pp.log1p(adata)
    logger.debug("Log1p transformation complete")

    # Identify HVGs
    logger.info("Identifying highly variable genes")
    # If batch_key is a list with multiple keys, create a combined batch column
    hvg_batch_key = None
    if batch_key is not None:
        if isinstance(batch_key, list):
            if len(batch_key) == 1:
                hvg_batch_key = batch_key[0]
            elif len(batch_key) > 1:
                # Combine multiple batch keys into a single column
                logger.info(f"Combining batch keys {batch_key} into a single column for HVG detection")
                adata.obs['_combined_batch'] = adata.obs[batch_key].apply(
                    lambda row: '_'.join(row.astype(str)), axis=1
                ).astype('category')
                hvg_batch_key = '_combined_batch'
        else:
            hvg_batch_key = batch_key
    
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=hvg_batch_key, subset=False)
    hvgs_count = int(adata.var['highly_variable'].sum())
    logger.info(f"Identified {hvgs_count} highly variable genes")

    # Dimensionality reduction and clustering
    logger.info("Computing PCA")
    sc.tl.pca(adata)

    if integration and batch_key is not None:
        logger.info(f"Performing Harmony integration using batch key(s): {batch_key}")
        sc.external.pp.harmony_integrate(adata, key=batch_key) 
        logger.info("Computing neighborhood graph")
        sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    else:
        logger.info("Computing neighborhood graph")
        sc.pp.neighbors(adata)

    logger.info("Computing UMAP embedding")
    sc.tl.umap(adata)

    # Build leiden cluster trees using different resolutions
    logger.info("Computing Leiden clusters at multiple resolutions")
    resolutions = np.arange(0, 1, 0.05).round(2)
    for resolution in tqdm(resolutions, desc="Computing Leiden clusters"):
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=2,
            key_added=f"leiden_{str(resolution).replace('.', '_')}",
        )

    ## prepare the cluster labels
    labels = [f"leiden_{str(resolution).replace('.', '_')}" for resolution in resolutions]
    
    # Save cluster trees
    logger.info("Generating cluster tree visualization")
    cluster_tree_file = os.path.join(fig_folder, f"cluster_trees.{nametag}.png")
    fig = clustree(
        adata,
        labels,
        edge_weight_threshold=0.1,
    )
    fig.set_size_inches(10, 13)
    fig.set_dpi(100)
    fig.savefig(cluster_tree_file, bbox_inches='tight')
    logger.info(f"Cluster tree saved: {cluster_tree_file}")

    ## build adjacency matrices for LLM to pick up a good resolution
    # Dictionary to store all the flow matrices
    flow_matrices = {}
    df = adata.obs
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
    
    flow_matrices_json_ready = json.dumps(flow_matrices_json_ready, indent=2)
    logger.debug(f"Flow matrices JSON ready (length: {len(flow_matrices_json_ready)} characters)")

    # Save adata
    subset_adata_file = os.path.join(data_folder, f"{nametag}_subset.h5ad")
    adata.write_h5ad(subset_adata_file, compression='gzip')
    logger.info(f"Subset AnnData saved: {subset_adata_file}")

    return {
        "subset_adata": subset_adata_file,
        "cluster_adjacency_matrix": flow_matrices_json_ready,
        "cluster_tree_file": cluster_tree_file
    }


def generate_umap_visualization(
    adata: sc.AnnData,
    cluster_key: str,
    output_file: str
) -> str:
    """
    Generate UMAP visualization colored by clusters and save to file.
    
    Args:
        adata: AnnData object with UMAP coordinates
        cluster_key: Column name in adata.obs for cluster labels
        output_file: Path to save the PNG image
    
    Returns:
        Base64-encoded PNG image string
    """
    logger.info(f"Generating UMAP colored by {cluster_key}")
    
    umap = sc.pl.umap(
        adata, color=cluster_key, 
        legend_loc="on data", return_fig=True,
        palette=sns.color_palette(cc.glasbey, n_colors=50)
    )
    umap.set_size_inches(10, 10)
    umap.set_dpi(100)
    umap.savefig(output_file)
    logger.info(f"UMAP saved: {output_file}")


def prepare_celltype_inputs(
    nametag: str, 
    adata: sc.AnnData,
    resolution: float,
    working_dir: str = "work",
    gene_sets: List[str] = None,
    cpus: int = 16, 
    top_genes: int = 20,
    gsea_databases: str = "MSigDB_Hallmark_2020,KEGG_2021_Human",
    top_pathways: int = 20
) -> Dict:
    """
    Prepare inputs for cell type annotation including marker genes, pathways, and UMAP.
    
    This function extracts all features needed for LLM-based cell type annotation:
    differential gene expression, pathway enrichment, cluster relationships, and
    visualization.
    
    Args:
        nametag: Dataset identifier used in output filenames.
                Example: "lev0_Epithelial_Cell"
        adata: AnnData object (already processed with clusters). Must contain
              leiden clustering results at the specified resolution.
        resolution: Clustering resolution to use for annotation. Should match one
                   of the resolutions computed during prepare_subset_dataset.
                   Example: 0.15
        working_dir: Working directory for outputs. Saves UMAP to 'figures' subdirectory.
        gene_sets: List of gene set database names for pathway enrichment.
                  Default: ['MSigDB_Hallmark_2020', 'GO_Biological_Process_2025']
        cpus: Number of CPUs for parallel processing during pathway enrichment.
        top_genes: Number of top marker genes to extract per cluster (default: 20).
        gsea_databases: Comma-separated list of gene set databases for GSEA.
        top_pathways: Number of top enriched pathways to extract per cluster (default: 20).
    
    Returns:
        Dictionary containing:
            - 'top_genes_by_score': Top genes by rank_genes_groups scores
            - 'top_genes_by_specificity': Top genes by pct_in - pct_out
            - 'top_pathways': Top enriched pathways per cluster
            - 'consolidator_neighbors': Cluster adjacency from PAGA
            - 'encoded_umap': Base64-encoded UMAP visualization
            - 'umap_file': Path to saved UMAP PNG file
    
    Example:
        >>> inputs = prepare_celltype_inputs(
        ...     nametag="lev0_test",
        ...     adata=adata,
        ...     resolution=0.15,
        ...     top_genes=10,
        ...     top_pathways=15
        ... )
        >>> print(inputs.keys())
        dict_keys(['top_genes_by_score', ...])
    
    Notes:
        - Uses Wilcoxon rank-sum test for differential gene expression
        - Computes both score-based and specificity-based marker rankings
        - Runs GSEA prerank for pathway enrichment (positive NES only)
        - Computes PAGA to identify spatially adjacent clusters
        - All numeric values include fold changes and statistical significance
    """
    if gene_sets is None:
        gene_sets = ['MSigDB_Hallmark_2020', 'GO_Biological_Process_2025']
    
    fig_folder = os.path.join(working_dir, "figures")
    os.makedirs(fig_folder, exist_ok=True)

    sel_cls = f"leiden_{str(resolution).replace('.', '_')}"
    logger.info(f"Preparing inputs for cluster key: {sel_cls}, top_genes={top_genes}, top_pathways={top_pathways}")
    
    # Generate UMAP plot
    umap_file = f"{fig_folder}/umap.{nametag}.png"
    generate_umap_visualization(adata, sel_cls, umap_file)

    # Rank genes
    logger.info(f"Ranking genes for clusters using method=wilcoxon")
    sc.tl.rank_genes_groups(adata, sel_cls, method="wilcoxon", use_raw=False)
    
    # Get marker genes
    from .marker_genes import get_top_marker_genes
    logger.info(f"Extracting top {top_genes} marker genes per cluster")
    top_genes_by_score, top_genes_by_specificity = get_top_marker_genes(adata, sel_cls, n_genes=top_genes)
    logger.info(f"Extracted marker genes...")

    # Get enriched pathways
    from .pathway_enrichment import get_top_enriched_pathways
    # Parse gene_sets from gsea_databases if not provided
    if gene_sets is None:
        gene_sets = [db.strip() for db in gsea_databases.split(',')]
    logger.info(f"Running pathway enrichment (GSEA), extracting top {top_pathways} pathways per cluster")
    pathways = get_top_enriched_pathways(adata, sel_cls, gene_sets, n_pathways=top_pathways, cpus=cpus)
    logger.info(f"Extracted pathways...")

    # Get cluster adjacency
    from .cluster_analysis import get_cluster_adjacency
    logger.info("Computing cluster adjacency (PAGA)")
    adjacency_dict = get_cluster_adjacency(adata, sel_cls)
    logger.info(f"Computed adjacency...")

    logger.info("All inputs prepared successfully")
    return {
        "top_genes_by_score": top_genes_by_score,
        "top_genes_by_specificity": top_genes_by_specificity,
        "top_pathways": pathways,
        "adjacency": adjacency_dict,
        "umap_file": umap_file
    }
