#!/usr/bin/env python3
"""
Cluster analysis utilities including adjacency computation.
"""

from typing import Dict, List
import scanpy as sc

from .logger import PipelineLogger

logger = PipelineLogger.get_logger(__name__)


def get_cluster_adjacency(adata: sc.AnnData, cluster_key: str) -> Dict[str, List[str]]:
    """
    Computes cluster adjacency using PAGA and formats it for downstream analysis.

    Args:
        adata: The annotated data object.
        cluster_key: The column name in adata.obs that contains cluster labels.

    Returns:
        A dictionary mapping each cluster to a list of its neighbors.
    """
    logger.info("Computing cluster adjacency using PAGA")
    
    sc.tl.paga(adata, groups=cluster_key)
    
    connectivity_matrix = adata.uns['paga']['connectivities']
    cluster_names = adata.obs[cluster_key].cat.categories
    
    adjacency_dict = {}
    for i, cluster_name in enumerate(cluster_names):
        # Find the indices of connected clusters in the sparse matrix
        connected_indices = connectivity_matrix[i].nonzero()[1]
        neighbor_names = [cluster_names[j] for j in connected_indices]
        adjacency_dict[str(cluster_name)] = neighbor_names
    
    logger.debug(f"Computed adjacency for {len(adjacency_dict)} clusters")
    return adjacency_dict
