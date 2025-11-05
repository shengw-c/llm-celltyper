#!/usr/bin/env python3
"""
Cluster analysis utilities including adjacency computation.
"""

from typing import Dict, List
import scanpy as sc
import pandas as pd
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
    
    connectivity_matrix = adata.uns['paga']['connectivities'].toarray().round(4)
    connectivity_matrix = pd.DataFrame(connectivity_matrix, 
                                       index=adata.obs[cluster_key].cat.categories,
                                       columns=adata.obs[cluster_key].cat.categories)
    connectivity_matrix_json = connectivity_matrix.to_dict(orient="split")
    
    return connectivity_matrix_json
