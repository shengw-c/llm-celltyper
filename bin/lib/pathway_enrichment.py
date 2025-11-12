#!/usr/bin/env python3
"""
Pathway enrichment analysis utilities.
"""

from typing import Dict, List, Any
import scanpy as sc
import gseapy

from .logger import PipelineLogger

logger = PipelineLogger.get_logger(__name__)


def get_top_enriched_pathways(
    adata: sc.AnnData, 
    cluster_key: str, 
    gene_sets: List[str], 
    n_pathways: int = 10, 
    cpus: int = 4
) -> List[Dict[str, Any]]:
    """
    Runs GSEA prerank for each cluster and returns the top N enriched pathways.

    Args:
        adata: The annotated data object.
        cluster_key: The column name in adata.obs that contains cluster labels.
        gene_sets: A list of gene sets to use for GSEA (e.g., ['MSigDB_Hallmark_2020']).
        n_pathways: The number of top pathways to return.
        cpus: Number of CPUs to use for the analysis.

    Returns:
        A list of dictionaries, one for each cluster, 
        structured as [{'cluster': ..., 'pathways': [...]}]
    """
    logger.info(f"Running GSEA for top {n_pathways} pathways per cluster")
    
    ranked_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
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
                gene_sets=gene_sets,
                threads=cpus,
                verbose=False  # Keep the output clean
            )
            
            # Extract the top N pathway terms
            top_pathways = gsea_result.res2d.query("NES>0").head(n_pathways)['Term'].tolist()
            pathway_results.append({'cluster': str(cls), 'pathways': top_pathways})
            
        except Exception as e:
            logger.warning(f"GSEA failed for cluster {cls}: {str(e)}")
            pathway_results.append({'cluster': str(cls), 'pathways': []})
    
    logger.info(f"Completed GSEA for {len(pathway_results)} clusters")
    return pathway_results
