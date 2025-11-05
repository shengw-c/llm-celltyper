#!/usr/bin/env python3
"""
Marker gene identification and extraction utilities.
"""

from typing import Dict
import scanpy as sc
import pandas as pd

from .logger import PipelineLogger

logger = PipelineLogger.get_logger(__name__)


def get_top_marker_genes(
    adata: sc.AnnData, 
    cluster_key: str, 
    n_genes: int = 10
) -> tuple[Dict[str, Dict[str, Dict[str, float]]], Dict[str, Dict[str, Dict[str, float]]]]:
    """
    Extracts top N marker genes per cluster, calculates expression percentages,
    and formats the data into nested dictionaries.

    Args:
        adata: The annotated data object.
        cluster_key: The column name in adata.obs that contains cluster labels.
        n_genes: The number of top genes to extract for each cluster.

    Returns:
        Tuple of two dictionaries:
            1. Genes ranked by rank_genes_groups scores
            2. Genes ranked by pct_diff (pct_in - pct_out)
    """
    logger.info(f"Extracting top {n_genes} marker genes for each cluster")
    
    # Get the full ranked genes DataFrame from the AnnData object
    ranked_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
    
    # --- Efficiently calculate expression percentage ('pct') ---
    # Get a list of all unique genes across all clusters
    all_genes = ranked_genes_df['names'].unique()
    
    # Create a DataFrame with boolean expression values (True if gene is expressed)
    logger.debug("Calculating expression percentages")
    expr_df = pd.DataFrame(
        (adata[:, all_genes].X > 0).toarray(),
        index=adata.obs[cluster_key],
        columns=all_genes
    )
    
    # For each cluster, calculate the percentage of cells expressing each gene
    clusters = ranked_genes_df.group.drop_duplicates()
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

    # Two ways to filter top N genes per cluster
    # 1. top N based on scores from rank_genes_groups
    top_genes_df1 = ranked_genes_df.groupby('group', observed=True).head(n_genes).copy()
    
    # 2. top N based on pct_in - pct_out, after filtering significant genes
    ranked_genes_df['pct_diff'] = ranked_genes_df['pct_in'] - ranked_genes_df['pct_out']
    top_genes_df2 = (
        ranked_genes_df
        .query("logfoldchanges>0.5 and pvals_adj<0.05")
        .sort_values(['group', 'pct_diff'], ascending=[True, False])
        .groupby('group', observed=True)
        .head(n_genes)
        .copy()
    )

    # --- Reformat into the required nested dictionary structure ---
    # Output 1 : based on rank_genes_groups scores
    final_gene_dict1 = {}
    cols_to_keep = ['names', 'logfoldchanges', 'pvals_adj', 'pct_in']

    ## round float values for cleaner output
    top_genes_df1[["logfoldchanges", "pct_in"]] = top_genes_df1[["logfoldchanges", "pct_in"]].round(4).astype(str).astype(float)
    
    for cluster, group_df in top_genes_df1.groupby('group', observed=True):
        cluster_dict = {}
        for record in group_df[cols_to_keep].to_dict(orient='records'):
            gene_symbol = record.pop('names')
            cluster_dict[gene_symbol] = record
        final_gene_dict1[str(cluster)] = cluster_dict

    # Output 2 : based on pct_in - pct_out
    final_gene_dict2 = {}
    cols_to_keep = ['names', 'logfoldchanges', 'pvals_adj', 'pct_diff']

    ## round float values for cleaner output
    top_genes_df2[["logfoldchanges", "pct_diff"]] = top_genes_df2[["logfoldchanges", "pct_diff"]].round(4).astype(str).astype(float)
    print(top_genes_df2[cols_to_keep].head())

    for cluster, group_df in top_genes_df2.groupby('group', observed=True):
        cluster_dict = {}
        for record in group_df[cols_to_keep].to_dict(orient='records'):
            gene_symbol = record.pop('names')
            cluster_dict[gene_symbol] = record
        final_gene_dict2[str(cluster)] = cluster_dict

    logger.info(f"Extracted markers for {len(final_gene_dict1)} clusters")
    return final_gene_dict1, final_gene_dict2
