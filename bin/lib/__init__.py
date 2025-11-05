#!/usr/bin/env python3
"""
Library modules for cell type annotation pipeline.
"""

from .logger import PipelineLogger
from .tree_utils import find_node_data, get_immediate_children, load_tree, get_tree_depth, TreeValidationError
from .data_processing import prepare_subset_dataset, prepare_celltype_inputs
from .marker_genes import get_top_marker_genes
from .pathway_enrichment import get_top_enriched_pathways
from .cluster_analysis import get_cluster_adjacency
from .llm_client import CellTypeAnnotationClient
from .annotator import annotate_cell_types, HierarchicalAnnotation
from .prompts import cluster_PROMPT, Celltyper_Instruction, Consolidator_Instruction

__all__ = [
    'PipelineLogger',
    'find_node_data',
    'get_immediate_children',
    'load_tree',
    'get_tree_depth',
    'TreeValidationError',
    'prepare_subset_dataset',
    'prepare_celltype_inputs',
    'get_top_marker_genes',
    'get_top_enriched_pathways',
    'get_cluster_adjacency',
    'CellTypeAnnotationClient',
    'annotate_cell_types',
    'HierarchicalAnnotation',
    'cluster_PROMPT',
    'Celltyper_Instruction',
    'Consolidator_Instruction',
]
