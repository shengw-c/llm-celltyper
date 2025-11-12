#!/usr/bin/env python3
"""
Wrapper script for running cell type annotation in Nextflow pipeline.
This script extracts a single major cell type and runs hierarchical annotation.
"""

import sys
import json
import argparse
import os
from pathlib import Path

# Add lib directory to path
sys.path.insert(0, str(Path(__file__).parent / "lib"))

from lib import (
    PipelineLogger,
    annotate_cell_types,
    HierarchicalAnnotation
)

# Module-level logger (will be reinitialized in run_annotation with cell_type)
logger = None


def get_major_cell_types(tree_file: str):
    """
    Extract all major (top-level) cell types from the hierarchy tree.
    
    Args:
        tree_file: Path to the cell type hierarchy JSON file.
    
    Returns:
        List of major cell type names.
    """
    with open(tree_file, 'r') as f:
        cell_tree = json.load(f)
    
    return list(cell_tree.keys())


def prepare_cell_type_config(tree_file: str, cell_type_name: str):
    """
    Prepare the expected_cells configuration for a specific major cell type.
    Returns the immediate children of the specified cell type for annotation.
    
    Args:
        tree_file: Path to the cell type hierarchy JSON file.
        cell_type_name: Name of the major cell type to process.
    
    Returns:
        Dictionary: Expected cells configuration containing the children of the specified cell type.
                   Format: {child_name: {'definition': str, 'has_children': bool}}
    """
    with open(tree_file, 'r') as f:
        cell_tree = json.load(f)
    
    # Normalize cell type name for lookup - try original, then with underscores, then with spaces
    # This handles different JSON formats (some use spaces, some use underscores)
    lookup_key = None
    if cell_type_name in cell_tree:
        lookup_key = cell_type_name
    elif cell_type_name.replace(' ', '_') in cell_tree:
        lookup_key = cell_type_name.replace(' ', '_')
    elif cell_type_name.replace('_', ' ') in cell_tree:
        lookup_key = cell_type_name.replace('_', ' ')
    
    if lookup_key is None:
        available_keys = ', '.join(list(cell_tree.keys())[:5])
        raise ValueError(
            f"Cell type '{cell_type_name}' not found in tree file. "
            f"Available keys (first 5): {available_keys}..."
        )
    
    cell_data = cell_tree[lookup_key]
    
    # Import get_immediate_children from lib
    from lib.tree_utils import get_immediate_children
    
    # Return the immediate children of this cell type
    return get_immediate_children(cell_data)


def prepare_all_cell_types_config(tree_file: str):
    """
    Prepare the expected_cells configuration for ALL major cell types.
    
    Args:
        tree_file: Path to the cell type hierarchy JSON file.
    
    Returns:
        Dictionary: Expected cells configuration for all top-level cell types.
    """
    with open(tree_file, 'r') as f:
        cell_tree = json.load(f)
    
    expected_cells = {}
    for cell_type_name, cell_data in cell_tree.items():
        expected_cells[cell_type_name] = {
            "has_children": bool(cell_data.get("children", {})),
            "definition": cell_data.get("definition", "")
        }
    
    return expected_cells


def run_annotation(
    input_file: str,
    tree_file: str,
    general_context: str,
    cell_type: str = None,
    all_types: bool = False,
    batch_key: str = None,
    integration: bool = False,
    cpus_per_task: int = 16,
    output_dir: str = ".",
            min_cells_for_subtype: int = 1000,
            condition_num: int = 1,
            max_resolution: float = 1.0,
            llm_model_general: str = "gemini-2.5-flash",    llm_model_complicated: str = "gemini-2.5.pro",
    llm_max_retries: int = 3,
    llm_nreplies: int = 1,
    llm_adaptive: bool = True,
    gsea_databases: str = "MSigDB_Hallmark_2020,KEGG_2021_Human",
    top_genes: int = 20,
    top_pathways: int = 20,
    nested: bool = True,
    log_level: str = "INFO"
):
    """
    Run hierarchical cell type annotation for one or all major cell types.
    
    Args:
        input_file: Path to QCed h5ad file.
        tree_file: Path to cell type hierarchy JSON file.
        general_context: Biological context description.
        cell_type: Specific cell type to annotate (required if all_types=False).
        all_types: If True, annotate ALL top-level cell types at once (for split mode).
        batch_key: Batch key for integration.
        integration: Whether to perform integration.
        cpus_per_task: Number of CPUs to use.
        output_dir: Output directory for results.
        min_cells_for_subtype: Minimum cells required for subtype annotation.
        condition_num: Number of expected condition-driven groups or artifacts.
        llm_model_general: Name of the LLM model for general purpose (default: "gemini-2.5-flash").
        llm_model_complicated: Name of the LLM model for complicated tasks (default: "gemini-2.5.pro").
        llm_nreplies: Number of replies to request from LLM (default: 1, used for cell type annotation).
        llm_max_retries: Maximum retries for LLM API calls (default: 3).
        llm_adaptive: If True, only replicate annotations for Low/Medium confidence clusters (default: True).
        gsea_databases: Comma-separated list of gene set databases.
        top_genes: Number of top marker genes per cluster.
        top_pathways: Number of top pathways per cluster.
        nested: Whether to perform recursive subtype annotation.
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR).
    """
    # Initialize logger with cell_type information
    global logger
    log_dir = os.path.join(output_dir, "logs")
    
    # Convert log level string to logging constant
    import logging
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    
    if all_types:
        logger = PipelineLogger.get_logger(__name__, log_dir=log_dir, cell_type="all_types", level=numeric_level)
    else:
        logger = PipelineLogger.get_logger(__name__, log_dir=log_dir, cell_type=cell_type, level=numeric_level)
    
    # Prepare cell type configuration based on mode
    if all_types:
        logger.info("=" * 80)
        logger.info("Starting annotation for ALL top-level cell types (SPLIT mode)")
        logger.info("=" * 80)
        expected_cells = prepare_all_cell_types_config(tree_file)
        nametag = "lev0_all_types"
        mode_description = f"all {len(expected_cells)} cell types"
    else:
        if not cell_type:
            raise ValueError("cell_type is required when all_types=False")
        logger.info("=" * 80)
        logger.info(f"Starting annotation for: {cell_type}")
        logger.info("=" * 80)
        expected_cells = prepare_cell_type_config(tree_file, cell_type)
        nametag = f"lev0_{cell_type.replace(' ', '_')}"
        mode_description = cell_type
    
    # Convert batch_key string to list if provided
    batch_keys = [batch_key] if batch_key and batch_key != "None" else None
    
    logger.info(f"Configuration:")
    logger.info(f"  Input: {input_file}")
    logger.info(f"  Mode: {mode_description}")
    if all_types:
        logger.info(f"  Cell types: {list(expected_cells.keys())}")
    logger.info(f"  Nametag: {nametag}")
    logger.info(f"  Batch keys: {batch_keys}")
    logger.info(f"  Integration: {integration}")
    logger.info(f"  CPUs: {cpus_per_task}")
    logger.info(f"  Output dir: {output_dir}")
    logger.info(f"  Min cells for subtype: {min_cells_for_subtype}")
    logger.info(f"  Condition num: {condition_num}")
    logger.info(f"  LLM model (general): {llm_model_general}")
    logger.info(f"  LLM model (complicated): {llm_model_complicated}")
    logger.info(f"  LLM max retries: {llm_max_retries}")
    logger.info(f"  LLM nreplies: {llm_nreplies}")
    logger.info(f"  LLM adaptive mode: {llm_adaptive}")
    logger.info(f"  GSEA databases: {gsea_databases}")
    logger.info(f"  Top genes: {top_genes}")
    logger.info(f"  Top pathways: {top_pathways}")
    logger.info(f"  Nested annotation: {nested}")
    
    # Run annotation
    try:
        collector = annotate_cell_types(
            expected_cells=expected_cells,
            tree_file=tree_file,
            general_context=general_context,
            nametag=nametag,
            input_file=input_file,
            batch_key=batch_keys,
            integration=integration,
            subset_cell_id_file=None,
            cpus_per_task=cpus_per_task,
            current_level=0,
            hierarchical_collector=None,
            min_cells_for_subtype=min_cells_for_subtype,
            condition_num=condition_num,
            max_resolution=max_resolution,
            llm_model_general=llm_model_general,
            llm_model_complicated=llm_model_complicated,
            llm_max_retries=llm_max_retries,
            llm_nreplies=llm_nreplies,
            llm_adaptive=llm_adaptive,
            gsea_databases=gsea_databases,
            top_genes=top_genes,
            top_pathways=top_pathways,
            nested=nested,
            output_dir=output_dir
        )
        
        # Export final hierarchical annotations (only for non-split mode)
        if not all_types or nested:
            final_ann_file = collector.export_final_annotations(
                output_file=f"{output_dir}/data/{nametag}_final_annotations.tsv"
            )
            logger.info(f"Exported final annotations: {final_ann_file}")
        
        if all_types and not nested:
            logger.info(f"✓ Split completed successfully")
            logger.info(f"  Split files are in: {output_dir}/split/")
        else:
            logger.info(f"✓ Annotation completed successfully for {mode_description}")
        return 0
        
    except Exception as e:
        logger.error(f"✗ Error during annotation for {mode_description}: {str(e)}")
        import traceback
        traceback.print_exc()
        logger.error(traceback.format_exc())
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Run hierarchical cell type annotation for a specific major cell type"
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Subcommand to list major cell types
    list_parser = subparsers.add_parser('list', help='List all major cell types')
    list_parser.add_argument('--tree', required=True, help='Path to cell type hierarchy JSON file')
    list_parser.add_argument('--nextflow', action='store_true', help='Format output for Nextflow')
    
    # Common arguments shared between split and annotate subcommands
    common_args = [
        ('--input', {'required': True, 'help': 'Path to QCed h5ad file'}),
        ('--tree', {'required': True, 'help': 'Path to cell type hierarchy JSON file'}),
        ('--context', {'required': True, 'help': 'Biological context description'}),
        ('--batch-key', {'default': None, 'help': 'Batch key for integration'}),
        ('--integration', {'action': 'store_true', 'help': 'Enable integration'}),
        ('--cpus', {'type': int, 'default': 16, 'help': 'Number of CPUs (default: 16)'}),
        ('--output-dir', {'default': 'work', 'help': 'Output directory (default: work)'}),
        ('--min-cells', {'type': int, 'default': 1000, 'help': 'Minimum cells for subtype annotation (default: 1000)'}),
        ('--condition-num', {'type': int, 'default': 1, 'help': 'Number of expected condition-driven groups or artifacts (default: 1)'}),
        ('--max-resolution', {'type': float, 'default': 1.0, 'help': 'Maximum resolution for Leiden clustering.'}),
        ('--llm-model-general', {'default': 'gemini-2.5-flash', 'help': 'LLM model name for general purpose (default: gemini-2.5-flash)'}),
        ('--llm-model-complicated', {'default': 'gemini-2.5.pro', 'help': 'LLM model name for complicated tasks (default: gemini-2.5.pro)'}),
        ('--llm-max-retries', {'type': int, 'default': 3, 'help': 'LLM max retries (default: 3)'}),
        ('--llm-nreplies', {'type': int, 'default': 1, 'help': 'Number of replies to request from LLM (default: 1)'}),
        ('--llm-adaptive', {'action': 'store_true', 'default': True, 'help': 'Enable adaptive replication (only replicate low/medium confidence clusters)'}),
        ('--no-llm-adaptive', {'action': 'store_false', 'dest': 'llm_adaptive', 'help': 'Disable adaptive replication (replicate all clusters)'}),
        ('--gsea-databases', {'default': 'MSigDB_Hallmark_2020,KEGG_2021_Human', 'help': 'Comma-separated list of gene set databases'}),
        ('--top-genes', {'type': int, 'default': 20, 'help': 'Number of top marker genes per cluster (default: 20)'}),
        ('--top-pathways', {'type': int, 'default': 20, 'help': 'Number of top pathways per cluster (default: 20)'}),
        ('--log-level', {'default': 'INFO', 'choices': ['DEBUG', 'INFO', 'WARNING', 'ERROR'], 'help': 'Logging level (default: INFO)'}),
    ]
    
    # Subcommand to split dataset by all top-level cell types
    split_parser = subparsers.add_parser('split', help='Split dataset by all top-level cell types')
    for arg_name, arg_params in common_args:
        split_parser.add_argument(arg_name, **arg_params)
    
    # Subcommand to run annotation
    annotate_parser = subparsers.add_parser('annotate', help='Run annotation for a cell type')
    annotate_parser.add_argument('--celltype', required=True, help='Major cell type to annotate')
    for arg_name, arg_params in common_args:
        annotate_parser.add_argument(arg_name, **arg_params)
    annotate_parser.add_argument('--no-nest', action='store_false', dest='nested',
                                 help='Disable nested/recursive annotation of subtypes (only annotate top-level)')
    
    args = parser.parse_args()
    
    if args.command == 'list':
        # List major cell types
        cell_types = get_major_cell_types(args.tree)
        if args.nextflow:
            for ct in cell_types:
                print(ct)
        else:
            print("Major cell types found:")
            for ct in cell_types:
                print(f"  - {ct}")
            print(f"\nTotal: {len(cell_types)} major cell types")
    
    elif args.command == 'split':
        # Run split mode - annotate all top-level types and split dataset
        return run_annotation(
            input_file=args.input,
            tree_file=args.tree,
            general_context=args.context,
            all_types=True,
            nested=False,
            batch_key=args.batch_key,
            integration=args.integration,
            cpus_per_task=args.cpus,
            output_dir=args.output_dir,
            min_cells_for_subtype=args.min_cells,
            condition_num=args.condition_num,
            max_resolution=args.max_resolution,
            llm_model_general=args.llm_model_general,
            llm_model_complicated=args.llm_model_complicated,
            llm_max_retries=args.llm_max_retries,
            llm_nreplies=args.llm_nreplies,
            llm_adaptive=args.llm_adaptive,
            gsea_databases=args.gsea_databases,
            top_genes=args.top_genes,
            top_pathways=args.top_pathways,
            log_level=args.log_level
        )
        
    elif args.command == 'annotate':
        # Run annotation
        return run_annotation(
            input_file=args.input,
            tree_file=args.tree,
            general_context=args.context,
            cell_type=args.celltype,
            all_types=False,
            batch_key=args.batch_key,
            integration=args.integration,
            cpus_per_task=args.cpus,
            output_dir=args.output_dir,
            min_cells_for_subtype=args.min_cells,
            condition_num=args.condition_num,
            max_resolution=args.max_resolution,
            llm_model_general=args.llm_model_general,
            llm_model_complicated=args.llm_model_complicated,
            llm_max_retries=args.llm_max_retries,
            llm_nreplies=args.llm_nreplies,
            llm_adaptive=args.llm_adaptive,
            gsea_databases=args.gsea_databases,
            top_genes=args.top_genes,
            top_pathways=args.top_pathways,
            nested=args.nested,
            log_level=args.log_level
        )
    else:
        parser.print_help()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
