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
    
    if cell_type_name not in cell_tree:
        raise ValueError(f"Cell type '{cell_type_name}' not found in tree file")
    
    cell_data = cell_tree[cell_type_name]
    
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
    llm_model: str = "gemini-2.0-flash-exp",
    llm_temperature: float = 0.1,
    llm_max_retries: int = 3,
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
        llm_model: Name of the LLM model to use.
        llm_temperature: Temperature for LLM generation.
        llm_max_retries: Maximum retries for LLM API calls.
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
    logger.info(f"  LLM model: {llm_model}")
    logger.info(f"  LLM temperature: {llm_temperature}")
    logger.info(f"  LLM max retries: {llm_max_retries}")
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
            llm_model=llm_model,
            llm_temperature=llm_temperature,
            llm_max_retries=llm_max_retries,
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
    
    # Subcommand to split dataset by all top-level cell types
    split_parser = subparsers.add_parser('split', help='Split dataset by all top-level cell types')
    split_parser.add_argument('--input', required=True, help='Path to QCed h5ad file')
    split_parser.add_argument('--tree', required=True, help='Path to cell type hierarchy JSON file')
    split_parser.add_argument('--context', required=True, help='Biological context description')
    split_parser.add_argument('--batch-key', default=None, help='Batch key for integration')
    split_parser.add_argument('--integration', action='store_true', help='Enable integration')
    split_parser.add_argument('--cpus', type=int, default=16, help='Number of CPUs (default: 16)')
    split_parser.add_argument('--output-dir', default='work', help='Output directory (default: work)')
    split_parser.add_argument('--min-cells', type=int, default=1000, 
                             help='Minimum cells for subtype annotation (default: 1000)')
    split_parser.add_argument('--condition-num', type=int, default=1,
                             help='Number of expected condition-driven groups or artifacts (default: 1)')
    
    # LLM model configuration for split
    split_parser.add_argument('--llm-model', default='gemini-2.0-flash-exp',
                             help='LLM model name (default: gemini-2.0-flash-exp)')
    split_parser.add_argument('--llm-temperature', type=float, default=0.1,
                             help='LLM temperature (default: 0.1)')
    split_parser.add_argument('--llm-max-retries', type=int, default=3,
                             help='LLM max retries (default: 3)')
    split_parser.add_argument('--gsea-databases', default='MSigDB_Hallmark_2020,KEGG_2021_Human',
                             help='Comma-separated list of gene set databases')
    split_parser.add_argument('--top-genes', type=int, default=20,
                             help='Number of top marker genes per cluster (default: 20)')
    split_parser.add_argument('--top-pathways', type=int, default=20,
                             help='Number of top pathways per cluster (default: 20)')
    split_parser.add_argument('--log-level', default='INFO',
                             choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                             help='Logging level (default: INFO)')
    
    # Subcommand to run annotation
    annotate_parser = subparsers.add_parser('annotate', help='Run annotation for a cell type')
    annotate_parser.add_argument('--input', required=True, help='Path to QCed h5ad file')
    annotate_parser.add_argument('--tree', required=True, help='Path to cell type hierarchy JSON file')
    annotate_parser.add_argument('--celltype', required=True, help='Major cell type to annotate')
    annotate_parser.add_argument('--context', required=True, help='Biological context description')
    annotate_parser.add_argument('--batch-key', default=None, help='Batch key for integration')
    annotate_parser.add_argument('--integration', action='store_true', help='Enable integration')
    annotate_parser.add_argument('--cpus', type=int, default=16, help='Number of CPUs (default: 16)')
    annotate_parser.add_argument('--output-dir', default='work', help='Output directory (default: work)')
    annotate_parser.add_argument('--min-cells', type=int, default=1000, 
                                 help='Minimum cells for subtype annotation (default: 1000)')
    annotate_parser.add_argument('--condition-num', type=int, default=1,
                                 help='Number of expected condition-driven groups or artifacts (default: 1)')
    
    # LLM model configuration
    annotate_parser.add_argument('--llm-model', default='gemini-2.0-flash-exp',
                                 help='LLM model name (default: gemini-2.0-flash-exp)')
    annotate_parser.add_argument('--llm-temperature', type=float, default=0.1,
                                 help='LLM temperature (default: 0.1)')
    annotate_parser.add_argument('--llm-max-retries', type=int, default=3,
                                 help='LLM max retries (default: 3)')
    
    # GSEA configuration
    annotate_parser.add_argument('--gsea-databases', default='MSigDB_Hallmark_2020,KEGG_2021_Human',
                                 help='Comma-separated list of gene set databases (default: MSigDB_Hallmark_2020,KEGG_2021_Human)')
    annotate_parser.add_argument('--top-genes', type=int, default=20,
                                 help='Number of top marker genes per cluster (default: 20)')
    annotate_parser.add_argument('--top-pathways', type=int, default=20,
                                 help='Number of top pathways per cluster (default: 20)')
    
    # Annotation control
    annotate_parser.add_argument('--no-nest', action='store_false', dest='nested',
                                 help='Disable nested/recursive annotation of subtypes (only annotate top-level)')
    annotate_parser.add_argument('--log-level', default='INFO',
                                 choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                                 help='Logging level (default: INFO)')
    
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
            llm_model=args.llm_model,
            llm_temperature=args.llm_temperature,
            llm_max_retries=args.llm_max_retries,
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
            llm_model=args.llm_model,
            llm_temperature=args.llm_temperature,
            llm_max_retries=args.llm_max_retries,
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
