#!/usr/bin/env python3
"""
Wrapper script for running cell type annotation in Nextflow pipeline.
This script wraps the existing annotate_cell_types function from lib/annotator.py
"""

import sys
import argparse
import os
import logging
from pathlib import Path
from typing import Optional

# Add lib directory to path
sys.path.insert(0, str(Path(__file__).parent / "lib"))

from lib.annotator import annotate_cell_types
from lib.logger import PipelineLogger

# Module-level logger
logger: Optional[logging.Logger] = None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run hierarchical cell type annotation"
    )
    
    # Required arguments
    parser.add_argument('--input', required=True, help='Path to QCed h5ad file')
    parser.add_argument('--context', required=True, help='Biological context description')
    
    # Optional arguments
    parser.add_argument('--batch-key', default=None, help='Batch key for integration')
    parser.add_argument('--integration', action='store_true', help='Enable integration')
    parser.add_argument('--cpus', type=int, default=16, help='Number of CPUs (default: 16)')
    parser.add_argument('--output-dir', default='work', help='Output directory (default: work)')
    parser.add_argument('--min-cells', type=int, default=1000, 
                       help='Minimum cells for subtype annotation (default: 1000)')
    parser.add_argument('--max-resolution', type=float, default=1.0, 
                       help='Maximum resolution for Leiden clustering (default: 1.0)')
    parser.add_argument('--transition-cutoff', type=float, default=0.1,
                       help='Stability threshold for clustering (default: 0.1)')
    parser.add_argument('--llm-model', default='gemini-2.5-flash',
                       help='LLM model name for annotation (default: gemini-2.5-flash)')
    parser.add_argument('--llm-max-retries', type=int, default=3,
                       help='LLM max retries (default: 3)')
    parser.add_argument('--gsea-databases', default='MSigDB_Hallmark_2020,KEGG_2021_Human',
                       help='Comma-separated list of gene set databases')
    parser.add_argument('--top-genes', type=int, default=20,
                       help='Number of top marker genes per cluster (default: 20)')
    parser.add_argument('--top-pathways', type=int, default=20,
                       help='Number of top pathways per cluster (default: 20)')
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level (default: INFO)')
    
    args = parser.parse_args()
    
    # Initialize logger
    global logger
    log_dir = os.path.join(args.output_dir, "logs")
    
    import logging
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logger = PipelineLogger.get_logger(__name__, log_dir=log_dir, level=numeric_level)
    
    # Convert batch_key string to list if provided
    batch_keys = [args.batch_key] if args.batch_key and args.batch_key != "None" else None
    
    logger.info("=" * 80)
    logger.info("Starting hierarchical cell type annotation")
    logger.info("=" * 80)
    logger.info(f"Configuration:")
    logger.info(f"  Input: {args.input}")
    logger.info(f"  Context: {args.context}")
    logger.info(f"  Batch keys: {batch_keys}")
    logger.info(f"  Integration: {args.integration}")
    logger.info(f"  CPUs: {args.cpus}")
    logger.info(f"  Output dir: {args.output_dir}")
    logger.info(f"  Min cells for subtype: {args.min_cells}")
    logger.info(f"  Max resolution: {args.max_resolution}")
    logger.info(f"  Transition cutoff: {args.transition_cutoff}")
    logger.info(f"  LLM model: {args.llm_model}")
    logger.info(f"  LLM max retries: {args.llm_max_retries}")
    logger.info(f"  GSEA databases: {args.gsea_databases}")
    logger.info(f"  Top genes: {args.top_genes}")
    logger.info(f"  Top pathways: {args.top_pathways}")
    
    # Check for mock mode environment variable
    if os.environ.get('MOCK_LLM_MODE') == '1':
        logger.info("  ⚠️  MOCK LLM MODE ENABLED - No actual API calls will be made")
        os.environ['MOCK_LLM'] = '1'  # Set environment variable for llm_client to use
    
    # Run annotation
    try:
        result = annotate_cell_types(
            general_context=args.context,
            input_file=args.input,
            batch_key=batch_keys,
            integration=args.integration,
            cpus_per_task=args.cpus,
            min_cells_for_subtype=args.min_cells,
            max_resolution=args.max_resolution,
            transition_cutoff=args.transition_cutoff,
            llm_model=args.llm_model,
            llm_max_retries=args.llm_max_retries,
            gsea_databases=args.gsea_databases,
            top_genes=args.top_genes,
            top_pathways=args.top_pathways,
            output_dir=args.output_dir
        )
        
        logger.info("✓ Annotation completed successfully")
        return 0
        
    except Exception as e:
        logger.error(f"✗ Error during annotation: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
