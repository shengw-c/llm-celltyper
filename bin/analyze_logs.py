#!/usr/bin/env python3
"""
Log analysis script that uses LLM to extract insights from pipeline logs.
"""

import sys
import json
import os
import argparse
from pathlib import Path
from typing import List, Dict
from google import genai
from google.genai import types

# Add lib directory to path
sys.path.insert(0, str(Path(__file__).parent / "lib"))

from lib import PipelineLogger

# Logger will be initialized in main() with the user-specified log level
logger = None


LOG_ANALYSIS_PROMPT = """You are an expert data scientist analyzing logs from a single-cell RNA-seq cell type annotation pipeline.

Your task is to extract key insights, warnings, errors, and performance metrics from the provided log files.

## Instructions

1. Read through all the provided log content carefully
2. Identify:
   - Successfully completed annotation runs and their details
   - Any errors or failures that occurred
   - Performance metrics (cell counts, gene counts, cluster counts, etc.)
   - Warnings or notable issues
   - Quality indicators (confidence scores, marker specificity, etc.)

3. Provide a structured JSON summary with the following format:

```json
{{
  "summary": {{
    "total_annotations": <number of annotation runs>,
    "successful_runs": <number of successful completions>,
    "failed_runs": <number of failures>,
    "total_cells_processed": <total number of cells>,
    "hierarchy_levels_completed": <deepest level reached>
  }},
  "runs": [
    {{
      "nametag": "<annotation run identifier>",
      "level": <hierarchy level number>,
      "status": "success|failed",
      "cell_count": <number of cells>,
      "cluster_count": <number of clusters>,
      "resolution": <clustering resolution>,
      "cell_types_identified": [<list of identified cell types>],
      "issues": [<list of any issues encountered>]
    }}
  ],
  "errors": [
    {{
      "message": "<error message>",
      "context": "<where it occurred>"
    }}
  ],
  "warnings": [
    {{
      "message": "<warning message>",
      "context": "<context>"
    }}
  ],
  "performance_notes": [
    "<any notable performance observations>"
  ],
  "recommendations": [
    "<suggestions for improvement or areas of concern>"
  ]
}}
```

Your response must be ONLY valid JSON with no additional text, markdown, or explanations.
"""


def read_log_content(log_files: List[str], max_chars: int = 500000) -> str:
    """
    Read and combine content from multiple log files.
    
    Args:
        log_files: List of log file paths
        max_chars: Maximum characters to read (to avoid token limits)
    
    Returns:
        Combined log content
    """
    combined_content = []
    total_chars = 0
    
    for log_file in log_files:
        try:
            with open(log_file, 'r') as f:
                content = f.read()
                if total_chars + len(content) > max_chars:
                    # Truncate to fit within limit
                    remaining = max_chars - total_chars
                    content = content[:remaining]
                    combined_content.append(f"\n{'='*80}\n")
                    combined_content.append(f"LOG: {log_file}\n")
                    combined_content.append(f"{'='*80}\n")
                    combined_content.append(content)
                    logger.warning(f"Truncated logs at {max_chars} characters")
                    break
                
                combined_content.append(f"\n{'='*80}\n")
                combined_content.append(f"LOG: {log_file}\n")
                combined_content.append(f"{'='*80}\n")
                combined_content.append(content)
                total_chars += len(content)
                
        except Exception as e:
            logger.error(f"Error reading log file {log_file}: {str(e)}")
    
    return ''.join(combined_content)


def analyze_logs_with_llm(log_content: str, model: str = "gemini-2.5-flash") -> Dict:
    """
    Use LLM to analyze log content and extract insights.
    
    Args:
        log_content: Combined log file content
        model: LLM model name to use
    
    Returns:
        Dictionary containing analysis results
    """
    logger.info(f"Analyzing logs with LLM (model: {model})")
    
    try:
        client = genai.Client()
        
        response = client.models.generate_content(
            model=model,
            contents=[
                {"text": LOG_ANALYSIS_PROMPT},
                {"text": f"\n\nLOG CONTENT:\n{log_content}"}
            ],
            config=types.GenerateContentConfig(
                temperature=0.1,
            )
        )
        
        # Clean and parse JSON response
        result = json.loads(response.text.replace("```json", "").replace("```", ""))
        logger.info("Log analysis completed successfully")
        return result
        
    except Exception as e:
        logger.error(f"Error during log analysis: {str(e)}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Analyze pipeline logs using LLM to extract insights"
    )
    parser.add_argument(
        '--log-files',
        help='Comma-separated list of log files to analyze'
    )
    parser.add_argument(
        '--output',
        default='work/log_analysis.json',
        help='Output file for analysis results (default: work/log_analysis.json)'
    )
    parser.add_argument(
        '--max-chars',
        type=int,
        default=500000,
        help='Maximum characters to read from logs (default: 500000)'
    )
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)'
    )
    parser.add_argument(
        '--llm-model',
        default='gemini-2.5-flash',
        help='LLM model to use for analysis (default: gemini-2.5-flash)'
    )
    
    args = parser.parse_args()
    
    # Initialize logger with user-specified log level
    global logger
    import logging
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logger = PipelineLogger.get_logger(__name__, level=numeric_level)
    
    logger.info("=" * 80)
    logger.info("Log Analysis")
    logger.info("=" * 80)
    
    # Collect log files
    log_files = args.log_files.split(',')
    
    if not log_files:
        logger.error("No log files found to analyze")
        return 1
    
    # Read log content
    logger.info("Reading log files...")
    log_content = read_log_content(log_files, args.max_chars)
    logger.info(f"Read {len(log_content)} characters from {len(log_files)} files")
    
    # Analyze with LLM
    try:
        analysis_result = analyze_logs_with_llm(log_content, model=args.llm_model)
        
        # Save results
        with open(args.output, 'w') as f:
            json.dump(analysis_result, f, indent=2)
        
        logger.info(f"✓ Analysis saved to: {args.output}")
        
        # Print summary
        summary = analysis_result.get('summary', {})
        logger.info("\nAnalysis Summary:")
        logger.info(f"  Total annotations: {summary.get('total_annotations', 'N/A')}")
        logger.info(f"  Successful runs: {summary.get('successful_runs', 'N/A')}")
        logger.info(f"  Failed runs: {summary.get('failed_runs', 'N/A')}")
        logger.info(f"  Total cells: {summary.get('total_cells_processed', 'N/A')}")
        
        return 0
        
    except Exception as e:
        logger.error(f"✗ Log analysis failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
