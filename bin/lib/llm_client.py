#!/usr/bin/env python3
"""
LLM client wrapper for cell type annotation using Gemini models.
Supports mock mode for testing without API calls.
"""

from google import genai
import json
import os
from typing import Dict, List, Optional, Union
from .logger import PipelineLogger

from .prompts import (
    Celltyper_Instruction,
    harmonize_instruction
)
from .schemas import (
    CellTypeHypotheses,
    CellTypeHarmonizations
)

logger = PipelineLogger.get_logger(__name__)


def _is_mock_mode() -> bool:
    """Check if mock LLM mode is enabled via environment variable."""
    return os.environ.get('MOCK_LLM', '').lower() in ('1', 'true', 'yes')


def _generate_mock_annotation(cluster_data: List[Dict]) -> List[Dict]:
    """Generate mock annotation responses for testing."""
    logger.info("🎭 MOCK MODE: Generating fake annotation responses")
    mock_responses = []
    
    for i, cluster in enumerate(cluster_data):
        cluster_id = cluster.get('cluster', 'unknown')
        # For mock mode, mark every 3rd cluster for splitting
        should_split = (i % 3 == 0)
        
        mock_responses.append({
            'cluster': cluster_id,
            'cell_type_hypothesis': f'CellType_{cluster_id}',
            'cell_type_description': f'Mock description for cluster {cluster_id}',
            'confidence': 'high',
            'justification': f'Mock annotation for cluster {cluster_id}',
            'split_status': 'Yes' if should_split else 'No',
            'next_round_context': f'Subdivide {cluster_id}' if should_split else '',
            'key_markers_cited': [],
            'key_pathways_cited': []
        })
    
    return mock_responses


def _generate_mock_harmonization(cell_type_annotations: List[Dict]) -> List[Dict]:
    """Generate mock harmonization responses for testing."""
    logger.info("🎭 MOCK MODE: Generating fake harmonization responses")
    mock_responses = []
    
    for annotation in cell_type_annotations:
        cluster_id = annotation.get('cluster', 'unknown')
        label = annotation.get('label', f'CellType_{cluster_id}')
        mock_responses.append({
            'cluster': cluster_id,  # Use cluster as per schema
            'harmonized_cell_type': label,  # Keep same in mock mode
            'harmonization_log': 'No change.'
        })
    
    return mock_responses

def llm_harmonize_cell_types(
    cell_type_annotations: List[Dict],
    custom_model: Optional[str] = None,
    max_retries: int = 3,
    response_mime_type: str = "application/json"
) -> List[Dict]:
    """
    Harmonizes cell type labels across clusters using Gemini LLM.
    
    Args:
        cell_type_annotations: List of annotation dictionaries for each cluster.
        custom_model: LLM model name (default: None).
        max_retries: Maximum retry attempts on failure.
        response_mime_type: Response format (default: "application/json").

    Returns:
        List of harmonization dictionaries with standardized labels.

    Raises:
        Exception: If all retry attempts fail.
    """
    # Check for mock mode
    if _is_mock_mode():
        return _generate_mock_harmonization(cell_type_annotations)
    
    logger.info(f"Requesting cell type harmonization from LLM using {custom_model}")
    client = genai.Client()

    for attempt in range(max_retries):
        try:
            prompt = harmonize_instruction.format(
                annotations_json=json.dumps(cell_type_annotations, indent=2)
            )

            logger.info(f"Generating harmonization with `{custom_model}` model...")
            harm_res = client.models.generate_content(
                model=custom_model,
                contents=prompt,
                config={
                    "response_mime_type": response_mime_type,
                    "response_json_schema": CellTypeHarmonizations.model_json_schema()
                }
            )
            logger.info(f"[Cell type harmonization] Harmonization done. ")
            logger.info(f"[Cell type harmonization done] Input tokens: {harm_res.usage_metadata.prompt_token_count}")
            logger.info(f"[Cell type harmonization done] Thoughts tokens: {harm_res.usage_metadata.thoughts_token_count}")
            logger.info(f"[Cell type harmonization done] Tokens used in this iteration: {harm_res.usage_metadata.total_token_count}")

            ## load results into json
            harm_res_json = json.loads(harm_res.text)
            return harm_res_json.get('harmonizations', [])
            
        except Exception as e:
            logger.warning(f"Error in cell type harmonization (attempt {attempt + 1}/{max_retries}): {str(e)}")
            if attempt == max_retries - 1:
                logger.error(f"Failed after {max_retries} attempts")
                raise
            
def llm_cell_typer(
    cluster_data: List[Dict],
    custom_model: Optional[str] = None,
    max_retries: int = 3,
    response_mime_type: str = "application/json"
) -> List[Dict]:
    """
    Performs cell type annotation using Gemini LLM based on marker genes and pathways.
    
    Args:
        cluster_data: List of cluster dictionaries with marker genes and pathways.
        custom_model: LLM model name (default: None).
        max_retries: Maximum retry attempts on failure.
        response_mime_type: Response format (default: "application/json").

    Returns:
        List of annotation dictionaries with cell type hypotheses, confidence,
        split status, and supporting evidence.

    Raises:
        Exception: If all retry attempts fail.
    """
    # Check for mock mode
    if _is_mock_mode():
        return _generate_mock_annotation(cluster_data)
    
    logger.info(f"Requesting cell type annotation from LLM using {custom_model}")
    client = genai.Client()

    for attempt in range(max_retries):
        try:
            logger.info(f"Generating annotation with `{custom_model}` model...")
            ## If there are too many sub-clusters, process in batches to avoid timeouts
            chunk_size = 50
            sub_cluster_batches = [cluster_data[i:i + chunk_size] for i in range(0, len(cluster_data), chunk_size)]
            if len(sub_cluster_batches) > 1:
                logger.info(f"Processing {len(cluster_data)} sub-clusters in {len(sub_cluster_batches)} batches of size {chunk_size}")
            ann_res = []
            total_token_count = 0
            for i in range(len(sub_cluster_batches)):
                logger.info(f"Processing batch {i+1}/{len(sub_cluster_batches)} with {len(sub_cluster_batches[i])} sub-clusters")
                batch = sub_cluster_batches[i]
                prompt = Celltyper_Instruction.format(
                    cluster_data_json=json.dumps(batch)
                )
                sub_ann_res = client.models.generate_content(
                    model=custom_model,
                    contents=prompt,
                    config={
                        "response_mime_type": response_mime_type,
                        "response_json_schema": CellTypeHypotheses.model_json_schema()
                    }
                )
                sub_ann_res_json = json.loads(sub_ann_res.text)
                ann_res.extend(sub_ann_res_json.get('hypotheses', []))
            
                logger.info(f"[Cell annotation {i+1}/{len(sub_cluster_batches)}] Annotation done. ")
                logger.info(f"[Cell annotation {i+1}/{len(sub_cluster_batches)}] Input tokens: {sub_ann_res.usage_metadata.prompt_token_count}")
                logger.info(f"[Cell annotation {i+1}/{len(sub_cluster_batches)}] Thoughts tokens: {sub_ann_res.usage_metadata.thoughts_token_count}")
                logger.info(f"[Cell annotation {i+1}/{len(sub_cluster_batches)}] Tokens used in this iteration: {sub_ann_res.usage_metadata.total_token_count}")
                total_token_count += sub_ann_res.usage_metadata.total_token_count
            
            logger.info(f"[Cell annotation done] Total tokens used across all batches: {total_token_count}")
            
            if len(sub_cluster_batches)>1:
                logger.info(f"Harmonizing cell type labels after batch annotation from {len(sub_cluster_batches)} batches...")
                ## Do harmonization after annotation
                to_harmonize = []
                for ann in ann_res:
                    to_harmonize.append({
                        "cluster": ann["cluster"],
                        "label": ann["cell_type_hypothesis"],
                        "definition": ann["cell_type_description"]
                    })
                
                logger.info("Harmonizing cell type labels after annotation...")
                harmonization_res = llm_harmonize_cell_types(
                    cell_type_annotations=to_harmonize,
                    custom_model=custom_model,
                    max_retries=max_retries
                )
                # Merge harmonization results back into annotation results
                harmonization_dict = {h['cluster']: h["harmonized_cell_type"] for h in harmonization_res}
                for ann in ann_res:
                    ann.update({
                        "cell_type_hypothesis": harmonization_dict.get(ann["cluster"], ann["cell_type_hypothesis"])
                    })

            return ann_res
            
        except Exception as e:
            logger.warning(f"Error in cell type annotation (attempt {attempt + 1}/{max_retries}): {str(e)}")
            if attempt == max_retries - 1:
                logger.error(f"Failed after {max_retries} attempts")
                raise

