#!/usr/bin/env python3
"""
LLM client wrapper for cell type annotation using Gemini models.
"""

from google import genai
from google.genai import types
import json
import re
import os
from typing import Dict, List, Optional, Union
import pandas as pd
from .logger import PipelineLogger

logger = PipelineLogger.get_logger(__name__)


def extract_json_from_response(text: str) -> str:
    """
    Extract JSON content from LLM response, handling various formats.
    
    Args:
        text: Raw text response from LLM
        
    Returns:
        Cleaned JSON string
    """
    # Remove markdown code blocks
    text = re.sub(r'```json\s*', '', text)
    text = re.sub(r'```\s*', '', text)
    
    # Try to find JSON object or array
    # Look for content between { } or [ ]
    json_match = re.search(r'(\{[\s\S]*\}|\[[\s\S]*\])', text)
    if json_match:
        return json_match.group(1).strip()
    
    # If no JSON structure found, return the original text for downstream error handling
    return text.strip()


class CellTypeAnnotationClient:
    """Wrapper for LLM-based cell type annotation."""
    
    def __init__(
        self, 
        model_name: str = "gemini-2.5-flash",
        response_mime_type: str = "application/json",
        max_retries: int = 3,
    ):
        """
        Initialize the Gemini API client.
        
        Args:
            model_name: Name of the LLM model to use
            response_mime_type: Expected MIME type of the API response
            max_retries: Number of retries for JSON parsing errors
        """
        self.client = genai.Client()
        self.model_name = model_name
        self.response_mime_type = response_mime_type
        self.max_retries = max_retries
        logger.info(f"Initialized Gemini API client with model: {model_name}")
        logger.info(f"Parameters: max_retries={max_retries}, response_mime_type={response_mime_type}")

    def select_cluster_resolution(self, prompt: str, custom_model: str = None) -> Dict:
        """
        Queries the Gemini model to determine the optimal clustering resolution 
        from a cluster tree image.
        
        Args:
            prompt: The formatted prompt for cluster resolution selection.
            custom_model: Optional custom model name to override default.
        
        Returns:
            Dictionary containing the selected resolution and justification.
            Example: {'resolution': 0.15, 'level_color': 'red', 'justification': '...'}
            
        Raises:
            Exception: If all retry attempts fail
        """     
        model_to_use = custom_model if custom_model else self.model_name
        logger.info(f"Requesting cluster resolution selection from LLM using {model_to_use}")
        for attempt in range(self.max_retries):
            try:
                input_token = self.client.models.count_tokens(
                    model=model_to_use,
                    contents=prompt
                )
                logger.info(f"[Cluster determination] Input prompt token count: {input_token.total_tokens}")
                cls_res = self.client.models.generate_content(
                    model=model_to_use,
                    contents=prompt,
                    config={
                        "response_mime_type": self.response_mime_type
                    }
                )
                logger.info(f"[Cluster determination] Thinking tokens used: {cls_res.usage_metadata.thoughts_token_count}")
                logger.info(f"[Cluster determination] Output tokens used: {cls_res.usage_metadata.candidates_token_count}")
                logger.info(f"[Cluster determination] Total tokens used: {cls_res.usage_metadata.total_token_count}")

                # Parse JSON response
                result = json.loads(cls_res.text)
                
                # Validate required fields
                if 'resolution' not in result:
                    raise ValueError("Response missing required 'resolution' field")
                
                logger.info(f"Selected resolution: {result.get('resolution', 'N/A')}")
                return result
                
            except json.JSONDecodeError as e:
                logger.warning(f"JSON parsing error (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                logger.debug(f"Raw response: {cls_res.text if 'cls_res' in locals() else 'N/A'}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed to parse JSON after {self.max_retries} attempts")
                    raise
            except Exception as e:
                logger.warning(f"Error in cluster resolution selection (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed after {self.max_retries} attempts")
                    raise

    def annotate_cell_types(self, prompt: str, nreplicate: int = 3, custom_model: str = None, custom_model_consolidation: Optional[str] = None) -> Dict:
        """
        Queries the Gemini model to perform cell type annotation based on marker genes 
        and UMAP visualization.
        
        Args:
            prompt: The formatted prompt containing marker genes, pathways, 
                   and candidate cell types.
            nreplicate: Number of replicate annotations to request for robustness.
            custom_model: Optional custom model name to override default.
            custom_model_consolidation: Optional custom model name for consolidation step, only effective if nreplicate>1.

        Returns:
            List of dictionaries containing cell type annotations for each cluster.
            Each dict has: {'cluster_id': str, 'cell_type_hypotheses': str, 
                           'justification': str, 'key_markers_cited': list, 
                           'confidence': str}
                           
        Raises:
            Exception: If all retry attempts fail
        """
        model_to_use = custom_model if custom_model else self.model_name
        logger.info(f"Requesting cell type annotation from LLM using {model_to_use}")
        
        for attempt in range(self.max_retries):
            try:
                input_token = self.client.models.count_tokens(
                    model=model_to_use,
                    contents=prompt
                )
                logger.info(f"[Cell annotation] Input prompt token count: {input_token.total_tokens}")
                res = []
                total_tokens = 0
                for i in range(nreplicate):
                    if nreplicate>1: logger.info(f"Generating annotation replicate {i + 1}/{nreplicate}")
                    ann_res = self.client.models.generate_content(
                        model=model_to_use,
                        contents=prompt,
                        config={
                            "response_mime_type": self.response_mime_type
                        }
                    )
                    res.append(ann_res)
                    total_tokens += ann_res.usage_metadata.total_token_count
                    if nreplicate>1: logger.info(f"[Cell annotation] {total_tokens} total tokens used until now.")
                logger.info(f"[Cell annotation] Annotation done. Total tokens used for annotation process: {total_tokens}")

                ## run consolidation if multiple replicates
                if nreplicate>1:
                    custom_model_consolidation = custom_model_consolidation if custom_model_consolidation else self.model_name
                    logger.info(f"[Cell annotation - Consolidation] Consolidating {nreplicate} annotation replicates with model {custom_model_consolidation}")
                    from .prompts import CONSENSUS_INSTRUCTION

                    ## combine all annotations into pandas DataFrame
                    res_df = []
                    for i in range(len(res)):
                        r = res[i]
                        res_df.append(pd.DataFrame(json.loads(r.text)).assign(rep = i))
                    res_df = pd.concat(res_df).sort_values("cluster_id")
                    ## convert it to dict 
                    res_df_dict = {}
                    for cluster, group_df in res_df.groupby("cluster_id", observed=True):
                        res_df_dict[str(cluster)] = group_df.drop(columns="cluster_id").to_dict(orient='records')

                    ## build prompt for consolidation
                    ann_consensus_contents = [
                        CONSENSUS_INSTRUCTION,
                        json.dumps(res_df_dict)
                    ]
                    input_token = self.client.models.count_tokens(
                        model=custom_model_consolidation,
                        contents=ann_consensus_contents
                    )
                    logger.info(f"[Cell annotation - Consolidation] Input prompt token count: {input_token.total_tokens}")
                    ann_consensus_res = self.client.models.generate_content(
                        model=custom_model_consolidation,
                        contents=ann_consensus_contents,
                        config={
                            "response_mime_type": self.response_mime_type
                        }
                    )
                    logger.info(f"[Cell annotation - Consolidation] Done. Tokens used: {ann_consensus_res.usage_metadata.total_token_count}")
                    
                    # Parse JSON response
                    result = json.loads(ann_consensus_res.text)
                else:
                    # Parse JSON response
                    result = json.loads(res[0].text)
                    res_df_dict = None
                logger.info(f"Annotation cell types for {len(result)} clusters done.")
                return result, res_df_dict
                
            except json.JSONDecodeError as e:
                logger.warning(f"JSON parsing error (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                logger.debug(f"Raw response: {res[0].text if 'res' in locals() and res[0] else 'N/A'}")
                logger.debug(f"Consolidation response: {ann_consensus_res.text if 'ann_consensus_res' in locals() else 'N/A'}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed to parse JSON after {self.max_retries} attempts")
                    raise
            except Exception as e:
                logger.warning(f"Error in cell type annotation (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed after {self.max_retries} attempts")
                    raise