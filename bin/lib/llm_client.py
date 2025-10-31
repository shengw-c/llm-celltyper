#!/usr/bin/env python3
"""
LLM client wrapper for cell type annotation using Gemini models.
"""

from google import genai
from google.genai import types
import json
import re
from typing import Dict, List, Union

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
        model_name: str = "gemini-2.0-flash-exp",
        temperature: float = 0.1,
        max_retries: int = 3
    ):
        """
        Initialize the Gemini API client.
        
        Args:
            model_name: Name of the LLM model to use
            temperature: Temperature for response generation (0.0-1.0)
            max_retries: Maximum number of retries for failed API calls
        """
        self.client = genai.Client()
        self.model_name = model_name
        self.temperature = temperature
        self.max_retries = max_retries
        logger.info(f"Initialized Gemini API client with model: {model_name}")
        logger.info(f"Parameters: temperature={temperature}, max_retries={max_retries}")

    def select_cluster_resolution(self, prompt: str, image_data: str) -> Dict:
        """
        Queries the Gemini model to determine the optimal clustering resolution 
        from a cluster tree image.
        
        Args:
            prompt: The formatted prompt for cluster resolution selection.
            image_data: Base64-encoded PNG image of the cluster tree.
        
        Returns:
            Dictionary containing the selected resolution and justification.
            Example: {'resolution': 0.15, 'level_color': 'red', 'justification': '...'}
            
        Raises:
            Exception: If all retry attempts fail
        """
        logger.info("Requesting cluster resolution selection from LLM")
        
        for attempt in range(self.max_retries):
            try:
                response = self.client.models.generate_content(
                    model=self.model_name,
                    contents=[
                        {"text": prompt},
                        {
                            "inline_data": {
                                "mime_type": "image/png",
                                "data": image_data
                            }
                        }
                    ],
                    config=types.GenerateContentConfig(
                        temperature=self.temperature,
                    )
                )
                
                # Clean and parse JSON response
                json_text = extract_json_from_response(response.text)
                result = json.loads(json_text)
                
                # Validate required fields
                if 'resolution' not in result:
                    raise ValueError("Response missing required 'resolution' field")
                
                logger.info(f"Selected resolution: {result.get('resolution', 'N/A')}")
                return result
                
            except json.JSONDecodeError as e:
                logger.warning(f"JSON parsing error (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                logger.debug(f"Raw response: {response.text if 'response' in locals() else 'N/A'}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed to parse JSON after {self.max_retries} attempts")
                    raise
            except Exception as e:
                logger.warning(f"Error in cluster resolution selection (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed after {self.max_retries} attempts")
                    raise
    
    def annotate_cell_types(self, prompt: str, image_data: str, candidate_cell_types: Dict = None) -> Union[List[Dict], Dict]:
        """
        Queries the Gemini model to perform cell type annotation based on marker genes 
        and UMAP visualization.
        
        Args:
            prompt: The formatted prompt containing marker genes, pathways, 
                   and candidate cell types.
            image_data: Base64-encoded PNG image of the UMAP plot.
            candidate_cell_types: Optional dictionary of candidate cell types to validate against.
                                If provided, will check that returned cell types match candidates.
        
        Returns:
            List of dictionaries containing cell type annotations for each cluster.
            Each dict has: {'cluster_id': str, 'cell_type_hypotheses': str, 
                           'justification': str, 'key_markers_cited': list, 
                           'confidence': str}
                           
        Raises:
            Exception: If all retry attempts fail
        """
        logger.info("Requesting cell type annotation from LLM")
        
        # Extract valid cell type names from candidates for validation
        valid_cell_types = set()
        if candidate_cell_types:
            valid_cell_types = set(candidate_cell_types.keys())
            valid_cell_types.add("Unknown")  # "Unknown" is always valid
            logger.debug(f"Valid candidate cell types: {valid_cell_types}")
        
        for attempt in range(self.max_retries):
            try:
                response = self.client.models.generate_content(
                    model=self.model_name,
                    contents=[
                        {"text": prompt},
                        {
                            "inline_data": {
                                "mime_type": "image/png",
                                "data": image_data
                            }
                        }
                    ],
                    config=types.GenerateContentConfig(
                        temperature=self.temperature,
                    )
                )
                
                # Clean and parse JSON response
                json_text = extract_json_from_response(response.text)
                result = json.loads(json_text)
                
                # Validate response format
                if isinstance(result, list):
                    # Validate each annotation has required fields
                    for ann in result:
                        if 'cluster_id' not in ann or 'cell_type_hypotheses' not in ann:
                            raise ValueError(f"Annotation missing required fields: {ann}")
                elif isinstance(result, dict):
                    # Single annotation, convert to list
                    if 'cluster_id' not in result or 'cell_type_hypotheses' not in result:
                        raise ValueError(f"Annotation missing required fields: {result}")
                else:
                    raise ValueError(f"Unexpected response type: {type(result)}")
                
                # Validate cell types against candidates if provided
                if candidate_cell_types and valid_cell_types:
                    annotations_to_check = result if isinstance(result, list) else [result]
                    
                    for ann in annotations_to_check:
                        cell_type = ann.get('cell_type_hypotheses', '')
                        
                        # Check if the cell type is in valid candidates or "Unknown"
                        is_valid = False
                        
                        if cell_type in valid_cell_types:
                            # Direct match
                            is_valid = True
                        
                        if not is_valid:
                            raise ValueError(
                                f"LLM returned invalid annotation with cell types, for example {cell_type} is not one of given candidate cell types: {valid_cell_types}"
                            )
                        
                
                num_annotations = len(result) if isinstance(result, list) else 1
                logger.info(f"Successfully received valid annotations for {num_annotations} clusters")
                return result
                
            except json.JSONDecodeError as e:
                logger.warning(f"JSON parsing error (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                logger.debug(f"Raw response: {response.text if 'response' in locals() else 'N/A'}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed to parse JSON after {self.max_retries} attempts")
                    raise
            except Exception as e:
                logger.warning(f"Error in cell type annotation (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed after {self.max_retries} attempts")
                    raise
