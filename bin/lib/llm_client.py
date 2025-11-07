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
import numpy as np
import pandas as pd
from .logger import PipelineLogger

from .prompts import (
    cluster_PROMPT,
    Celltyper_Instruction,
    Consolidate_Instruction
)

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

    def select_cluster_resolution(self, except_cell_types: int, conditions: int, cluster_adjacency_matrix_dict: Dict, custom_model: str = None, transition_cutoff: float = 0.1) -> Dict:
        """
        Queries the Gemini model to determine the optimal clustering resolution 
        from a cluster tree image.
        
        Args:
            except_cell_types: Number of cell types to exclude from clustering
            conditions: Number of conditions to consider for clustering
            cluster_adjacency_matrix_dict: Dictionary containing the cluster adjacency matrix data.
            custom_model: Optional custom model name to override default.
            transition_cutoff: Cutoff value to consider a transition between clusters.
        
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
                ## try to do manual determination before sending to LLM
                logger.info(f"Try manual determination of cluster resolution")
                target_num_clusters = min(max(2, len(except_cell_types)) + conditions, max(2, len(except_cell_types)) * 2)
                stables = {}
                substables = {}
                for key in cluster_adjacency_matrix_dict.keys():
                    dat = np.array(cluster_adjacency_matrix_dict[key]["data"])
                    split_resolutions = ((dat>=transition_cutoff).sum(1)>1).any()
                    merge_resoultions = ((dat>=transition_cutoff).sum(0)>1).any()
                    if not (split_resolutions or merge_resoultions):
                        if dat.shape[0]>= target_num_clusters: 
                            stables[key] = dat.shape[0]
                    if not merge_resoultions:
                        if dat.shape[0]>= target_num_clusters:
                            substables[key] = dat.shape[0]
                valid_stables = {key:value for key, value in stables.items() if value>=target_num_clusters}
                valid_substables = {key:value for key, value in substables.items() if value>=target_num_clusters}
                if valid_stables:
                    valid_stables = {float(key.split("_to_")[0].replace("leiden_", "").replace("_", ".")):value for key,value in valid_stables.items()}
                    sel_resolution = min(valid_stables.keys())

                    logger.info(f"Manual determination successful. Selected resolution: {sel_resolution} with {valid_stables[sel_resolution]} clusters.")
                    return {
                        "resolution": sel_resolution,
                        "justification": f"Manual determination found stable resolution with {valid_stables[sel_resolution]} clusters, meeting target of {target_num_clusters}."
                    }
                elif valid_substables:
                    valid_substables = {float(key.split("_to_")[0].replace("leiden_", "").replace("_", ".")):value for key,value in valid_substables.items()}
                    sel_resolution = min(valid_substables.keys())

                    logger.info(f"Manual determination successful (substable). Selected resolution: {sel_resolution} with {valid_substables[sel_resolution]} clusters.")
                    return {
                        "resolution": sel_resolution,
                        "justification": f"Manual determination found substable resolution with {valid_substables[sel_resolution]} clusters, meeting target of {target_num_clusters}."
                    }
                else:
                    logger.info(f"Manual determination failed. Proceeding with LLM query.")
                    prompt = [
                        cluster_PROMPT.format(
                            cell_type_num=len(except_cell_types),
                            condition_num=conditions,
                            transition_cutoff=transition_cutoff
                        ),
                        json.dumps(cluster_adjacency_matrix_dict)
                    ]
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
    
    def procCellAnnotationResults(self, res: types.GenerateContentResponse) -> List:

        result = json.loads(res.text)
        cluster_states = {}  # cluster_id -> {'high_conf': bool, 'annotations': [ann]}
        clusters_high_conf, clusters_to_replicate = [], []

        for cluster_ann in result:
            cluster_id = cluster_ann.get('cluster_id')
            confidence = cluster_ann.get('confidence', 'Low')
            cell_type = cluster_ann.get('cell_type_hypothesis', 'Unknown')
            alternative_hypothesis = cluster_ann.get('alternative_hypothesis', None)

            # Check if High confidence AND not Unknown
            if confidence == 'High' and not (cell_type == 'Unknown' and alternative_hypothesis is None):
                clusters_high_conf.append(cluster_id)
            else:
                clusters_to_replicate.append(cluster_id)
        return clusters_high_conf, clusters_to_replicate

    def annotate_cell_types(self, 
                            context: str,
                            candidate_cell_types: List[dict],
                            marker_genes: Dict[str, List[str]],
                            pathways: Dict[str, List[str]],
                            cluster_adjacency: Dict[str, List[str]],
                            parent_cell_type: Optional[str] = None,
                            nreplicate: int = 3,
                            custom_model: str = None, 
                            custom_model_consolidation: Optional[str] = None, 
                            adaptive: bool = True) -> Dict:
        """
        Queries the Gemini model to perform cell type annotation based on marker genes 
        and UMAP visualization.
        
        Args:
            context: General context information for the annotation task
            parent_cell_type: Parent cell type for hierarchical annotation (if any)
            candidate_cell_types: List of candidate cell types with definitions
            marker_genes: Dictionary mapping cluster IDs to lists of marker genes
            pathways: List of pathway enrichment results per cluster
            cluster_adjacency: Cluster adjacency information
            nreplicate: Number of replicate annotations to request for robustness.
            custom_model: Optional custom model name to override default.
            custom_model_consolidation: Optional custom model name for consolidation step, only effective if nreplicate>1.
            adaptive: If True, only run replicates for clusters with Low/Medium confidence or Unknown annotations.
                     If False, run nreplicate for all clusters.

        Returns:
            Tuple of (result, res_df_dict) where:
            - result: List of dictionaries containing cell type annotations for each cluster.
                Each dict has: {'cluster_id': str, 'cell_type_hypotheses': str, 
                'justification': str, 'key_markers_cited': list, 
                'confidence': str}
            - res_df_dict: Dictionary of replicate annotations by cluster (None if no replicates were run).
                {'cluster_id': str, 'cell_type': str, 
                'justification': str, 'key_marker_genes': list, 
                'confidence': str}

        Raises:
            Exception: If all retry attempts fail
        """
        model_to_use = custom_model if custom_model else self.model_name
        logger.info(f"Requesting cell type annotation from LLM using {model_to_use}")
        if adaptive and nreplicate > 1:
            logger.info(f"Adaptive mode enabled: will only replicate Low/Medium confidence or Unknown annotations")
        
        for attempt in range(self.max_retries):
            try:
                initial_prompt = Celltyper_Instruction.format(
                    expr_context=context,
                    candidate_cell_types=json.dumps(candidate_cell_types),
                    marker_genes_json=json.dumps(marker_genes),
                    pathway_json=json.dumps(pathways),
                    cluster_adjacency_json=json.dumps(cluster_adjacency)
                )
                input_token = self.client.models.count_tokens(
                    model=model_to_use,
                    contents=initial_prompt
                )
                logger.info(f"[Cell annotation] Input prompt token count: {input_token.total_tokens}")
                
                # First pass: always run once to get initial annotations
                total_tokens = 0
                logger.info(f"Generating initial annotation pass")
                first_ann_res = self.client.models.generate_content(
                    model=model_to_use,
                    contents=initial_prompt,
                    config={
                        "response_mime_type": self.response_mime_type
                    }
                )
                total_tokens += first_ann_res.usage_metadata.total_token_count
                logger.info(f"[Cell annotation] Initial pass done. Tokens used: {first_ann_res.usage_metadata.total_token_count}")
                
                # Adaptive logic: determine which clusters need replication
                if adaptive and nreplicate > 1:
                    clusters_high_conf, clusters_to_replicate = self.procCellAnnotationResults(first_ann_res)

                    logger.info(f"[Cell annotation - Adaptive] {len(clusters_high_conf)} clusters with High confidence initially: {clusters_high_conf}")
                    logger.info(f"[Cell annotation - Adaptive] {len(clusters_to_replicate)} clusters need replication: {clusters_to_replicate}")
                    
                    # Store all results (first pass for all clusters)
                    res = [first_ann_res.text]
                    
                    # Dynamic early stopping: replicate until all clusters reach High confidence or max replicates
                    remaining_nreplicate = nreplicate - 1
                    logger.info(f"[Cell annotation - Adaptive] Starting adaptive replication using model {model_to_use}...")
                    while len(clusters_to_replicate) > 0 and remaining_nreplicate > 0:
                        logger.info(f"[Cell annotation - Adaptive] Replicating {len(clusters_to_replicate)} clusters: round {nreplicate - remaining_nreplicate + 1}/{nreplicate}")
                        remaining_nreplicate -= 1

                        ## update prompts with only clusters needing replication
                        new_top_genes_by_specificity = {
                            key:value 
                            for key,value in marker_genes.items() 
                            if key in clusters_to_replicate
                        }
                        new_pathways = [
                            item
                            for item in pathways
                            if item['cluster'] in clusters_to_replicate
                        ]
                        
                        new_prompt = Celltyper_Instruction.format(
                            expr_context=context,
                            candidate_cell_types=json.dumps(candidate_cell_types),
                            marker_genes_json=json.dumps(new_top_genes_by_specificity),
                            pathway_json=json.dumps(new_pathways),
                            cluster_adjacency_json=json.dumps(cluster_adjacency)
                        )
                        new_ann_res = self.client.models.generate_content(
                            model=model_to_use,
                            contents=new_prompt,
                            config={
                                "response_mime_type": self.response_mime_type
                            }
                        )
                        res.append(new_ann_res.text)
                        total_tokens += new_ann_res.usage_metadata.total_token_count
                        logger.info(f"[Cell annotation] {total_tokens} total tokens used until now.")

                        # Parse new annotations and update cluster states
                        current_high_conf, current_to_replicate = self.procCellAnnotationResults(new_ann_res)
                        clusters_to_replicate = current_to_replicate

                ## output total tokens used after adaptive replication
                logger.info(f"[Cell annotation - Adaptive] Adaptive replication done. Total tokens used: {total_tokens}")
                # After adaptive replication, check if any clusters still need replication
                if clusters_to_replicate:
                    logger.info(f"[Cell annotation - Adaptive] {len(clusters_to_replicate)} clusters still have not reached High confidence: {clusters_to_replicate}")

                ## combine all annotations into pandas DataFrame
                res_df = []
                for i in range(len(res)):
                    r = res[i]
                    res_df.append(pd.DataFrame(json.loads(r)).assign(rep = i))
                res_df = pd.concat(res_df).sort_values("cluster_id")
                res_df.confidence = pd.Categorical(res_df.confidence, categories=["High", "Medium", "Low"])

                ## part1: Only keep one high confidence annotation for each cluster
                ## Keep all annotations for clusters that never reached high confidence
                cluster_done = res_df.query("confidence=='High'").sort_values("rep").groupby("cluster_id").head(1)
                cluster_done = cluster_done.assign(cell_type = lambda x: np.where(x.cell_type_hypothesis=='Unknown', x.alternative_hypothesis, x.cell_type_hypothesis))
                ## remove any clusters that are still Unknown after replicated calls
                cluster_done = cluster_done[~cluster_done.cluster_id.isin(clusters_to_replicate)]

                ## part2: Keep all annotations for clusters that never reached high confidence
                cluster_toassign_df = res_df[~res_df.cluster_id.isin(cluster_done.cluster_id)]
                cluster_toassign = {}
                for cluster_id, group_df in cluster_toassign_df.groupby("cluster_id"):
                    cluster_toassign[cluster_id] = {f"justification{ix+1}":text for ix, text in enumerate(group_df.justification.tolist())}

                
                ## run consolidation if multiple replicates
                if cluster_toassign:
                    custom_model_consolidation = custom_model_consolidation if custom_model_consolidation else self.model_name
                    logger.info(f"[Cell annotation - Consolidation] Consolidating {len(cluster_toassign)} annotation replicates with model {custom_model_consolidation}")

                    logger.info(f"[Cell annotation - Consolidation] Clusters to consolidate: {json.dumps(cluster_toassign, indent=2)}")
                    
                    ## call consolidation model
                    con_res = self.client.models.generate_content(
                        model=model_to_use,
                        contents=Consolidate_Instruction.format(
                            expert_justifications=json.dumps(cluster_toassign),
                            candidate_cell_types=json.dumps(candidate_cell_types),
                            parent_cell_type=parent_cell_type
                        ),
                        config={
                            "response_mime_type": self.response_mime_type
                        }
                    )
                    logger.info(f"[Cell annotation - Consolidation] Consolidation done. Tokens used: {con_res.usage_metadata.total_token_count}")

                    ## make results for the consolidated clusters
                    cluster_toassign_df = (
                        cluster_toassign_df
                        .groupby("cluster_id", as_index=False).agg({'key_markers_cited': "sum"})
                        .merge(
                            pd.DataFrame(json.loads(con_res.text))
                        )
                    )

                    ## combine with done clusters
                    res_df = pd.concat([cluster_done, cluster_toassign_df], ignore_index=True)[["cluster_id", "cell_type", "key_markers_cited", "justification", "cell_status", "confidence"]]
                else:
                    ## all clusters are assigned with high confidence
                    res_df = cluster_done[["cluster_id", "cell_type", "key_markers_cited", "justification", "cell_status", "confidence"]]

                res_df_dict = res_df.to_dict(orient="records")
                logger.info(f"Annotation cell types for {len(res_df_dict)} clusters done.")
                return res, res_df_dict
                
            except json.JSONDecodeError as e:
                logger.warning(f"JSON parsing error (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                logger.debug(f"Raw response: {res[0].text if 'res' in locals() and res[0] else 'N/A'}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed to parse JSON after {self.max_retries} attempts")
                    raise
            except Exception as e:
                logger.warning(f"Error in cell type annotation (attempt {attempt + 1}/{self.max_retries}): {str(e)}")
                if attempt == self.max_retries - 1:
                    logger.error(f"Failed after {self.max_retries} attempts")
                    raise