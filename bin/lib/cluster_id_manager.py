#!/usr/bin/env python3
"""
Cluster ID management utility for generating traceable random cluster identifiers.
"""

import random
import string
import json
import os
from typing import Dict, Optional, List
from .logger import PipelineLogger

logger = PipelineLogger.get_logger(__name__)


class ClusterIDManager:
    """Manages generation and tracking of cluster IDs across iterations."""
    
    def __init__(self, mapping_file: str):
        """
        Initialize the cluster ID manager.
        
        Args:
            mapping_file: Path to JSON file storing cluster ID mappings.
        """
        self.mapping_file = mapping_file
        self.mappings = self._load_mappings()
        self.used_ids = set(self.mappings.keys())
    
    def _load_mappings(self) -> Dict:
        """Load existing cluster ID mappings from file."""
        if os.path.exists(self.mapping_file):
            # Check if file is not empty
            if os.path.getsize(self.mapping_file) > 0:
                with open(self.mapping_file, 'r') as f:
                    return json.load(f)
        return {}
    
    def _save_mappings(self):
        """Save cluster ID mappings to file."""
        with open(self.mapping_file, 'w') as f:
            json.dump(self.mappings, f, indent=2)
        logger.info(f"Saved cluster ID mappings to: {self.mapping_file}")
    
    def generate_id(self) -> str:
        """Generate a unique 6-character alphanumeric ID (public method)."""
        while True:
            # Generate 6-character code (uppercase letters + digits)
            new_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
            if new_id not in self.used_ids:
                self.used_ids.add(new_id)
                return new_id
    
    def generate_multiple_ids(self, count: int) -> List[str]:
        """
        Generate multiple unique random IDs at once.
        
        Args:
            count: Number of IDs to generate.
        
        Returns:
            List of unique random IDs.
        """
        return [self.generate_id() for _ in range(count)]
    
    def register_cluster(
        self,
        cluster_id: str,
        iteration: int,
        parent_id: Optional[str] = None,
        annotation: Optional[Dict] = None
    ):
        """
        Register a single cluster with metadata.
        
        Args:
            cluster_id: The random cluster ID.
            iteration: Current iteration number.
            parent_id: Random ID of parent cluster (if subcluster).
            annotation: Optional dictionary containing LLM annotation results
                       (cell_type_hypothesis, confidence, etc.)
        """
        self.mappings[cluster_id] = {
            "iteration": iteration,
            "parent_id": parent_id,
        }
        
        # Add annotation if provided
        if annotation:
            self.mappings[cluster_id]["annotation"] = annotation
            
        self._save_mappings()
    
    def register_clusters_batch(
        self,
        cluster_ids: List[str],
        iteration: int,
        parent_id: Optional[str] = None,
        annotations: Optional[Dict[str, Dict]] = None
    ):
        """
        Register multiple clusters at once.
        
        Args:
            cluster_ids: List of random cluster IDs.
            iteration: Current iteration number.
            parent_id: Random ID of parent cluster (if subclusters).
            annotations: Optional dictionary mapping cluster_id to annotation dict
        """
        for cluster_id in cluster_ids:
            self.mappings[cluster_id] = {
                "iteration": iteration,
                "parent_id": parent_id,
            }
            
            # Add annotation if provided for this cluster
            if annotations and cluster_id in annotations:
                self.mappings[cluster_id]["annotation"] = annotations[cluster_id]
                
        self._save_mappings()
        logger.info(f"Registered {len(cluster_ids)} clusters for iteration {iteration}")
    
    def update_cluster_annotations(self, annotations: List[Dict]):
        """
        Update clusters with LLM annotation results.
        
        Args:
            annotations: List of annotation dictionaries from LLM containing:
                        cluster, cell_type_hypothesis, confidence, etc.
        """
        for ann in annotations:
            cluster_id = ann.get('cluster')
            if cluster_id and cluster_id in self.mappings:
                # Store relevant annotation fields
                self.mappings[cluster_id]["annotation"] = {
                    "cell_type_hypothesis": ann.get('cell_type_hypothesis'),
                    "cell_type_description": ann.get('cell_type_description'),
                    "confidence": ann.get('confidence'),
                    "split_status": ann.get('split_status'),
                    "key_markers_cited": ann.get('key_markers_cited', []),
                    "key_pathways_cited": ann.get('key_pathways_cited', [])
                }
        self._save_mappings()
        logger.info(f"Updated annotations for {len(annotations)} clusters")
    
    def get_lineage(self, cluster_id: str) -> List[Dict]:
        """
        Get the full lineage (ancestry) of a cluster.
        
        Args:
            cluster_id: Random cluster ID.
        
        Returns:
            List of cluster info dictionaries from root to current cluster.
        """
        lineage = []
        current_id = cluster_id
        
        while current_id:
            if current_id not in self.mappings:
                logger.warning(f"Cluster ID {current_id} not found in mappings")
                break
            
            cluster_info = self.mappings[current_id].copy()
            cluster_info['cluster_id'] = current_id
            lineage.append(cluster_info)
            current_id = cluster_info.get('parent_id')
        
        return list(reversed(lineage))  # Root to current
    
    def export_readable_table(self, output_file: str):
        """
        Export cluster ID mappings as a readable CSV table.
        
        Args:
            output_file: Path to output CSV file.
        """
        import pandas as pd
        
        rows = []
        for cluster_id, info in self.mappings.items():
            lineage = self.get_lineage(cluster_id)
            rows.append({
                'cluster_id': cluster_id,
                'iteration': info['iteration'],
                'parent_id': info.get('parent_id', ''),
                'depth': len(lineage) - 1,
                'lineage': ' → '.join([c['cluster_id'] for c in lineage])
            })
        
        df = pd.DataFrame(rows)
        df = df.sort_values(['iteration', 'cluster_id'])
        df.to_csv(output_file, index=False)
        logger.info(f"Exported cluster ID table to: {output_file}")
        logger.info(f"Exported cluster ID table to: {output_file}")
