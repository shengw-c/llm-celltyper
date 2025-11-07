#!/usr/bin/env python3
"""
Script to aggregate hierarchical annotation JSON files and create consolidated outputs.

This script reads all *_cell_annotation.json files and *_subset.tsv files from the
annotation pipeline, builds a hierarchical tree structure, and generates:
1. A consolidated JSON with nested annotations and unique IDs
2. A TSV mapping cells to their cluster unique IDs and annotations at each level

Usage:
    python aggregate_annotations.py \\
        --responses-dir DIR \\
        --data-dir DIR \\
        --output-json FILE \\
        --output-tsv FILE
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict
import argparse


class AnnotationNode:
    """Represents a node in the annotation hierarchy tree."""
    
    def __init__(self, nametag: str, data: Dict):
        self.nametag = nametag
        self.data = data
        self.children: List[AnnotationNode] = []
        
    def to_dict(self) -> Dict:
        """Convert node to dictionary format with nested children."""
        result = dict(self.data)
        
        # Add unique IDs and process children recursively
        if 'final_annotations' in result:
            self._add_unique_ids_and_children(result, self.nametag, self.children)
        
        return result
    
    def to_flat_dict(self) -> Dict[str, Dict]:
        """
        Convert tree to a flat dictionary mapping unique_id -> annotation.
        
        Returns:
            Dictionary where keys are unique_ids and values are annotation objects
        """
        flat_dict = {}
        
        # First, add unique IDs using the tree structure
        temp_result = dict(self.data)
        if 'final_annotations' in temp_result:
            self._add_unique_ids_and_children(temp_result, self.nametag, self.children)
        
        # Now flatten the tree structure
        def flatten_annotations(annotations: List[Dict], parent_cell_type: str = "Root"):
            """Recursively flatten annotations and their children."""
            for annot in annotations:
                unique_id = annot.get('unique_id', '')
                if unique_id:
                    # Create a copy without the 'children' field for the flat structure
                    flat_annot = {k: v for k, v in annot.items() if k != 'children'}
                    
                    # Add parent cell type (passed from the recursion)
                    flat_annot['parent_cell_type'] = parent_cell_type
                    
                    flat_dict[unique_id] = flat_annot
                    
                    # Process children recursively, passing this annotation's cell type as parent
                    if 'children' in annot and annot['children']:
                        current_cell_type = annot.get('cell_type', 'Unknown')
                        flatten_annotations(annot['children'], current_cell_type)
        
        if 'final_annotations' in temp_result:
            flatten_annotations(temp_result['final_annotations'], "Root")
        
        return flat_dict
    
    def _add_unique_ids_and_children(self, result: Dict, current_nametag: str, 
                                      child_nodes: List['AnnotationNode'], parent_unique_id: str = ""):
        """
        Recursively add unique IDs and nest children in the annotation response.
        
        Args:
            result: The dictionary containing final_annotations to process
            current_nametag: The nametag of the current level
            child_nodes: List of child AnnotationNode objects
            parent_unique_id: The unique_id prefix from parent (e.g., "0" or "0.1")
        """
        if 'final_annotations' not in result:
            return
        
        # Add unique IDs to annotations at this level
        for annotation in result['final_annotations']:
            cluster_id = annotation['cluster_id']
            
            # Create unique hierarchical ID
            if parent_unique_id:
                unique_id = f"{parent_unique_id}.{cluster_id}"
            else:
                unique_id = cluster_id
            
            annotation['unique_id'] = unique_id
            annotation['hierarchy_path'] = current_nametag
        
        # Process child nodes
        for child_node in child_nodes:
            # Extract the cell type this child represents
            child_cell_type = self._extract_cell_type_from_nametag(child_node.nametag, current_nametag)
            
            # Find parent annotation(s) that match this cell type
            for annotation in result['final_annotations']:
                if annotation['cell_type'] == child_cell_type:
                    # This annotation should contain the child's annotations
                    parent_unique_id_for_children = annotation['unique_id']
                    
                    # Get child data and process it recursively
                    child_result = {'final_annotations': list(child_node.data.get('final_annotations', []))}
                    
                    # Recursively add unique IDs to child annotations
                    self._add_unique_ids_and_children(
                        child_result, 
                        child_node.nametag, 
                        child_node.children,
                        parent_unique_id_for_children
                    )
                    
                    # Add the processed children to this annotation
                    annotation['children'] = child_result['final_annotations']
                    break
    
    def _extract_cell_type_from_nametag(self, child_nametag: str, parent_nametag: str) -> str:
        """
        Extract the cell type name from a child nametag relative to parent.
        
        Examples:
            child: lev0_Immune_Cell, parent: lev0_all_types -> "Immune Cell"
            child: lev0_Immune_Cell_Myeloid_Cell, parent: lev0_Immune_Cell -> "Myeloid Cell"
        """
        if parent_nametag.endswith('_all_types'):
            # Child is like lev0_Immune_Cell, extract "Immune_Cell" and convert to "Immune Cell"
            parts = child_nametag.split('_', 1)  # Split on first underscore only
            if len(parts) > 1:
                cell_type_underscored = parts[1]
                return cell_type_underscored.replace('_', ' ')
        else:
            # Child has parent as prefix, extract the additional part
            prefix = parent_nametag + '_'
            if child_nametag.startswith(prefix):
                remaining = child_nametag[len(prefix):]
                return remaining.replace('_', ' ')
        
        return ""


def parse_nametag(filename: str) -> Optional[str]:
    """
    Extract nametag from filename.
    
    Example: 
        lev0_Immune_Cell_Myeloid_Cell_cell_annotation.json 
        -> lev0_Immune_Cell_Myeloid_Cell
    """
    if not filename.endswith('_cell_annotation.json'):
        return None
    
    # Remove the suffix
    nametag = filename.replace('_cell_annotation.json', '')
    return nametag


def get_hierarchy_path(nametag: str) -> List[str]:
    """
    Extract hierarchy path from nametag.
    
    Example:
        lev0_Immune_Cell_Myeloid_Cell_Macrophage
        -> ['lev0_all_types', 'lev0_Immune_Cell', 'lev0_Immune_Cell_Myeloid_Cell', 
            'lev0_Immune_Cell_Myeloid_Cell_Macrophage']
    """
    parts = nametag.split('_')
    
    # First part is the level (e.g., lev0)
    level_prefix = parts[0]
    
    # Remaining parts are the cell type hierarchy
    # We need to rebuild the path incrementally
    if len(parts) == 2:  # e.g., lev0_all_types
        return [nametag]
    
    # Build the path by accumulating parts
    path = []
    current_path = level_prefix
    
    for i in range(1, len(parts)):
        current_path += '_' + parts[i]
        path.append(current_path)
    
    # Prepend the root (all_types) if this is not already it
    if nametag != f'{level_prefix}_all_types':
        path.insert(0, f'{level_prefix}_all_types')
    
    return path


def cell_type_to_nametag_part(cell_type: str) -> str:
    """
    Convert a cell type name to its nametag representation.
    
    Example:
        "Airway Epithelium" -> "Airway_Epithelium"
        "T Cell" -> "T_Cell"
    """
    return cell_type.replace(' ', '_').replace('-', '-')


def find_parent_by_expected_types(child_nametag: str, all_annotation_data: Dict[str, Dict]) -> Optional[str]:
    """
    Find the parent of a child nametag by checking which annotation's
    expected_cell_types contains the child's cell type.
    
    The naming convention is:
    - Root: lev0_all_types
    - Level 1 children: lev0_{CellType} (e.g., lev0_Immune_Cell)
    - Deeper levels: lev0_{Parent}_{Child} (e.g., lev0_Immune_Cell_Myeloid_Cell)
    
    Examples:
        lev0_Immune_Cell -> parent is lev0_all_types (expected_cell_types contains "Immune Cell")
        lev0_Immune_Cell_Myeloid_Cell -> parent is lev0_Immune_Cell (expected_cell_types contains "Myeloid Cell")
        lev0_Epithelial_Cell_Airway_Epithelium -> parent is lev0_Epithelial_Cell
    """
    # Special case: root
    if child_nametag.endswith('_all_types'):
        return None
    
    # Extract level prefix
    parts = child_nametag.split('_')
    level_prefix = parts[0]
    
    # Find all potential parents
    # Sort by length descending to find the closest (most specific) parent first
    potential_parents = sorted(all_annotation_data.keys(), key=len, reverse=True)
    
    for parent_nametag in potential_parents:
        if parent_nametag == child_nametag:
            continue
            
        parent_data = all_annotation_data[parent_nametag]
        
        # Check if this parent's expected_cell_types contains a type that matches the child
        if 'expected_cell_types' not in parent_data:
            continue
        
        # Determine what part of the child represents the cell type we're looking for
        # Case 1: Child is lev0_CellType, parent should be lev0_all_types
        # Case 2: Child is lev0_Parent_CellType, parent should be lev0_Parent
        
        # Check if child's nametag contains any of the parent's expected types
        for expected_type in parent_data['expected_cell_types']:
            expected_nametag = cell_type_to_nametag_part(expected_type)
            
            # For root parent (all_types), check if child is lev{level}_{expected_type}
            if parent_nametag.endswith('_all_types'):
                expected_child = f"{level_prefix}_{expected_nametag}"
                if child_nametag == expected_child or child_nametag.startswith(expected_child + '_'):
                    return parent_nametag
            
            # For non-root parents, check if child starts with parent and then has expected type
            elif child_nametag.startswith(parent_nametag + '_'):
                remaining = child_nametag[len(parent_nametag) + 1:]
                if remaining == expected_nametag or remaining.startswith(expected_nametag + '_'):
                    return parent_nametag
    
    return None


def load_annotation_files(input_dir: str) -> Dict[str, Dict]:
    """Load all annotation JSON files from the input directory."""
    annotation_data = {}
    input_path = Path(input_dir)
    
    for json_file in input_path.glob('*_cell_annotation.json'):
        nametag = parse_nametag(json_file.name)
        if nametag:
            with open(json_file, 'r') as f:
                data = json.load(f)
                annotation_data[nametag] = data
                print(f"Loaded: {nametag}")
    
    return annotation_data


def build_tree(annotation_data: Dict[str, Dict]) -> Optional[AnnotationNode]:
    """
    Build a tree structure from annotation data.
    
    Returns the root node of the tree.
    """
    # Create nodes for all annotations
    nodes: Dict[str, AnnotationNode] = {}
    
    for nametag, data in annotation_data.items():
        nodes[nametag] = AnnotationNode(nametag, data)
    
    # Build parent-child relationships
    root = None
    orphans = []
    
    for nametag, node in nodes.items():
        # Find parent using expected_cell_types method
        parent_nametag = find_parent_by_expected_types(nametag, annotation_data)
        
        if parent_nametag is None:
            # This is a root node
            if root is not None:
                print(f"Warning: Multiple root nodes found: {root.nametag} and {nametag}")
            root = node
        elif parent_nametag in nodes:
            # Add as child to parent
            nodes[parent_nametag].children.append(node)
        else:
            # Parent not found - orphan node
            orphans.append(nametag)
            print(f"Warning: Parent '{parent_nametag}' not found for '{nametag}'")
    
    if orphans:
        print(f"\nFound {len(orphans)} orphan nodes (no parent found):")
        for orphan in orphans[:10]:  # Show first 10
            print(f"  - {orphan}")
        if len(orphans) > 10:
            print(f"  ... and {len(orphans) - 10} more")
    
    return root


def build_hierarchy_level_mapping(tree_json_path: str) -> Dict[str, int]:
    """
    Build a mapping from cell type names to their hierarchy levels using the tree JSON.
    
    Args:
        tree_json_path: Path to the tree JSON file defining the hierarchy
        
    Returns:
        Dictionary mapping cell type names to their hierarchy levels
    """
    if not tree_json_path or not Path(tree_json_path).exists():
        return {}
    
    with open(tree_json_path, 'r') as f:
        tree_data = json.load(f)
    
    level_mapping = {}
    
    def traverse_tree(node: Dict, current_level: int = 0, parent_path: str = ""):
        """Recursively traverse the tree and assign levels."""
        for cell_type, data in node.items():
            # Store the mapping
            level_mapping[cell_type] = current_level
            
            # Process children if they exist
            if isinstance(data, dict) and 'children' in data and data['children']:
                traverse_tree(data['children'], current_level + 1, cell_type)
    
    # Start traversal from the root level
    traverse_tree(tree_data)
    
    return level_mapping


def calculate_hierarchy_level_from_unique_id(unique_id: str) -> int:
    """
    Calculate the hierarchy level from a unique_id.
    
    The unique_id format is like "0", "0.1", "0.1.2", "4.7", etc.
    The level is the number of dots plus 1 (but we use 0-indexed levels).
    Actually, the level is simply the number of dots.
    
    Args:
        unique_id: The unique identifier (e.g., "2.0.1" means level 2 in the hierarchy)
        
    Returns:
        The hierarchy level (0-indexed)
    """
    if not unique_id or unique_id == '':
        return 0
    
    # Count the number of dots to determine depth
    # "0" = level 0
    # "0.1" = level 1
    # "0.1.2" = level 2
    return unique_id.count('.')


def process_tsv_files(tsv_dir: str, consolidated_data: Dict[str, Dict], output_tsv: str, tree_json_path: str = None):
    """
    Process TSV files to map cell IDs to unique cluster IDs and create combined output.
    
    Args:
        tsv_dir: Directory containing *_subset.tsv files
        consolidated_data: Flat dictionary mapping unique_id -> annotation data
        output_tsv: Path to output TSV file
    """
    import pandas as pd
    import numpy as np
    from pathlib import Path
    
    # Build a mapping from (nametag, cluster_id) -> unique_id and annotation
    cluster_mapping = {}
    
    # Extract mappings from flat consolidated data
    for unique_id, annot in consolidated_data.items():
        cluster_id = annot.get('cluster_id', '')
        cell_type = annot.get('cell_type', '')
        hierarchy_path = annot.get('hierarchy_path', '')
        
        # Store mapping using (hierarchy_path, cluster_id) as key
        if hierarchy_path and cluster_id:
            key = (hierarchy_path, cluster_id)
            cluster_mapping[key] = {
                'unique_id': unique_id,
                'cell_type': cell_type,
                'hierarchy_path': hierarchy_path
            }
    
    # Build hierarchy level mapping from tree JSON if provided
    hierarchy_levels = {}
    if tree_json_path:
        print(f"  Loading hierarchy levels from: {tree_json_path}")
        hierarchy_levels = build_hierarchy_level_mapping(tree_json_path)
        print(f"  Built hierarchy mapping for {len(hierarchy_levels)} cell types")
    
    print(f"  Built mapping for {len(cluster_mapping)} cluster annotations")
    
    # Process all TSV files
    tsv_path = Path(tsv_dir)
    tsv_files = list(tsv_path.glob('*_subset.tsv'))
    
    if not tsv_files:
        print(f"  No TSV files found in {tsv_dir}")
        return
    
    print(f"  Found {len(tsv_files)} TSV files to process")
    
    # Dictionary to store cell annotations by cell_id
    cell_annotations = {}
    
    for tsv_file in tsv_files:
        # Extract nametag from filename
        nametag = tsv_file.stem.replace('_subset', '')
        
        # Read TSV file
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Handle missing first column header (some files have empty first column name)
        if len(df.columns) >= 3 and (df.columns[0] == '' or df.columns[0].strip() == '' or df.columns[0].startswith('Unnamed:')):
            # Fix missing first column header
            new_columns = list(df.columns)
            new_columns[0] = 'cellid'
            df.columns = new_columns
        
        # Column names: cellid, leiden_X_XX, ann
        # We need to use the second column (leiden) as cluster_id
        if len(df.columns) < 3:
            print(f"  Warning: Skipping {tsv_file.name} - insufficient columns")
            continue
        
        cell_id_col = df.columns[0]
        cluster_col = df.columns[1]
        ann_col = df.columns[2]
        
        # Process each row
        for _, row in df.iterrows():
            cell_id = row[cell_id_col]
            cluster_id = str(row[cluster_col])  # Convert to string to match JSON
            annotation = row[ann_col]
            
            # Look up unique_id and cell type
            key = (nametag, cluster_id)
            if key in cluster_mapping:
                mapping_info = cluster_mapping[key]
                unique_id = mapping_info['unique_id']
                cell_type = mapping_info['cell_type']
                hierarchy_path = mapping_info['hierarchy_path']
                
                # Determine annotation level
                # The most reliable way is to use the unique_id, which encodes the hierarchy depth
                level = calculate_hierarchy_level_from_unique_id(unique_id)
                
                # Fallback: if unique_id is not available or invalid, try other methods
                if level is None or level < 0:
                    if hierarchy_levels and cell_type in hierarchy_levels:
                        # Use tree-based level mapping
                        level = hierarchy_levels[cell_type]
                    else:
                        # Last resort: use a simple heuristic
                        # Root level has "_all_types"
                        if hierarchy_path.endswith('_all_types'):
                            level = 0
                        else:
                            # Count hierarchy segments by finding known cell types in the path
                            # This is a rough estimate
                            level = 1
                
                # Store in cell_annotations
                if cell_id not in cell_annotations:
                    cell_annotations[cell_id] = {}
                
                # Store annotation at this level
                level_key = f'ann_level{level}'
                cell_annotations[cell_id][level_key] = cell_type
                
                # Update unique_id if this is deeper (more specific)
                if 'unique_id' not in cell_annotations[cell_id]:
                    cell_annotations[cell_id]['unique_id'] = unique_id
                else:
                    # Keep the more specific (longer) unique_id
                    current_id = cell_annotations[cell_id]['unique_id']
                    if len(unique_id.split('.')) > len(current_id.split('.')):
                        cell_annotations[cell_id]['unique_id'] = unique_id
    
    print(f"  Mapped {len(cell_annotations)} cells to annotations")
    
    # Convert to DataFrame
    rows = []
    for cell_id, annotations in cell_annotations.items():
        row = {'cell_id': cell_id}
        row.update(annotations)
        rows.append(row)
    
    result_df = pd.DataFrame(rows)
    
    # Sort columns: cell_id, then ann_level0, ann_level1, ..., unique_id
    # Use numeric sorting for level columns and ensure sequential levels
    level_cols = [c for c in result_df.columns if c.startswith('ann_level')]
    if level_cols:
        # Extract level numbers and find the max level
        level_numbers = [int(c.replace('ann_level', '')) for c in level_cols]
        max_level = max(level_numbers)
        
        # Create complete sequential list of level columns (even if some are missing)
        complete_level_cols = [f'ann_level{i}' for i in range(max_level + 1)]
        
        # Add missing level columns with NaN values
        for col in complete_level_cols:
            if col not in result_df.columns:
                result_df[col] = None
        
        # Use the complete sequential list
        level_cols = complete_level_cols
    
    ordered_cols = ['cell_id'] + level_cols + ['unique_id']
    result_df = result_df[ordered_cols]
    
    # Add ann_finest column (last non-NaN annotation)
    ann_cols = [c for c in level_cols]
    if ann_cols:
        result_df['ann_finest'] = result_df[ann_cols].ffill(axis=1).iloc[:, -1]
        ordered_cols = ['cell_id'] + level_cols + ['ann_finest', 'unique_id']
        result_df = result_df[ordered_cols]
    
    # Sort by cell_id
    result_df = result_df.sort_values('cell_id')
    
    # Save to file
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"\n✓ TSV consolidation complete!")
    print(f"  Output file: {output_tsv}")
    print(f"  Total cells: {len(result_df)}")
    print(f"  Annotation levels: {len(level_cols)}")


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate hierarchical annotation JSON files into a single consolidated structure.'
    )
    parser.add_argument(
        '--responses-dir',
        required=True,
        help='Input directory containing annotation JSON files (*_cell_annotation.json)'
    )
    parser.add_argument(
        '--data-dir',
        required=True,
        help='Input directory containing TSV annotation files (*_subset.tsv)'
    )
    parser.add_argument(
        '--output-json',
        required=True,
        help='Output JSON file path for consolidated annotations'
    )
    parser.add_argument(
        '--output-tsv',
        required=True,
        help='Output TSV file path for cell-to-cluster mapping'
    )
    parser.add_argument(
        '--tree-json',
        required=False,
        help='Optional tree JSON file defining the hierarchical structure for level calculation'
    )
    parser.add_argument(
        '--pretty',
        action='store_true',
        help='Pretty print the output JSON (with indentation)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("Aggregating Annotations")
    print("="*80)
    print(f"Responses directory: {args.responses_dir}")
    print(f"Data directory     : {args.data_dir}")
    print(f"Output JSON        : {args.output_json}")
    print(f"Output TSV         : {args.output_tsv}")
    print("="*80)
    
    # Load annotation files
    print(f"\n[1/4] Loading annotation files from: {args.responses_dir}")
    annotation_data = load_annotation_files(args.responses_dir)
    
    if not annotation_data:
        print("ERROR: No annotation files found!")
        return 1
    
    print(f"      ✓ Loaded {len(annotation_data)} annotation files")
    
    # Build tree structure
    print("\n[2/4] Building tree structure...")
    root = build_tree(annotation_data)
    
    if root is None:
        print("ERROR: No root node found!")
        return 1
    
    print(f"      ✓ Root node: {root.nametag}")
    print(f"      ✓ Total nodes in tree: {len(annotation_data)}")
    
    # Convert tree to flat dictionary
    print("\n[3/4] Converting tree to flat dictionary...")
    consolidated_data_flat = root.to_flat_dict()
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write JSON output
    print(f"      ✓ Writing consolidated annotations to: {args.output_json}")
    with open(args.output_json, 'w') as f:
        if args.pretty:
            json.dump(consolidated_data_flat, f, indent=2)
        else:
            json.dump(consolidated_data_flat, f)
    
    total_clusters = len(consolidated_data_flat)
    print(f"      ✓ Total annotations: {total_clusters}")
    
    # Process TSV files
    print(f"\n[4/4] Processing TSV files from: {args.data_dir}")
    if not Path(args.data_dir).exists():
        print(f"WARNING: Data directory not found: {args.data_dir}")
        print("Skipping TSV processing.")
        return 1
    
    process_tsv_files(args.data_dir, consolidated_data_flat, args.output_tsv, args.tree_json)
    
    print("\n" + "="*80)
    print("✓ Aggregation complete!")
    print("="*80)
    print(f"  JSON output: {args.output_json}")
    print(f"  TSV output : {args.output_tsv}")
    print("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
