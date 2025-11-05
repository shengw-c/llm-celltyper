#!/usr/bin/env python3
"""
Utility functions for working with cell type hierarchy trees.
"""

from typing import Dict, Union, Optional
import json


class TreeValidationError(Exception):
    """Raised when tree structure is invalid or node cannot be found."""
    pass


def find_node_data(tree_dict: Dict, cell_name: str, raise_on_missing: bool = False) -> Optional[Dict]:
    """
    Recursively searches the entire tree for a cell by its name.
    
    Args:
        tree_dict (Dict): The dictionary to search in (e.g., the full tree
                         or a 'children' sub-dictionary).
        cell_name (str): The name of the cell to find.
        raise_on_missing (bool): If True, raise TreeValidationError when node not found.
                                If False, return None (default behavior).

    Returns:
        Optional[Dict]: The data dictionary for the found cell, or None if not found.
        
    Raises:
        TreeValidationError: If raise_on_missing=True and node is not found.
    """
    # Check if the cell_name is a key in the current dictionary
    if cell_name in tree_dict:
        return tree_dict[cell_name]
    
    # If not, iterate through all nodes in this dictionary
    # and search their children recursively
    for node_name, node_data in tree_dict.items():
        children_dict = node_data.get("children", {})
        
        # Recursively search within the children dictionary
        found = find_node_data(children_dict, cell_name, raise_on_missing=False)
        
        # If the recursive call found it, return the result
        if found:
            return found
            
    # If we've searched everywhere and not found it
    if raise_on_missing:
        raise TreeValidationError(
            f"Cell type '{cell_name}' not found in hierarchy tree. "
            f"Available root nodes: {list(tree_dict.keys())}"
        )
    return None


def get_immediate_children(node_data: Optional[Dict], include_parent: bool = False, parent_name: str = "") -> Dict[str, Dict[str, Union[str, bool]]]:
    """
    Extracts only the direct children of a node from the cell type hierarchy tree.
    
    Args:
        node_data (Optional[Dict]): The data dictionary for a single cell type node.
                                   Expected to contain a 'children' key with nested cell type data.
        include_parent (bool): If True, include the parent node in the output (default: False).
        parent_name (str): The name of the parent node (used if include_parent is True).
    
    Returns:
        Dict[str, Dict[str, Union[str, bool]]]: A dictionary mapping child cell type names to their metadata.
                                               Format: {cell_name: {'definition': str, 'has_children': bool}}
                                               Returns empty dict if node_data is None or has no children.
                                               
    Raises:
        TreeValidationError: If node_data is malformed (missing 'definition' field).
    """
    if node_data is None:
        return {}
        
    children_dict = node_data.get("children", node_data)

    res_data = {}
    for key, values in children_dict.items():
        # Validate that child has required fields
        if "definition" not in values:
            raise TreeValidationError(
                f"Child node '{key}' is missing required 'definition' field"
            )
        
        res_data[key] = {
            "definition": values["definition"],
            "has_children": bool(values.get("children", {})),
        }
    if include_parent:
        res_data[parent_name] = {
            "definition": node_data["definition"],
            "has_children": True,
        }
    return res_data


def load_tree(tree_file: str) -> Dict:
    """
    Load cell type hierarchy tree from JSON file.
    
    Args:
        tree_file: Path to JSON file containing tree structure
        
    Returns:
        Dictionary containing the tree structure
    """
    with open(tree_file, 'r') as f:
        return json.load(f)


def get_tree_depth(node_data: Dict, current_depth: int = 0) -> int:
    """
    Calculate the maximum depth of a tree node.
    
    Args:
        node_data: Dictionary containing node information
        current_depth: Current depth level (default 0)
        
    Returns:
        Maximum depth of the tree from this node
    """
    if not node_data or "children" not in node_data or not node_data["children"]:
        return current_depth
    
    max_depth = current_depth
    for child_name, child_data in node_data["children"].items():
        child_depth = get_tree_depth(child_data, current_depth + 1)
        max_depth = max(max_depth, child_depth)
    
    return max_depth
