#!/usr/bin/env python3
"""
Generate an interactive HTML summary of cell type annotation results.

This module collects all LLM responses, UMAP plots, cluster tree visualizations,
log analysis, and generates a comprehensive HTML report with interactive navigation.
"""

import json
import os
import glob
import sys
from typing import Dict, List, Any, Optional
import base64
from pathlib import Path

# Add lib directory to path
sys.path.insert(0, str(Path(__file__).parent / "lib"))

from lib import PipelineLogger

logger = PipelineLogger.get_logger(__name__)


def encode_image_to_base64(image_path: str) -> str:
    """
    Encode an image file to base64 string.
    
    Args:
        image_path: Path to the image file.
    
    Returns:
        Base64-encoded image string.
    """
    if not os.path.exists(image_path):
        return ""
    
    with open(image_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def load_log_analysis(work_dir: str = "work") -> Optional[Dict]:
    """
    Load the log analysis JSON if available.
    
    Args:
        work_dir: Working directory containing log_analysis.json
    
    Returns:
        Dictionary containing log analysis or None if not found
    """
    log_analysis_file = os.path.join(work_dir, "log_analysis.json")
    
    if os.path.exists(log_analysis_file):
        logger.info(f"Loading log analysis from {log_analysis_file}")
        with open(log_analysis_file, 'r') as f:
            return json.load(f)
    else:
        logger.warning(f"Log analysis file not found: {log_analysis_file}")
        return None


def collect_annotation_results(work_dir: str = "work") -> Dict[str, List[Dict[str, Any]]]:
    """
    Collect all annotation results from the work directory.
    
    Args:
        work_dir: Path to the work directory containing responses and figures.
    
    Returns:
        Dictionary containing organized annotation results by hierarchy level.
    """
    logger.info("Collecting annotation results")
    
    responses_dir = os.path.join(work_dir, "responses")
    figures_dir = os.path.join(work_dir, "figures")
    
    if not os.path.exists(responses_dir):
        logger.warning(f"Responses directory not found: {responses_dir}")
        return {}
    
    # Find all cluster selection and annotation files
    cluster_files = sorted(glob.glob(os.path.join(responses_dir, "*_cluster_selection.json")))
    annotation_files = sorted(glob.glob(os.path.join(responses_dir, "*_cell_annotation.json")))
    
    logger.info(f"Found {len(cluster_files)} cluster selection files")
    logger.info(f"Found {len(annotation_files)} annotation files")
    
    results = {}
    
    for cluster_file in cluster_files:
        nametag = os.path.basename(cluster_file).replace("_cluster_selection.json", "")
        
        # Load cluster selection response
        with open(cluster_file, 'r') as f:
            cluster_data = json.load(f)
        
        # Find matching annotation file
        annotation_file = os.path.join(responses_dir, f"{nametag}_cell_annotation.json")
        if os.path.exists(annotation_file):
            with open(annotation_file, 'r') as f:
                annotation_data = json.load(f)
        else:
            annotation_data = {}
        
        # Collect associated images
        cluster_tree_path = cluster_data.get("cluster_tree_file", "")
        umap_path = annotation_data.get("umap_file", "")
        
        # Make paths absolute if they're relative
        # The JSON files have paths like "./figures/..." which are relative to where they were created
        # We need to resolve them from the current work_dir location
        if cluster_tree_path and not os.path.isabs(cluster_tree_path):
            # Remove leading ./ if present
            clean_path = cluster_tree_path.lstrip('./')
            # Try several possible locations
            possible_paths = [
                os.path.join(work_dir, clean_path),  # work_dir/figures/...
                os.path.join(work_dir, '..', clean_path),  # work_dir/../figures/...
                os.path.join(work_dir, '..', 'figures', os.path.basename(cluster_tree_path)),  # work_dir/../figures/basename
                os.path.join(work_dir, '..', 'figures', 'figures', os.path.basename(cluster_tree_path)),  # Handle nested figures/figures
            ]
            for abs_path in possible_paths:
                if os.path.exists(abs_path):
                    cluster_tree_path = abs_path
                    break
        
        if umap_path and not os.path.isabs(umap_path):
            # Remove leading ./ if present
            clean_path = umap_path.lstrip('./')
            # Try several possible locations
            possible_paths = [
                os.path.join(work_dir, clean_path),  # work_dir/figures/...
                os.path.join(work_dir, '..', clean_path),  # work_dir/../figures/...
                os.path.join(work_dir, '..', 'figures', os.path.basename(umap_path)),  # work_dir/../figures/basename
                os.path.join(work_dir, '..', 'figures', 'figures', os.path.basename(umap_path)),  # Handle nested figures/figures
            ]
            for abs_path in possible_paths:
                if os.path.exists(abs_path):
                    umap_path = abs_path
                    break
        
        results[nametag] = {
            "cluster_selection": cluster_data,
            "annotation": annotation_data,
            "cluster_tree_image": encode_image_to_base64(cluster_tree_path) if cluster_tree_path else "",
            "umap_image": encode_image_to_base64(umap_path) if umap_path else "",
            "cluster_tree_path": cluster_tree_path,
            "umap_path": umap_path,
            "level": cluster_data.get("level", 0)
        }
    
    logger.info(f"Collected results for {len(results)} annotation runs")
    return results


def build_tree_recursive(nametag: str, results: Dict[str, Any], all_nametags: List[str]) -> Dict[str, Any]:
    """
    Recursively build a tree node with its children.
    
    Args:
        nametag: Current node's nametag
        results: All results data
        all_nametags: List of all available nametags
    
    Returns:
        Dictionary representing the tree node with nested children
    """
    # Find immediate children (nametags that start with this nametag + "_")
    children = {}
    for other_tag in all_nametags:
        if other_tag != nametag and other_tag.startswith(nametag + "_"):
            # Check if this is an IMMEDIATE child (no intermediate parent)
            is_immediate = True
            for potential_parent in all_nametags:
                if (potential_parent != nametag and 
                    potential_parent != other_tag and
                    other_tag.startswith(potential_parent + "_") and
                    potential_parent.startswith(nametag + "_")):
                    # Found an intermediate parent
                    is_immediate = False
                    break
            
            if is_immediate:
                # Extract just the child's own name (not the full path)
                child_name = other_tag.replace(nametag + "_", "").replace("_", " ")
                children[other_tag] = build_tree_recursive(other_tag, results, all_nametags)
                children[other_tag]["display_name"] = child_name
    
    return {
        "nametag": nametag,
        "data": results.get(nametag, {}),
        "children": children
    }


def organize_results_hierarchy(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Organize results into hierarchical structure based on nametag patterns.
    
    Groups results by:
    - First round (lev0_all_types)
    - Major cell types (lev0_MajorType)
    - Subtypes under each major type (lev0_MajorType_SubType_etc)
    
    Args:
        results: Flat dictionary of annotation results
    
    Returns:
        Hierarchically organized dictionary with 'overview', 'first_round', and 'cell_types' keys
    """
    hierarchy = {
        "overview": None,
        "first_round": None,
        "cell_types": {}
    }
    
    all_nametags = list(results.keys())
    
    # Separate first round
    for nametag, data in results.items():
        if nametag == "lev0_all_types":
            hierarchy["first_round"] = {nametag: data}
            break
    
    # Find major cell types (those that have no parent except possibly lev0_all_types)
    major_types = []
    for nametag in all_nametags:
        if nametag == "lev0_all_types":
            continue
        
        # Check if this has a parent other than lev0_all_types
        has_parent = False
        for other_tag in all_nametags:
            if other_tag != nametag and other_tag != "lev0_all_types" and nametag.startswith(other_tag + "_"):
                has_parent = True
                break
        
        if not has_parent:
            major_types.append(nametag)
    
    # Build tree for each major type
    for major_nametag in sorted(major_types):
        major_type_name = "_".join(major_nametag.split("_")[1:])  # Remove "lev0_"
        hierarchy["cell_types"][major_type_name] = build_tree_recursive(major_nametag, results, all_nametags)
        hierarchy["cell_types"][major_type_name]["display_name"] = major_type_name.replace("_", " ")
    
    return hierarchy


def generate_log_analysis_section(log_analysis: Optional[Dict]) -> str:
    """
    Generate HTML section for log analysis if available.
    
    Args:
        log_analysis: Dictionary containing log analysis results
    
    Returns:
        HTML string for log analysis section
    """
    if not log_analysis:
        return ""
    
    summary = log_analysis.get('summary', {})
    runs = log_analysis.get('runs', [])
    errors = log_analysis.get('errors', [])
    warnings = log_analysis.get('warnings', [])
    recommendations = log_analysis.get('recommendations', [])
    
    html = """
                <div class="section">
                    <h2>📊 Pipeline Analysis Summary</h2>
                    <div class="summary-stats">
                        <div class="stat-card">
                            <div class="stat-value">{total_annotations}</div>
                            <div class="stat-label">Total Annotations</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{successful_runs}</div>
                            <div class="stat-label">Successful Runs</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{failed_runs}</div>
                            <div class="stat-label">Failed Runs</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-value">{total_cells}</div>
                            <div class="stat-label">Total Cells Processed</div>
                        </div>
                    </div>
""".format(
        total_annotations=summary.get('total_annotations', 0),
        successful_runs=summary.get('successful_runs', 0),
        failed_runs=summary.get('failed_runs', 0),
        total_cells=summary.get('total_cells_processed', 0)
    )
    
    if errors:
        html += """
                    <h3>⚠️ Errors</h3>
                    <div class="error-list">
"""
        for error in errors:
            html += f"""
                        <div class="error-item">
                            <strong>Message:</strong> {error.get('message', 'N/A')}<br>
                            <strong>Context:</strong> {error.get('context', 'N/A')}
                        </div>
"""
        html += """
                    </div>
"""
    
    if warnings:
        html += """
                    <h3>⚡ Warnings</h3>
                    <div class="warning-list">
"""
        for warning in warnings:
            html += f"""
                        <div class="warning-item">
                            <strong>Message:</strong> {warning.get('message', 'N/A')}<br>
                            <strong>Context:</strong> {warning.get('context', 'N/A')}
                        </div>
"""
        html += """
                    </div>
"""
    
    if recommendations:
        html += """
                    <h3>💡 Recommendations</h3>
                    <ul class="recommendation-list">
"""
        for rec in recommendations:
            html += f"""
                        <li>{rec}</li>
"""
        html += """
                    </ul>
"""
    
    html += """
                </div>
"""
    
    return html


def render_tree_node(node: Dict[str, Any], depth: int = 0, first_active_ref: List[bool] = None) -> str:
    """
    Recursively render a tree node and its children as HTML buttons.
    
    Args:
        node: Tree node dictionary with 'nametag', 'display_name', and 'children'
        depth: Current depth in the tree (0 for major types)
        first_active_ref: Mutable list containing a single boolean for tracking first active state
    
    Returns:
        HTML string for the node and all its children
    """
    if first_active_ref is None:
        first_active_ref = [True]  # Will be set on first use
    
    html = ""
    nametag = node["nametag"]
    display_name = node.get("display_name", nametag)
    
    # Determine button class and styling based on depth
    if depth == 0:
        # Major cell type - use sidebar-item
        # Only set active if first_active_ref is True
        active_class = "active" if first_active_ref[0] else ""
        html += f'                <button class="sidebar-item {active_class}" onclick="openTab(event, \'{nametag}\')">↳ {display_name}</button>\n'
        if first_active_ref[0]:
            first_active_ref[0] = False  # Mark that we've used the first active
    else:
        # Child node - use sidebar-subitem with indentation
        padding = 20 + depth * 15
        arrow = "↳ " + "  " * (depth - 1)
        html += f'                <button class="sidebar-subitem" style="padding-left: {padding}px;" onclick="openTab(event, \'{nametag}\')">{arrow}{display_name}</button>\n'
    
    # Recursively render children
    if node.get("children"):
        for child_nametag in sorted(node["children"].keys()):
            child_node = node["children"][child_nametag]
            html += render_tree_node(child_node, depth + 1, first_active_ref)
    
    return html


def generate_html_summary(
    results: Dict[str, Any], 
    log_analysis: Optional[Dict],
    output_file: str = "work/annotation_summary.html"
) -> None:
    """
    Generate an interactive HTML summary of all annotation results.
    
    Args:
        results: Annotation results from collect_annotation_results.
        log_analysis: Log analysis results (optional).
        output_file: Path to save the HTML file.
    """
    logger.info(f"Generating HTML summary: {output_file}")
    
    # Start HTML with improved styling and floating sidebar
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cell Type Annotation Summary</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
            padding-left: 280px; /* Make room for sidebar */
        }
        
        /* Floating Sidebar */
        .sidebar {
            position: fixed;
            left: 20px;
            top: 20px;
            width: 250px;
            max-height: calc(100vh - 40px);
            background: white;
            border-radius: 12px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow-y: auto;
            z-index: 1000;
        }
        
        .sidebar-header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 12px 12px 0 0;
            text-align: center;
        }
        
        .sidebar-header h2 {
            font-size: 1.3em;
            margin-bottom: 5px;
        }
        
        .sidebar-header p {
            font-size: 0.85em;
            opacity: 0.9;
        }
        
        .sidebar-nav {
            padding: 10px;
        }
        
        .sidebar-item {
            padding: 12px 15px;
            cursor: pointer;
            background: transparent;
            border: none;
            font-size: 0.95em;
            font-weight: 500;
            color: #6c757d;
            border-left: 3px solid transparent;
            transition: all 0.3s;
            width: 100%;
            text-align: left;
            border-radius: 6px;
            margin: 3px 0;
            display: block;
        }
        
        .sidebar-item:hover {
            color: #495057;
            background: rgba(102, 126, 234, 0.1);
            border-left-color: rgba(102, 126, 234, 0.3);
        }
        
        .sidebar-item.active {
            color: #667eea;
            background: rgba(102, 126, 234, 0.15);
            border-left-color: #667eea;
            font-weight: 600;
        }
        
        .sidebar-section {
            margin: 10px 0;
        }
        
        .sidebar-section-header {
            padding: 10px 15px;
            font-weight: 600;
            color: #495057;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            background: rgba(102, 126, 234, 0.05);
            border-left: 3px solid #667eea;
            margin: 5px 0;
        }
        
        .sidebar-subitem {
            padding: 10px 15px 10px 30px;
            cursor: pointer;
            background: transparent;
            border: none;
            font-size: 0.9em;
            color: #6c757d;
            border-left: 3px solid transparent;
            transition: all 0.3s;
            width: 100%;
            text-align: left;
            border-radius: 6px;
            margin: 2px 0;
            display: block;
        }
        
        .sidebar-subitem:hover {
            color: #495057;
            background: rgba(102, 126, 234, 0.08);
            border-left-color: rgba(102, 126, 234, 0.3);
        }
        
        .sidebar-subitem.active {
            color: #667eea;
            background: rgba(102, 126, 234, 0.12);
            border-left-color: #667eea;
            font-weight: 600;
        }
        
        /* Scrollbar styling for sidebar */
        .sidebar::-webkit-scrollbar {
            width: 6px;
        }
        
        .sidebar::-webkit-scrollbar-track {
            background: #f1f1f1;
            border-radius: 0 12px 12px 0;
        }
        
        .sidebar::-webkit-scrollbar-thumb {
            background: #667eea;
            border-radius: 3px;
        }
        
        .sidebar::-webkit-scrollbar-thumb:hover {
            background: #764ba2;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        .header p {
            font-size: 1.1em;
            opacity: 0.9;
        }
        
        .content {
            padding: 30px;
        }
        
        .tab-content {
            display: none;
        }
        
        .tab-content.active {
            display: block;
            animation: fadeIn 0.3s;
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .section {
            margin-bottom: 30px;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }
        
        .section h2 {
            color: #667eea;
            margin-bottom: 15px;
            font-size: 1.5em;
        }
        
        .section h3 {
            color: #495057;
            margin: 15px 0 10px 0;
            font-size: 1.2em;
        }
        
        .info-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin: 15px 0;
        }
        
        .info-item {
            background: white;
            padding: 12px;
            border-radius: 6px;
            border-left: 3px solid #667eea;
        }
        
        .info-label {
            font-weight: 600;
            color: #6c757d;
            font-size: 0.85em;
            text-transform: uppercase;
            margin-bottom: 5px;
        }
        
        .info-value {
            color: #495057;
            font-size: 1.1em;
        }
        
        .images-row {
            display: flex;
            gap: 20px;
            margin: 20px 0;
            flex-wrap: wrap;
        }
        
        .images-row .image-container {
            flex: 1;
            min-width: 400px;
            margin: 0;
        }
        
        .image-container {
            margin: 20px 0;
            text-align: center;
        }
        
        .image-container img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            cursor: pointer;
            transition: transform 0.3s;
        }
        
        .image-container img:hover {
            transform: scale(1.02);
        }
        
        .annotations-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        
        .annotations-table thead {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }
        
        .annotations-table th {
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }
        
        .annotations-table td {
            padding: 12px 15px;
            border-bottom: 1px solid #dee2e6;
        }
        
        .annotations-table tr:hover {
            background: #f8f9fa;
        }
        
        .confidence-badge {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
        }
        
        .confidence-high {
            background: #d4edda;
            color: #155724;
        }
        
        .confidence-medium {
            background: #fff3cd;
            color: #856404;
        }
        
        .confidence-low {
            background: #f8d7da;
            color: #721c24;
        }
        
        .marker-list {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
        }
        
        .marker-tag {
            background: #e7f3ff;
            color: #004085;
            padding: 4px 10px;
            border-radius: 4px;
            font-size: 0.85em;
            font-family: 'Courier New', monospace;
        }
        
        .justification {
            line-height: 1.6;
            color: #495057;
            padding: 15px;
            background: white;
            border-left: 4px solid #667eea;
            border-radius: 4px;
            margin: 10px 0;
        }
        
        .summary-stats {
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            margin: 20px 0;
        }
        
        .stat-card {
            flex: 1;
            min-width: 200px;
            padding: 20px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            text-align: center;
        }
        
        .stat-value {
            font-size: 2.5em;
            font-weight: 700;
            color: #667eea;
            margin-bottom: 5px;
        }
        
        .stat-label {
            color: #6c757d;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        
        .error-list, .warning-list {
            margin: 15px 0;
        }
        
        .error-item, .warning-item {
            background: white;
            padding: 15px;
            margin: 10px 0;
            border-radius: 6px;
            border-left: 4px solid #dc3545;
        }
        
        .warning-item {
            border-left-color: #ffc107;
        }
        
        .recommendation-list {
            background: white;
            padding: 20px 40px;
            border-radius: 6px;
            line-height: 1.8;
        }
        
        .recommendation-list li {
            margin: 10px 0;
        }
        
        .modal {
            display: none;
            position: fixed;
            z-index: 2000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background: rgba(0,0,0,0.9);
            cursor: pointer;
        }
        
        .modal img {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            max-width: 95%;
            max-height: 95%;
        }
        
        .close {
            position: absolute;
            top: 20px;
            right: 40px;
            color: white;
            font-size: 40px;
            font-weight: bold;
            cursor: pointer;
        }
        
        /* Responsive design */
        @media (max-width: 768px) {
            body {
                padding-left: 20px;
            }
            
            .sidebar {
                left: -250px;
                transition: left 0.3s;
            }
            
            .sidebar.open {
                left: 20px;
            }
            
            .sidebar-toggle {
                position: fixed;
                left: 20px;
                top: 20px;
                z-index: 999;
                background: white;
                border: none;
                padding: 10px 15px;
                border-radius: 8px;
                box-shadow: 0 4px 12px rgba(0,0,0,0.2);
                cursor: pointer;
                font-size: 1.2em;
            }
        }
    </style>
</head>
<body>
    <!-- Floating Sidebar -->
    <div class="sidebar" id="sidebar">
        <div class="sidebar-header">
            <h2>🔬 Navigation</h2>
            <p>Annotation Results</p>
        </div>
        <div class="sidebar-nav" id="sidebarNav">
"""
    
    # Organize results hierarchically
    hierarchy = organize_results_hierarchy(results)
    
    # Track whether we've set the first active tab
    # True = first active tab not yet set, False = already set
    first_tab_active = [True]  # Use list to make it mutable across function calls
    
    # Add overview tab if log analysis exists
    if log_analysis:
        html_content += '            <button class="sidebar-item active" onclick="openTab(event, \'overview\')">📊 Overview</button>\n'
        first_tab_active[0] = False  # Mark as used
    
    # Add first round annotation if it exists
    if hierarchy["first_round"]:
        html_content += '            <div class="sidebar-section">\n'
        html_content += '                <div class="sidebar-section-header">🔬 FIRST ROUND</div>\n'
        for nametag in hierarchy["first_round"].keys():
            active_class = "active" if first_tab_active[0] else ""
            html_content += f'                <button class="sidebar-item {active_class}" onclick="openTab(event, \'{nametag}\')">All Cell Types</button>\n'
            if first_tab_active[0]:
                first_tab_active[0] = False  # Mark as used
        html_content += '            </div>\n'
    
    # Generate hierarchical sidebar items for each major cell type
    for major_type, cell_data in sorted(hierarchy["cell_types"].items()):
        html_content += '            <div class="sidebar-section">\n'
        # Clean up major type name for display
        display_name = cell_data.get("display_name", major_type.replace('_', ' '))
        html_content += f'                <div class="sidebar-section-header">🧬 {display_name.upper()}</div>\n'
        
        # Recursively render the entire tree for this major type
        html_content += render_tree_node(cell_data, depth=0, first_active_ref=first_tab_active)
        
        html_content += '            </div>\n'
    
    html_content += """        </div>
    </div>
    
    <div class="container">
        <div class="header">
            <h1>🔬 Cell Type Annotation Summary</h1>
            <p>Hierarchical LLM-based Cell Type Annotation Results</p>
        </div>
        
        <div class="content">
"""
    
    # Add overview tab content if log analysis exists
    if log_analysis:
        html_content += """            <div id="overview" class="tab-content active">
"""
        html_content += generate_log_analysis_section(log_analysis)
        html_content += """            </div>
"""
        first_content_active = False
    else:
        first_content_active = True
    
    # Generate content for each annotation tab
    for i, (nametag, data) in enumerate(sorted(results.items())):
        active_class = "active" if first_content_active and i == 0 else ""
        cluster_data = data.get("cluster_selection", {})
        annotation_data = data.get("annotation", {})
        
        cluster_response = cluster_data.get("cluster_selection_response", {})
        annotations = annotation_data.get("annotation_response", [])
        level = data.get("level", 0)
        
        html_content += f"""            <div id="{nametag}" class="tab-content {active_class}">
                <div class="summary-stats">
                    <div class="stat-card">
                        <div class="stat-value">{level}</div>
                        <div class="stat-label">Hierarchy Level</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">{cluster_response.get('resolution', 'N/A')}</div>
                        <div class="stat-label">Clustering Resolution</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">{len(annotations)}</div>
                        <div class="stat-label">Clusters Annotated</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">{len(cluster_data.get('expected_cell_types', []))}</div>
                        <div class="stat-label">Expected Cell Types</div>
                    </div>
                </div>
                
                <div class="section">
                    <h2>📊 Cluster Selection</h2>
                    <div class="info-grid">
                        <div class="info-item">
                            <div class="info-label">Timestamp</div>
                            <div class="info-value">{cluster_data.get('timestamp', 'N/A')}</div>
                        </div>
                        <div class="info-item">
                            <div class="info-label">Resolution</div>
                            <div class="info-value">{cluster_response.get('resolution', 'N/A')}</div>
                        </div>
                        <div class="info-item">
                            <div class="info-label">Level Color</div>
                            <div class="info-value">{cluster_response.get('level_color', 'N/A')}</div>
                        </div>
                    </div>
                    <div class="justification">
                        <strong>Justification:</strong> {cluster_response.get('justification', 'N/A')}
                    </div>
"""
        
        # Add cluster tree and UMAP images side by side
        if data.get("cluster_tree_image") or data.get("umap_image"):
            html_content += """                    <div class="images-row">
"""
            if data.get("cluster_tree_image"):
                html_content += f"""                        <div class="image-container">
                            <h3>Cluster Tree</h3>
                            <img src="data:image/png;base64,{data['cluster_tree_image']}" 
                                 alt="Cluster Tree" 
                                 onclick="openModal(this.src)">
                        </div>
"""
            if data.get("umap_image"):
                html_content += f"""                        <div class="image-container">
                            <h3>UMAP Visualization</h3>
                            <img src="data:image/png;base64,{data['umap_image']}" 
                                 alt="UMAP" 
                                 onclick="openModal(this.src)">
                        </div>
"""
            html_content += """                    </div>
"""
        
        html_content += """                </div>
                
                <div class="section">
                    <h2>🏷️ Cell Type Annotations</h2>
                    <div class="info-grid">
                        <div class="info-item">
                            <div class="info-label">Context</div>
                            <div class="info-value">""" + annotation_data.get('general_context', 'N/A') + """</div>
                        </div>
                        <div class="info-item">
                            <div class="info-label">Timestamp</div>
                            <div class="info-value">""" + annotation_data.get('timestamp', 'N/A') + """</div>
                        </div>
                    </div>
"""
        
        # Add annotations table
        if annotations:
            html_content += """                    <table class="annotations-table">
                        <thead>
                            <tr>
                                <th>Cluster ID</th>
                                <th>Cell Type</th>
                                <th>Confidence</th>
                                <th>Key Markers</th>
                                <th>Justification</th>
                            </tr>
                        </thead>
                        <tbody>
"""
            
            for ann in annotations:
                confidence = ann.get('confidence', 'Unknown')
                confidence_class = f"confidence-{confidence.lower()}"
                markers = ann.get('key_markers_cited', [])
                markers_html = ''.join([f'<span class="marker-tag">{m}</span>' for m in markers])
                
                html_content += f"""                            <tr>
                                <td><strong>{ann.get('cluster_id', 'N/A')}</strong></td>
                                <td>{ann.get('cell_type', 'N/A')}</td>
                                <td><span class="confidence-badge {confidence_class}">{confidence}</span></td>
                                <td><div class="marker-list">{markers_html}</div></td>
                                <td>{ann.get('justification', 'N/A')}</td>
                            </tr>
"""
            
            html_content += """                        </tbody>
                    </table>
"""
        
        html_content += """                </div>
            </div>
"""
    
    html_content += """        </div>
    </div>
    
    <div id="imageModal" class="modal" onclick="closeModal()">
        <span class="close">&times;</span>
        <img id="modalImage" src="">
    </div>
    
    <script>
        function openTab(evt, tabName) {
            // Hide all tab contents
            const tabContents = document.getElementsByClassName("tab-content");
            for (let i = 0; i < tabContents.length; i++) {
                tabContents[i].classList.remove("active");
            }
            
            // Remove active class from all sidebar items
            const sidebarItems = document.getElementsByClassName("sidebar-item");
            for (let i = 0; i < sidebarItems.length; i++) {
                sidebarItems[i].classList.remove("active");
            }
            
            // Remove active class from all sidebar subitems
            const sidebarSubitems = document.getElementsByClassName("sidebar-subitem");
            for (let i = 0; i < sidebarSubitems.length; i++) {
                sidebarSubitems[i].classList.remove("active");
            }
            
            // Show the selected tab content
            document.getElementById(tabName).classList.add("active");
            
            // Mark the clicked sidebar item as active
            evt.currentTarget.classList.add("active");
            
            // Scroll to top of content
            window.scrollTo({top: 0, behavior: 'smooth'});
        }
        
        function openModal(src) {
            document.getElementById("imageModal").style.display = "block";
            document.getElementById("modalImage").src = src;
        }
        
        function closeModal() {
            document.getElementById("imageModal").style.display = "none";
        }
        
        // Close modal with Escape key
        document.addEventListener('keydown', function(event) {
            if (event.key === 'Escape') {
                closeModal();
            }
        });
        
        // Mobile sidebar toggle (optional)
        function toggleSidebar() {
            const sidebar = document.getElementById('sidebar');
            sidebar.classList.toggle('open');
        }
    </script>
</body>
</html>
"""
    
    # Write HTML file
    output_dir = os.path.dirname(output_file)
    if output_dir:  # Only create directory if path has a directory component
        os.makedirs(output_dir, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logger.info(f"✓ HTML summary generated: {output_file}")


def create_annotation_summary(work_dir: str = "work", output_file: str = "work/annotation_summary.html") -> str:
    """
    Main function to create a complete annotation summary HTML report.
    
    Args:
        work_dir: Path to the work directory containing responses and figures.
        output_file: Path to save the HTML file.
    
    Returns:
        Path to the generated HTML file.
    """
    logger.info("Creating annotation summary")
    
    # Collect annotation results
    results = collect_annotation_results(work_dir)
    
    if not results:
        logger.warning("No annotation results found!")
        return ""
    
    logger.info(f"Found {len(results)} annotation runs")
    
    # Load log analysis if available
    log_analysis = load_log_analysis(work_dir)
    
    # Generate HTML summary
    generate_html_summary(results, log_analysis, output_file)
    
    return output_file


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate HTML summary of annotations")
    parser.add_argument('--work-dir', default='work', help='Work directory (default: work)')
    parser.add_argument('--output', default='work/annotation_summary.html', 
                       help='Output HTML file (default: work/annotation_summary.html)')
    
    args = parser.parse_args()
    
    html_file = create_annotation_summary(args.work_dir, args.output)
    if html_file:
        logger.info(f"✓ Summary report generated successfully!")
        logger.info(f"  Open in browser: file://{os.path.abspath(html_file)}")
