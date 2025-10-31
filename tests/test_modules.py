#!/usr/bin/env python3
"""
Test script to verify the modularized pipeline structure.
"""

import sys
import os
from pathlib import Path

# Add bin directory to path so we can import from bin.lib
# Since we're in tests/, we need to go up one level to find bin/
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

def test_imports():
    """Test that all modules can be imported."""
    print("Testing module imports...")
    
    try:
        from lib.logger import PipelineLogger
        print("  ✓ logger module")
        
        from lib.tree_utils import find_node_data, get_immediate_children, load_tree
        print("  ✓ tree_utils module")
        
        from lib.data_processing import prepare_subset_dataset, prepare_celltype_inputs
        print("  ✓ data_processing module")
        
        from lib.marker_genes import get_top_marker_genes
        print("  ✓ marker_genes module")
        
        from lib.pathway_enrichment import get_top_enriched_pathways
        print("  ✓ pathway_enrichment module")
        
        from lib.cluster_analysis import get_cluster_adjacency
        print("  ✓ cluster_analysis module")
        
        from lib.llm_client import CellTypeAnnotationClient
        print("  ✓ llm_client module")
        
        from lib.annotator import annotate_cell_types, HierarchicalAnnotation
        print("  ✓ annotator module")
        
        from lib.prompts import cluster_PROMPT, Celltyper_Instruction
        print("  ✓ prompts module")
        
        print("\n✓ All module imports successful!")
        return True
        
    except Exception as e:
        print(f"\n✗ Import failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def test_logger():
    """Test the logging system."""
    print("\nTesting logging system...")
    
    try:
        from lib.logger import PipelineLogger
        
        # Create a test logger
        logger = PipelineLogger.get_logger("test", log_dir="test_logs", log_file="test.log")
        logger.info("Test info message")
        logger.warning("Test warning message")
        logger.debug("Test debug message")
        
        # Check if log file was created
        if os.path.exists("test_logs/test.log"):
            print("  ✓ Log file created")
            
            # Read and verify log content
            with open("test_logs/test.log", 'r') as f:
                content = f.read()
                if "Test info message" in content and "Test warning message" in content:
                    print("  ✓ Log messages written correctly")
                else:
                    print("  ✗ Log messages not found in file")
                    return False
            
            # Clean up
            import shutil
            if os.path.exists("test_logs"):
                shutil.rmtree("test_logs", ignore_errors=True)
            print("  ✓ Test logs cleaned up")
            
        else:
            print("  ✗ Log file not created")
            return False
        
        print("\n✓ Logging system test passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Logging test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def test_tree_utils():
    """Test tree utility functions."""
    print("\nTesting tree utility functions...")
    
    try:
        from lib.tree_utils import find_node_data, get_immediate_children
        
        # Create a test tree
        test_tree = {
            "Cell Type A": {
                "definition": "Test cell type A",
                "children": {
                    "Cell Type A1": {
                        "definition": "Test cell type A1",
                        "children": {}
                    },
                    "Cell Type A2": {
                        "definition": "Test cell type A2",
                        "children": {}
                    }
                }
            },
            "Cell Type B": {
                "definition": "Test cell type B",
                "children": {}
            }
        }
        
        # Test find_node_data
        node_a = find_node_data(test_tree, "Cell Type A")
        if node_a and node_a["definition"] == "Test cell type A":
            print("  ✓ find_node_data works for root level")
        else:
            print("  ✗ find_node_data failed for root level")
            return False
        
        node_a1 = find_node_data(test_tree, "Cell Type A1")
        if node_a1 and node_a1["definition"] == "Test cell type A1":
            print("  ✓ find_node_data works for nested nodes")
        else:
            print("  ✗ find_node_data failed for nested nodes")
            return False
        
        # Test get_immediate_children
        children = get_immediate_children(node_a)
        if len(children) == 2 and "Cell Type A1" in children and "Cell Type A2" in children:
            print("  ✓ get_immediate_children works correctly")
        else:
            print("  ✗ get_immediate_children failed")
            return False
        
        print("\n✓ Tree utility functions test passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Tree utils test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def test_hierarchical_collector():
    """Test the HierarchicalAnnotation collector."""
    print("\nTesting HierarchicalAnnotation collector...")
    
    try:
        from lib.annotator import HierarchicalAnnotation
        import pandas as pd
        import tempfile
        import os
        
        # Create a temporary directory for testing
        with tempfile.TemporaryDirectory() as tmpdir:
            collector = HierarchicalAnnotation(output_dir=tmpdir)
            
            # Add some test annotations with proper cell IDs as index
            df1 = pd.DataFrame({
                'ann': ['Cell Type A', 'Cell Type B', 'Cell Type C']
            }, index=['cell_1', 'cell_2', 'cell_3'])
            collector.add_annotation("test_lev0", df1, level=0)
            
            # Level 1 annotations for a subset of cells
            df2 = pd.DataFrame({
                'ann': ['Cell Type A1', 'Cell Type A2']
            }, index=['cell_1', 'cell_2'])
            collector.add_annotation("test_lev1", df2, level=1)
            
            # Test export
            output_file = collector.export_final_annotations(
                output_file=os.path.join(tmpdir, "final_annotations.tsv")
            )
            
            if os.path.exists(output_file):
                # Read and verify the exported file
                exported_df = pd.read_csv(output_file, sep='\t')
                
                # Check for level columns (ann_level0, ann_level1, etc.)
                level_cols = [col for col in exported_df.columns if col.startswith('ann_level')]
                if len(level_cols) >= 2:
                    print(f"  ✓ Annotation level columns present: {', '.join(level_cols)}")
                else:
                    print(f"  ✗ Expected at least 2 level columns, got {len(level_cols)}")
                    return False
                
                # Check for ann_finest column
                if 'ann_finest' in exported_df.columns:
                    print("  ✓ ann_finest column present")
                else:
                    print("  ✗ ann_finest column missing")
                    return False
                
                # Should have 3 cells total (from level 0)
                if len(exported_df) == 3:
                    print("  ✓ All cells exported correctly")
                else:
                    print(f"  ✗ Expected 3 cells, got {len(exported_df)}")
                    return False
                
                print("  ✓ Final annotations file created and verified")
            else:
                print("  ✗ Final annotations file not created")
                return False
        
        print("\n✓ HierarchicalAnnotation collector test passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Hierarchical collector test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("=" * 80)
    print("Testing Modularized Pipeline Structure")
    print("=" * 80)
    
    tests = [
        ("Module Imports", test_imports),
        ("Logging System", test_logger),
        ("Tree Utilities", test_tree_utils),
        ("Hierarchical Collector", test_hierarchical_collector),
    ]
    
    results = {}
    for test_name, test_func in tests:
        results[test_name] = test_func()
    
    print("\n" + "=" * 80)
    print("Test Summary")
    print("=" * 80)
    
    for test_name, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{test_name:.<50} {status}")
    
    print("=" * 80)
    
    all_passed = all(results.values())
    if all_passed:
        print("\n✓ All tests passed!")
        return 0
    else:
        print("\n✗ Some tests failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
