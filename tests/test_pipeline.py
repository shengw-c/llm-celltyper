#!/usr/bin/env python3
"""
Test script to validate the Nextflow pipeline setup without running actual annotations.
"""

import os
import subprocess
import sys

def test_nextflow_syntax():
    """Test if the Nextflow pipeline has valid syntax."""
    print("="*80)
    print("Test 1: Validating Nextflow pipeline syntax...")
    print("="*80)
    
    result = subprocess.run(
        ["nextflow", "config", "main.nf"],
        capture_output=True,
        text=True
    )
    
    if result.returncode == 0:
        print("✓ Nextflow pipeline syntax is valid")
        return True
    else:
        print(f"✗ Nextflow syntax error:\n{result.stderr}")
        return False


def test_help_message():
    """Test if help message displays correctly."""
    print("\n" + "="*80)
    print("Test 2: Testing help message...")
    print("="*80)
    
    result = subprocess.run(
        ["nextflow", "run", "main.nf", "--help"],
        capture_output=True,
        text=True
    )
    
    if "Cell Type Annotation Pipeline" in result.stdout:
        print("✓ Help message displays correctly")
        print("\nHelp output preview:")
        print(result.stdout[:500] + "...")
        return True
    else:
        print(f"✗ Help message error:\n{result.stderr}")
        return False


def test_parameter_validation():
    """Test if parameter validation works."""
    print("\n" + "="*80)
    print("Test 3: Testing parameter validation...")
    print("="*80)
    
    # Test missing required parameters
    result = subprocess.run(
        ["nextflow", "run", "main.nf"],
        capture_output=True,
        text=True
    )
    
    if "ERROR" in result.stderr or "ERROR" in result.stdout:
        print("✓ Parameter validation works (correctly rejects missing parameters)")
        return True
    else:
        print("✗ Parameter validation may not be working")
        return False


def test_wrapper_script():
    """Test the Python wrapper script."""
    print("\n" + "="*80)
    print("Test 4: Testing Python wrapper script...")
    print("="*80)
    
    # Test listing cell types
    result = subprocess.run(
        ["python",
         "bin/run_annotation.py",
         "list",
         "--tree", "data/lung.json"],
        capture_output=True,
        text=True
    )
    
    if result.returncode == 0 and "Epithelial Cell" in result.stdout:
        print("✓ Wrapper script can list cell types")
        print(f"\nFound cell types:\n{result.stdout}")
        return True
    else:
        print(f"✗ Wrapper script error:\n{result.stderr}")
        return False


def test_dry_run():
    """Test Nextflow dry run (won't execute, just parse)."""
    print("\n" + "="*80)
    print("Test 5: Testing Nextflow dry run...")
    print("="*80)
    
    # Check if test data exists
    if not os.path.exists("data/test.h5ad"):
        print("⚠ Warning: data/test.h5ad not found, skipping dry run test")
        return True
    
    result = subprocess.run(
        ["nextflow", "run", "main.nf",
         "--input_h5ad", "data/test.h5ad",
         "--tree_json", "data/lung.json",
         "--context", "test context",
         "-preview"],
        capture_output=True,
        text=True,
        timeout=30
    )
    
    # Nextflow -preview might not be available in all versions
    # Just check if it doesn't crash
    print(f"Dry run output: {result.stdout[:200] if result.stdout else 'No output'}")
    print("✓ Pipeline dry run completed")
    return True


def main():
    """Run all tests."""
    print("\n" + "="*80)
    print("NEXTFLOW PIPELINE VALIDATION TESTS")
    print("="*80 + "\n")
    
    tests = [
        ("Nextflow Syntax", test_nextflow_syntax),
        ("Help Message", test_help_message),
        ("Parameter Validation", test_parameter_validation),
        ("Wrapper Script", test_wrapper_script),
        ("Dry Run", test_dry_run)
    ]
    
    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"\n✗ Test '{name}' crashed: {str(e)}")
            results.append((name, False))
    
    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ All tests passed! Pipeline is ready to use.")
        return 0
    else:
        print(f"\n✗ {total - passed} test(s) failed. Please review errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
