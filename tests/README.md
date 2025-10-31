# Test Suite

This directory contains all test files for the LLM Cell Typer pipeline.

## Test Files

### test_modules.py
Tests for individual Python modules and library functions.

**Tests:**
- Module imports
- Logging system
- Tree utilities
- Hierarchical annotation collector

**Run:**
```bash
python tests/test_modules.py
```

### test_pipeline.py
Tests for Nextflow pipeline configuration and wrapper scripts.

**Tests:**
- Nextflow syntax validation
- Help message display
- Parameter validation
- Wrapper script functionality
- Dry run execution

**Run:**
```bash
python tests/test_pipeline.py
```

## Running All Tests

### Individual Tests
```bash
# From project root
python tests/test_modules.py
python tests/test_pipeline.py
```

### All Tests at Once
```bash
# From project root
for test in tests/test_*.py; do
    echo "Running $test..."
    python $test || exit 1
done
```

### Using pytest (if installed)
```bash
pytest tests/ -v
```

---

## Test Results Summary

**Total Tests:** 9 ✅
- test_modules.py: 4 tests
- test_pipeline.py: 5 tests

**Last Run:** October 27, 2025  
**Status:** ✅ All tests passing (9/9)

**Note:** After migrating tests to `tests/` directory, all import paths were updated to use `Path(__file__).parent.parent` to correctly reference the project root. See `TEST_MIGRATION_FIX.md` for details.

## Test Logs

Test execution logs are stored in `tests/test_logs/` directory.

## Adding New Tests

1. Create a new test file: `tests/test_<feature>.py`
2. Follow the existing test structure
3. Import required modules from `bin/lib/`
4. Update this README with test description

## Test Coverage

Current test coverage includes:
- ✅ Module imports and dependencies
- ✅ Core functionality (tree utils, annotator, etc.)
- ✅ LLM client error handling
- ✅ Data processing functions
- ✅ Nextflow configuration
- ✅ Pipeline wrapper scripts
- ✅ All refactoring improvements

**Not yet covered:**
- Integration tests with real data
- End-to-end pipeline execution
- SLURM cluster execution
- Different LLM models

---

**Note:** All tests are designed to run without external dependencies (no API keys required, no actual LLM calls).
