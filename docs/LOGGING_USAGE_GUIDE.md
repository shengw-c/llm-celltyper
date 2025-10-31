# How to Use the New Logging Features

## Overview

Two logging improvements have been implemented:

1. **Timestamps in Console Output**: All console log messages now show timestamps in the format `YYYY-MM-DD HH:MM:SS`
2. **Log Level Control**: New `--log-level` parameter to control logging verbosity

## 1. Console Timestamps

### Before
```
INFO | Starting annotation for: T cell
INFO | Configuration:
INFO |   Input: data.h5ad
```

### After
```
2025-10-29 14:30:15 | INFO | Starting annotation for: T cell
2025-10-29 14:30:15 | INFO | Configuration:
2025-10-29 14:30:15 | INFO |   Input: data.h5ad
```

**Benefit**: Easier to track when events occur and correlate with log files.

## 2. Log Level Control

### Available Log Levels

- `DEBUG`: Most verbose - shows all log messages including detailed debug information
- `INFO`: Default - shows informational messages, warnings, and errors
- `WARNING`: Shows only warnings and errors
- `ERROR`: Least verbose - shows only error messages

### Usage with run_annotation.py

#### Annotate Mode

```bash
# Default (INFO level)
./bin/run_annotation.py annotate \
  --input data/test.h5ad \
  --tree data/lung.json \
  --celltype "T cell" \
  --context "Human lung tissue"

# Debug mode (verbose)
./bin/run_annotation.py annotate \
  --input data/test.h5ad \
  --tree data/lung.json \
  --celltype "T cell" \
  --context "Human lung tissue" \
  --log-level DEBUG

# Warning mode (less verbose)
./bin/run_annotation.py annotate \
  --input data/test.h5ad \
  --tree data/lung.json \
  --celltype "T cell" \
  --context "Human lung tissue" \
  --log-level WARNING

# Error mode (minimal output)
./bin/run_annotation.py annotate \
  --input data/test.h5ad \
  --tree data/lung.json \
  --celltype "T cell" \
  --context "Human lung tissue" \
  --log-level ERROR
```

#### Split Mode

```bash
# Default (INFO level)
./bin/run_annotation.py split \
  --input data/test.h5ad \
  --tree data/lung.json \
  --context "Human lung tissue"

# Debug mode
./bin/run_annotation.py split \
  --input data/test.h5ad \
  --tree data/lung.json \
  --context "Human lung tissue" \
  --log-level DEBUG
```

### Usage with analyze_logs.py

```bash
# Default (INFO level)
./bin/analyze_logs.py \
  --log-files logs/annotator_*.log \
  --output analysis.json

# Debug mode
./bin/analyze_logs.py \
  --log-files logs/annotator_*.log \
  --output analysis.json \
  --log-level DEBUG

# Warning mode
./bin/analyze_logs.py \
  --log-files logs/annotator_*.log \
  --output analysis.json \
  --log-level WARNING
```

## Output Examples

### DEBUG Level Output
Shows everything including internal processing details:

```
2025-10-29 14:30:15 | DEBUG | Loading data from test.h5ad
2025-10-29 14:30:16 | DEBUG | Data shape: (50000, 2000)
2025-10-29 14:30:16 | DEBUG | Found 2000 highly variable genes
2025-10-29 14:30:16 | INFO | Starting annotation for: T cell
2025-10-29 14:30:16 | DEBUG | Clustering with resolution 0.5
2025-10-29 14:30:17 | DEBUG | Found 8 clusters
2025-10-29 14:30:17 | INFO | Computing marker genes...
2025-10-29 14:30:18 | DEBUG | Running differential expression for cluster 0
...
```

### INFO Level Output (Default)
Shows normal operational information:

```
2025-10-29 14:30:16 | INFO | Logger initialized. Log file: logs/T_cell_20251029_143016.log
2025-10-29 14:30:16 | INFO | ================================================================================
2025-10-29 14:30:16 | INFO | Starting annotation for: T cell
2025-10-29 14:30:16 | INFO | ================================================================================
2025-10-29 14:30:16 | INFO | Configuration:
2025-10-29 14:30:16 | INFO |   Input: test.h5ad
2025-10-29 14:30:16 | INFO |   Nametag: lev0_T_cell
2025-10-29 14:30:17 | INFO | Computing marker genes...
2025-10-29 14:30:20 | INFO | Running pathway enrichment...
2025-10-29 14:30:25 | INFO | Calling LLM for annotation...
```

### WARNING Level Output
Shows only warnings and errors:

```
2025-10-29 14:30:25 | WARNING | Cluster 3 has low marker specificity (score: 0.45)
2025-10-29 14:30:30 | WARNING | LLM API rate limit approached, waiting 5s
2025-10-29 14:30:45 | ERROR | Failed to annotate cluster 5: API timeout
```

### ERROR Level Output
Shows only errors:

```
2025-10-29 14:30:45 | ERROR | Failed to annotate cluster 5: API timeout
2025-10-29 14:31:00 | ERROR | Maximum retries exceeded for cluster 5
```

## When to Use Each Level

### Use DEBUG when:
- Developing or testing new features
- Investigating issues or bugs
- Understanding the detailed execution flow
- Troubleshooting unexpected behavior
- Performance profiling

### Use INFO (default) when:
- Running normal analysis
- Monitoring progress
- General production use
- Getting overview of what's happening

### Use WARNING when:
- Running production pipelines
- Only want to see potential issues
- Reducing log noise
- Automated batch processing

### Use ERROR when:
- Only need to know about failures
- Running in automated systems
- Minimal output required
- Checking if workflow completed successfully

## File Logs vs Console Output

**Important**: Regardless of the console log level, all log files (in the `logs/` directory) always contain **full DEBUG-level** information. This ensures you never lose important debugging information.

```
Console:     Controlled by --log-level parameter
Log Files:   Always contain ALL messages (DEBUG level)
```

## Testing the New Features

Run the test script to verify the changes:

```bash
python test_logging_simple.py
```

This will demonstrate:
- Timestamps appearing in console output
- Different log levels filtering messages appropriately
- All messages being saved to log files

## Backward Compatibility

All changes are backward compatible:
- If `--log-level` is not specified, the default is `INFO` (same behavior as before)
- Existing scripts and pipelines work without modification
- Only the console format changed (timestamps added)
- Log file format unchanged

## Migration Guide

No migration required! However, you can improve your workflows:

### Before:
```bash
./bin/run_annotation.py annotate --input data.h5ad --tree tree.json --celltype "T cell" --context "lung"
```

### After (optional improvements):
```bash
# Enable debug for troubleshooting
./bin/run_annotation.py annotate --input data.h5ad --tree tree.json --celltype "T cell" --context "lung" --log-level DEBUG

# Reduce output for production
./bin/run_annotation.py annotate --input data.h5ad --tree tree.json --celltype "T cell" --context "lung" --log-level WARNING
```

## Tips and Best Practices

1. **Development**: Use `--log-level DEBUG` to see everything
2. **Production**: Use default `INFO` or `--log-level WARNING`
3. **Troubleshooting**: Temporarily enable `DEBUG` without code changes
4. **Automated Pipelines**: Use `--log-level ERROR` to reduce noise
5. **Always check log files**: Console may filter messages, but files have everything

## Example Nextflow Integration

If using with Nextflow, you can add the parameter to your process calls:

```groovy
process ANNOTATE {
    input:
    path input_file
    val cell_type
    val log_level
    
    script:
    """
    ./bin/run_annotation.py annotate \\
        --input ${input_file} \\
        --tree ${tree_file} \\
        --celltype "${cell_type}" \\
        --context "${context}" \\
        --log-level ${log_level}
    """
}
```

Then configure in `nextflow.config`:

```groovy
params {
    log_level = 'INFO'  // or 'DEBUG', 'WARNING', 'ERROR'
}
```
