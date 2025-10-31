# LLM Model & SLURM Configuration Guide

This guide covers the enhanced configuration capabilities for managing different LLM models and running the pipeline on SLURM clusters.

---

## Table of Contents
1. [LLM Model Configuration](#llm-model-configuration)
2. [SLURM Cluster Execution](#slurm-cluster-execution)
3. [Execution Profiles](#execution-profiles)
4. [Configuration Examples](#configuration-examples)
5. [Troubleshooting](#troubleshooting)

---

## LLM Model Configuration

### Supported Parameters

The pipeline now supports comprehensive LLM configuration through `nextflow.config`:

```groovy
params {
    // LLM Model Configuration
    llm_model = "gemini-2.0-flash-exp"  // Model name
    llm_temperature = 0.1                // Temperature (0.0-1.0)
    llm_max_retries = 3                  // Max retry attempts
    llm_timeout = 120                    // Timeout in seconds
}
```

### Supported Models

#### Google Gemini (Default - Currently Supported)
- `gemini-2.0-flash-exp` (default, fastest)
- `gemini-1.5-pro` (more capable)
- `gemini-1.5-flash` (balanced)

#### Future Support (Requires Adapter Implementation)
- OpenAI: `gpt-4`, `gpt-4-turbo`, `gpt-3.5-turbo`
- Anthropic: `claude-3-opus`, `claude-3-sonnet`

### Command-Line Override

You can override LLM settings when running the pipeline:

```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --llm_model "gemini-1.5-pro" \
  --llm_temperature 0.2 \
  --llm_max_retries 5 \
  --llm_timeout 180
```

### Python API

When using the Python modules directly:

```python
from lib.llm_client import CellTypeAnnotationClient

# Initialize with custom parameters
client = CellTypeAnnotationClient(
    model_name="gemini-2.0-flash-exp",
    temperature=0.1,
    max_retries=3,
    timeout=120
)

# Or via annotate_cell_types
from lib.annotator import annotate_cell_types

collector = annotate_cell_types(
    expected_cells=expected_cells,
    tree_file="data/lung.json",
    general_context="lung tissue",
    nametag="test",
    input_file="data/test.h5ad",
    llm_model="gemini-1.5-pro",
    llm_temperature=0.15,
    llm_max_retries=5,
    llm_timeout=180
)
```

### CLI Wrapper

Using the `run_annotation.py` script:

```bash
python bin/run_annotation.py annotate \
  --input data/test.h5ad \
  --tree data/lung.json \
  --celltype "Epithelial Cell" \
  --context "lung tissue from healthy adult" \
  --llm-model gemini-1.5-pro \
  --llm-temperature 0.15 \
  --llm-max-retries 5 \
  --llm-timeout 180
```

---

## SLURM Cluster Execution

### Quick Start

Run on SLURM cluster with default settings:

```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile slurm
```

### SLURM Profile Features

The `slurm` profile in `nextflow.config` provides:

- ✅ Automatic job submission to SLURM scheduler
- ✅ Dynamic resource allocation based on task attempt
- ✅ Queue selection (short, normal, long)
- ✅ Error handling with retry logic
- ✅ Rate limiting for job submission
- ✅ Proper cleanup and error reporting

### SLURM Configuration Details

```groovy
profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'normal'
            
            // Dynamic QoS based on retry attempt
            clusterOptions = { task.attempt == 1 ? '--qos=normal' : '--qos=high' }
            
            // Error handling
            errorStrategy = { task.exitStatus in [104,134,137,139,140,143,247] ? 'retry' : 'finish' }
            maxRetries = 3
            
            // Process-specific resources
            withName: 'EXTRACT_CELL_TYPES' {
                cpus = 1
                memory = { 4.GB * task.attempt }
                time = { 30.m * task.attempt }
                queue = 'short'
            }
            
            withName: 'ANNOTATE_CELL_TYPE' {
                cpus = { params.cpus }
                memory = { 16.GB * task.attempt }
                time = { 8.h * task.attempt }
                queue = { task.attempt > 1 ? 'long' : 'normal' }
                clusterOptions = "--ntasks=1 --cpus-per-task=${params.cpus}"
            }
            
            withName: 'GENERATE_SUMMARY' {
                cpus = 2
                memory = { 8.GB * task.attempt }
                time = { 1.h * task.attempt }
                queue = 'normal'
            }
        }
        
        executor {
            queueSize = 50              // Max 50 parallel jobs
            submitRateLimit = '10 sec'  // Max 1 job every 10 seconds
            pollInterval = '30 sec'     // Check job status every 30 seconds
        }
    }
}
```

### Resource Scaling

Resources automatically scale with retry attempts:

| Attempt | Memory | Time | Queue |
|---------|--------|------|-------|
| 1 | 16 GB | 8 h | normal |
| 2 | 32 GB | 16 h | long |
| 3 | 48 GB | 24 h | long |

### Customizing SLURM Settings

#### Option 1: Modify nextflow.config

Edit `nextflow.config` to customize for your cluster:

```groovy
profiles {
    slurm {
        process {
            queue = 'your_partition_name'
            clusterOptions = '--account=your_account --qos=your_qos'
        }
    }
}
```

#### Option 2: Create Custom Config

Create `custom.config`:

```groovy
process {
    executor = 'slurm'
    queue = 'compute'
    clusterOptions = '--account=mylab --qos=normal --time=12:00:00'
}
```

Run with custom config:

```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile slurm \
  -c custom.config
```

---

## Execution Profiles

### Available Profiles

| Profile | Description | Use Case |
|---------|-------------|----------|
| `standard` | Local execution (default) | Development, testing on workstation |
| `slurm` | SLURM cluster | Production runs on HPC |
| `slurm_gpu` | SLURM with GPU | Future GPU-accelerated models |
| `cluster_hpc` | High-performance cluster | Large-scale parallel execution |
| `test` | Minimal resources | Quick testing, CI/CD |
| `production` | Enhanced monitoring | Production deployments |
| `conda` | Conda environment | Isolated dependencies |

### Profile Details

#### Standard (Local)
```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json
  # No profile flag needed (default)
```

**Resources:**
- Executor: local
- Max parallel: 10 jobs
- Memory: 8-16 GB per task

#### SLURM
```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile slurm
```

**Features:**
- Auto job submission
- Dynamic resource scaling
- Queue management (short/normal/long)
- Max 50 parallel jobs

#### SLURM GPU
```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile slurm_gpu
```

**Features:**
- GPU allocation (`--gres=gpu:1`)
- Queue: gpu
- Higher memory (32 GB base)
- Faster timeout (4h base)

#### Test
```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile test
```

**Resources:**
- CPUs: 2
- Memory: 4 GB
- Time: 30m
- Min cells: 100 (faster testing)
- LLM retries: 1

#### Production
```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile production
```

**Features:**
- Max retries: 5
- LLM timeout: 180s
- Enhanced error handling
- Full tracing/reporting enabled
- Higher memory allocation

### Combining Profiles

You can combine multiple profiles:

```bash
# SLURM + Conda environment
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile slurm,conda

# Test + Production monitoring
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile test,production
```

---

## Configuration Examples

### Example 1: Basic SLURM Run

```bash
#!/bin/bash
#SBATCH --job-name=celltyper
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load nextflow

nextflow run main.nf \
  --input_h5ad /data/samples/lung_healthy.h5ad \
  --tree_json /data/references/lung_hierarchy.json \
  --context "lung tissue from healthy adult donor" \
  --batch_key donor_id \
  --integration \
  --cpus 16 \
  -profile slurm \
  -resume
```

### Example 2: Custom LLM Settings

```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --llm_model "gemini-1.5-pro" \
  --llm_temperature 0.05 \
  --llm_max_retries 10 \
  --llm_timeout 300 \
  --min_cells 500 \
  -profile slurm
```

### Example 3: High-Performance Setup

```bash
nextflow run main.nf \
  --input_h5ad data/large_dataset.h5ad \
  --tree_json data/lung.json \
  --context "lung tissue, 500k cells" \
  --cpus 32 \
  --min_cells 5000 \
  -profile cluster_hpc \
  -with-trace \
  -with-timeline \
  -with-report
```

### Example 4: Testing Configuration

```bash
# Quick test with minimal resources
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --cpus 2 \
  -profile test \
  -resume
```

---

## Troubleshooting

### SLURM Issues

#### Jobs Not Submitting

**Problem:** Pipeline starts but no jobs appear in SLURM queue

**Solutions:**
```bash
# Check SLURM configuration
sinfo  # Verify partition names
squeue -u $USER  # Check your queue

# Verify nextflow can submit jobs
nextflow run main.nf --help  # Should work

# Check partition names in nextflow.config
grep "queue" nextflow.config

# Update partition names
# Edit nextflow.config:
#   queue = 'your_partition_name'
```

#### Out of Memory Errors

**Problem:** Jobs fail with OOM (Out of Memory) errors

**Solutions:**
```bash
# Increase base memory in nextflow.config
# Edit process.withName: 'ANNOTATE_CELL_TYPE':
#   memory = { 32.GB * task.attempt }  # Increased from 16

# Or override via command line
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  -profile slurm \
  -c <(echo "process.memory = '32 GB'")
```

#### Time Limit Exceeded

**Problem:** Jobs killed due to walltime limits

**Solutions:**
```bash
# Increase time limit in nextflow.config
# Edit process.withName: 'ANNOTATE_CELL_TYPE':
#   time = { 24.h * task.attempt }

# Or use long queue
#   queue = 'long'
```

#### Account/QoS Issues

**Problem:** Jobs rejected due to account or QoS

**Solutions:**
```bash
# Check your accounts
sacctmgr show user $USER

# Update clusterOptions in nextflow.config
# Edit profiles.slurm.process:
#   clusterOptions = '--account=YOUR_ACCOUNT --qos=YOUR_QOS'
```

### LLM Issues

#### API Timeout

**Problem:** LLM calls timing out

**Solutions:**
```bash
# Increase timeout
nextflow run main.nf \
  --llm_timeout 300 \  # 5 minutes
  --llm_max_retries 10
```

#### Rate Limiting

**Problem:** Too many API calls causing rate limits

**Solutions:**
```bash
# Reduce parallel jobs
# Edit nextflow.config profiles.slurm.executor:
#   queueSize = 20  # Reduced from 50
#   submitRateLimit = '30 sec'  # Slower submission
```

#### Model Not Found

**Problem:** Specified model doesn't exist

**Solutions:**
```bash
# Verify model name (currently only Gemini supported)
nextflow run main.nf \
  --llm_model "gemini-2.0-flash-exp"  # Use exact name

# Check available models in llm_client.py
```

### General Pipeline Issues

#### Resume Not Working

**Problem:** `-resume` not resuming from cached results

**Solutions:**
```bash
# Clean work directory
nextflow clean -f

# Use explicit resume
nextflow run main.nf -resume [previous_run_name]

# Check work directory permissions
ls -la work/
```

#### No Output Generated

**Problem:** Pipeline completes but no results

**Solutions:**
```bash
# Check error strategy (should not be 'ignore' for all)
grep errorStrategy nextflow.config

# Check exit codes
find work/ -name ".exitcode" -exec cat {} \;

# Review logs
find work/ -name "*.log" -exec tail {} \;
```

---

## Advanced Configuration

### Custom Profile for Your Cluster

Create a custom profile in `nextflow.config`:

```groovy
profiles {
    mylab_cluster {
        process {
            executor = 'slurm'
            queue = 'lab_queue'
            clusterOptions = '--account=mylab --partition=compute --qos=normal'
            
            withName: 'ANNOTATE_CELL_TYPE' {
                cpus = 32
                memory = '128 GB'
                time = '48 h'
                clusterOptions = '--account=mylab --partition=himem --qos=high'
            }
        }
        
        executor {
            queueSize = 100
            submitRateLimit = '20/1min'
        }
        
        params {
            llm_max_retries = 10
            llm_timeout = 300
        }
    }
}
```

Use it:
```bash
nextflow run main.nf -profile mylab_cluster
```

### Environment Variables

Set LLM API keys via environment:

```bash
# For Google Gemini
export GOOGLE_API_KEY="your_api_key_here"

# For OpenAI (future)
export OPENAI_API_KEY="your_api_key_here"

# For Anthropic (future)
export ANTHROPIC_API_KEY="your_api_key_here"

nextflow run main.nf --input_h5ad data/test.h5ad --tree_json data/lung.json
```

---

## Summary

### Key Features

✅ **Flexible LLM Configuration**
- Support for multiple models
- Configurable temperature, retries, timeout
- Command-line and config file options

✅ **SLURM Integration**
- Multiple profiles (slurm, slurm_gpu, cluster_hpc)
- Dynamic resource scaling
- Intelligent queue selection
- Robust error handling

✅ **Easy Customization**
- Override any setting via command line
- Create custom profiles
- Combine multiple profiles
- Environment variable support

### Quick Reference

```bash
# Local run
nextflow run main.nf --input_h5ad X --tree_json Y

# SLURM run
nextflow run main.nf --input_h5ad X --tree_json Y -profile slurm

# Custom LLM
nextflow run main.nf --input_h5ad X --tree_json Y \
  --llm_model gemini-1.5-pro --llm_temperature 0.2

# Test mode
nextflow run main.nf --input_h5ad X --tree_json Y -profile test

# Production with monitoring
nextflow run main.nf --input_h5ad X --tree_json Y \
  -profile slurm,production
```

---

**Last Updated:** October 27, 2025  
**Version:** 2.0.0
