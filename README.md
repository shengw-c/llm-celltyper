# Cell Type Annotation Pipeline v2

Nextflow pipeline for hierarchical cell type annotation of single-cell RNA-seq data using LLM-based iterative annotation with **random cluster ID tracking**.

## Overview

This pipeline performs automated, hierarchical cell type annotation on single-cell RNA-seq data using:
- Iterative clustering with automatic resolution optimization
- **Random 6-character alphanumeric cluster IDs** for privacy and traceability
- LLM-based cell type identification
- Marker gene and pathway enrichment analysis
- Hierarchical annotation tracking across multiple levels
- Complete lineage tracking for all clusters

## New Feature: Random Cluster IDs

The pipeline now assigns **random 6-character alphanumeric IDs** (e.g., `A3X9K2`, `B7M4N1`) to clusters instead of sequential numeric IDs. This provides:

- **Privacy**: Cluster IDs are not predictable or sequential
- **Traceability**: Full lineage tracking in `cluster_id_mappings.json`
- **Uniqueness**: Each cluster gets a globally unique identifier
- **Hierarchical tracking**: Parent-child relationships are preserved

### Cluster ID Outputs

The pipeline generates additional files for cluster tracking:

```
results/
├── annotation/
│   ├── data/
│   │   ├── cluster_id_mappings.json      # Complete ID tracking with metadata
│   │   ├── cluster_id_lineage.csv        # Human-readable lineage table
│   │   └── hierarchical_annotation_complete.csv
```

**cluster_id_mappings.json** structure:
```json
{
  "A3X9K2": {
    "iteration": 1,
    "parent_id": null,
    "parent_original_id": null
  },
  "B7M4N1": {
    "iteration": 2,
    "parent_id": "A3X9K2",
    "parent_original_id": "A3X9K2"
  }
}
```

**cluster_id_lineage.csv** structure:
```csv
cluster_id,iteration,parent_id,depth,lineage
A3X9K2,1,,0,A3X9K2
B7M4N1,2,A3X9K2,1,A3X9K2 → B7M4N1
```

## Quick Start

### Prerequisites

- Nextflow (>= 21.04.0)
- Python 3.8+
- Required Python packages (install via pip or conda)
- Google API key for Gemini models (if not using mock mode)

### Installation

```bash
# Clone the repository
cd /home/sheng/projects/llm_celltyper_v2

# Set up Python environment (if needed)
# Install dependencies from requirements.txt

# Set Google API key (if using real LLM)
export GOOGLE_API_KEY="your-api-key-here"
```

### Basic Usage

```bash
# Run with test data and mock LLM mode (no API calls)
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --context "lung tissue, healthy adult samples" \
  --mock_llm \
  -profile test

# Run with real LLM annotation
nextflow run main.nf \
  --input_h5ad data/your_data.h5ad \
  --context "your biological context" \
  --batch_key donor_id \
  --integration \
  --cpus 16
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input_h5ad` | Path to QCed h5ad file with single-cell data |

### Optional Parameters

#### General Configuration
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--context` | "single-cell RNA-seq data" | Biological context description |
| `--batch_key` | null | Batch key column name in adata.obs for integration |
| `--integration` | false | Enable Harmony integration |
| `--cpus` | 16 | Number of CPUs per task |
| `--outdir` | results | Output directory |
| `--min_cells` | 500 | Minimum cells required for subtype annotation |

#### LLM Configuration
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--llm_model_general` | "gemini-2.5-flash" | LLM model for general tasks |
| `--llm_model_complicated` | "gemini-2.5-pro" | LLM model for complex tasks |
| `--llm_max_retries` | 3 | Maximum retries for LLM API calls |
| `--mock_llm` | false | **Enable mock mode (no API calls)** |

#### Clustering Configuration
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_resolution` | 1.0 | Maximum resolution for Leiden clustering |
| `--transition_cutoff` | 0.1 | Stability threshold for clustering |

#### Feature Extraction
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--top_genes` | 20 | Number of top marker genes per cluster |
| `--top_pathways` | 20 | Number of top pathways per cluster |
| `--gsea_databases` | "MSigDB_Hallmark_2020,KEGG_2021_Human" | GSEA databases |

#### Logging
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--log_level` | INFO | Logging level (DEBUG, INFO, WARNING, ERROR) |

## Execution Profiles

### Standard (Local)
```bash
nextflow run main.nf --input_h5ad data/test.h5ad -profile standard
```

### SLURM Cluster
```bash
nextflow run main.nf --input_h5ad data/test.h5ad -profile slurm
```

### Test (Mock LLM)
```bash
nextflow run main.nf --input_h5ad data/test.h5ad -profile test
```

The test profile automatically enables:
- Mock LLM mode (no API calls)
- Reduced resources (2 CPUs, 4 GB memory)
- Minimal cells threshold (100)

## Mock Mode

Mock mode allows testing the pipeline without making actual LLM API calls. This is useful for:
- Pipeline development and testing
- Validating data processing steps
- Testing on systems without API access

### Enabling Mock Mode

**Option 1: Command-line flag**
```bash
nextflow run main.nf --input_h5ad data/test.h5ad --mock_llm
```

**Option 2: Test profile**
```bash
nextflow run main.nf --input_h5ad data/test.h5ad -profile test
```

**Option 3: Environment variable**
```bash
export MOCK_LLM_MODE=1
nextflow run main.nf --input_h5ad data/test.h5ad
```

In mock mode, the LLM functions will return dummy annotations instead of calling the API.

## Output Structure

```
results/
├── annotation/
│   ├── data/
│   │   ├── iter1_cluster_data.json
│   │   ├── iter2_cluster_data.json
│   │   ├── ...
│   │   ├── cluster_id_mappings.json      # Complete ID tracking with metadata
│   │   └── cluster_id_lineage.csv        # Human-readable lineage table
│   ├── responses/
│   │   └── cell_annotation.json
│   └── figures/
│       ├── umap_iter1.png
│       └── ...
├── logs/
│   ├── annotator_*.log
│   └── data_processing_*.log
├── pipeline_trace.txt
├── pipeline_timeline.html
├── pipeline_report.html
└── pipeline_dag.html
```

### Key Output Files

- **hierarchical_annotation_complete.csv**: Final cell type annotations with hierarchical structure
- **cell_annotation.json**: LLM annotation responses for all iterations
- **iter*_cluster_data.json**: Cluster information (markers, pathways) for each iteration
- **pipeline_report.html**: Nextflow execution report

## Pipeline Workflow

1. **Data Loading**: Load h5ad file and perform QC
2. **Initial Clustering**: Leiden clustering with resolution optimization
3. **Feature Extraction**: Identify marker genes and enriched pathways
4. **LLM Annotation**: Cell type identification using LLM
5. **Iteration**: For clusters marked for splitting:
   - Re-cluster at higher resolution
   - Extract features
   - Annotate subtypes
6. **Termination**: Stop when no clusters need further subdivision
7. **Output**: Export hierarchical annotations

## Examples

### Example 1: Quick Test
```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --context "test dataset" \
  --mock_llm \
  --cpus 4 \
  -profile test
```

### Example 2: Production Run
```bash
nextflow run main.nf \
  --input_h5ad data/lung_data.h5ad \
  --context "lung tissue from healthy adult donors" \
  --batch_key donor_id \
  --integration \
  --cpus 32 \
  --min_cells 1000 \
  --max_resolution 1.5 \
  --llm_model_general "gemini-2.5-flash" \
  --outdir results/lung_annotation \
  -profile slurm
```

### Example 3: Custom GSEA Databases
```bash
nextflow run main.nf \
  --input_h5ad data/immune_cells.h5ad \
  --context "peripheral blood mononuclear cells" \
  --gsea_databases "MSigDB_Hallmark_2020,GO_Biological_Process_2025,KEGG_2021_Human" \
  --top_genes 30 \
  --top_pathways 30
```

## Troubleshooting

### Common Issues

**1. API Key Error**
```
Error: GOOGLE_API_KEY not set
```
Solution: Set the environment variable or use `--mock_llm`

**2. Memory Issues**
```
Error: Java heap space
```
Solution: Increase memory in nextflow.config or use `-profile` with more resources

**3. Clustering Fails**
```
Error: No stable resolution found
```
Solution: Adjust `--max_resolution` or `--transition_cutoff` parameters

## Architecture

This pipeline uses the **v2 iterative architecture**:
- Processes the entire dataset in a single workflow
- Uses iterative deepening for hierarchical annotation
- No tree-based splitting (unlike v1)
- Automatic determination of annotation levels based on split decisions

## Citation

If you use this pipeline, please cite [your publication here].

## License

[Your license here]

## Contact

For questions or issues, please contact [your contact info].
