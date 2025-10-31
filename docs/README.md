# LLM Cell Typer - Documentation

**Hierarchical, LLM-based cell type annotation for single-cell RNA-seq data**

---

## Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Architecture](#architecture)
4. [Pipeline Components](#pipeline-components)
5. [Configuration](#configuration)
6. [Output Files](#output-files)
7. [Advanced Usage](#advanced-usage)
8. [Troubleshooting](#troubleshooting)
9. [API Reference](#api-reference)

---

## Overview

This pipeline performs automated, hierarchical cell type annotation on single-cell RNA-seq data using Large Language Models (Gemini). It combines computational biology techniques with AI-powered decision making to provide accurate, reproducible cell type annotations.

### Key Features

- **Automated Annotation**: LLM-based cell type assignment with confidence scores
- **Hierarchical Processing**: Recursive annotation from major cell types to subtypes
- **Parallel Execution**: Process multiple cell types simultaneously via Nextflow
- **Interactive Reports**: Beautiful HTML summaries with visualizations
- **Complete Audit Trail**: All LLM decisions saved to JSON
- **Resource Control**: Configurable CPU and memory usage
- **Error Resilience**: Continues processing even if individual annotations fail

### Workflow

```
Input h5ad → Extract Major Cell Types → Annotate in Parallel → Recursive Subtyping → HTML Report
```

---

## Quick Start

### Prerequisites

```bash
# Required
python >= 3.8
nextflow >= 21.04.0
scanpy, pandas, numpy, etc. (see requirements.txt)

# Environment
export GOOGLE_API_KEY="your-gemini-api-key"
```

### Installation

```bash
# Clone repository
git clone <repo-url>
cd llm_celltyper

# Install dependencies
pip install -r requirements.txt

# Verify installation
python test_modules.py
python test_pipeline.py
```

### Basic Usage

```bash
# Run complete pipeline
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --context "lung tissue from healthy adult" \
  --cpus 16

# View results
firefox results/annotation_summary.html
```

### Python API

```python
from bin.lib import annotate_cell_types

collector = annotate_cell_types(
    expected_cells={"Epithelial Cell": {"has_children": True, "definition": "..."}},
    tree_file="data/lung.json",
    general_context="lung tissue",
    nametag="lev0",
    input_file="data/test.h5ad",
    batch_key=["donor_id"],
    integration=True,
    cpus_per_task=16
)

# Export results
collector.export_final_annotations("final_annotations.tsv")
```

---

## Architecture

### Project Structure

```
llm_celltyper/
├── main.nf                      # Nextflow pipeline
├── nextflow.config              # Pipeline configuration
├── bin/
│   ├── run_annotation.py        # CLI wrapper
│   ├── generate_summary.py      # HTML report generator
│   ├── analyze_logs.py          # LLM-based log analysis
│   └── lib/                     # Core modules
│       ├── __init__.py
│       ├── annotator.py         # Main annotation logic
│       ├── data_processing.py   # Data preprocessing
│       ├── marker_genes.py      # Marker gene analysis
│       ├── pathway_enrichment.py
│       ├── cluster_analysis.py
│       ├── llm_client.py        # Gemini API wrapper
│       ├── tree_utils.py        # Tree operations
│       ├── logger.py            # Logging system
│       └── prompts.py           # LLM prompts
├── data/
│   ├── test.h5ad               # Example dataset
│   └── lung.json               # Cell type hierarchy
├── docs/                        # Documentation
├── work/                        # Temporary outputs
│   ├── responses/              # LLM response JSONs
│   ├── figures/                # Visualizations
│   ├── data/                   # Processed h5ad files
│   └── logs/                   # Execution logs
└── results/                     # Final outputs
    └── annotation_summary.html
```

### Module Responsibilities

| Module | Purpose |
|--------|---------|
| `annotator.py` | Hierarchical annotation orchestration |
| `data_processing.py` | Normalization, clustering, UMAP |
| `marker_genes.py` | Differential expression analysis |
| `pathway_enrichment.py` | GSEA pathway analysis |
| `cluster_analysis.py` | PAGA cluster adjacency |
| `llm_client.py` | Gemini API communication |
| `tree_utils.py` | Cell type hierarchy navigation |
| `logger.py` | Centralized logging |
| `prompts.py` | LLM prompt templates |

---

## Pipeline Components

### 1. Data Preprocessing (`prepare_subset_dataset`)

**Input**: Raw h5ad file with counts layer  
**Output**: Normalized, clustered AnnData object

**Steps**:
1. Load and subset data by cell IDs (optional)
2. Normalize to total counts, log1p transform
3. Identify 2000 highly variable genes
4. PCA dimensionality reduction (50 components)
5. Optional Harmony integration for batch correction
6. UMAP embedding
7. Leiden clustering at 20 resolutions (0.0 - 0.95)
8. Generate cluster tree visualization

**Key Parameters**:
- `batch_key`: List of batch variables for correction
- `integration`: Enable Harmony integration
- `cpus`: Parallel processing threads

### 2. Cluster Resolution Selection (`select_cluster_resolution`)

**Input**: Cluster tree image, expected number of cell types  
**Output**: Optimal resolution value

**Process**:
1. LLM analyzes cluster tree stability (plateaus)
2. Selects resolution matching expected cell type count
3. Prioritizes stability over exact match
4. Returns resolution with justification

### 3. Feature Extraction (`prepare_celltype_inputs`)

**Input**: Clustered AnnData, selected resolution  
**Output**: Marker genes, pathways, UMAP, adjacency

**Features Computed**:
- **Marker Genes**: Top 20 by fold change and specificity (pct_diff)
- **Pathways**: Top 20 enriched pathways via GSEA
- **UMAP**: Visualization colored by clusters
- **Adjacency**: PAGA cluster neighbors

### 4. Cell Type Annotation (`annotate_cell_types`)

**Input**: Features + candidate cell types  
**Output**: Cluster → cell type mappings

**LLM Decision Process**:
1. Evaluates marker genes (expression + specificity)
2. Considers pathway enrichment
3. Checks spatial relationships (UMAP proximity)
4. Assigns cell type with confidence (High/Medium/Low)
5. Provides scientific justification

### 5. Hierarchical Recursion

**Process**:
1. For each annotated cell type with children:
   - Check cell count threshold (default: 1000 cells)
   - Extract cell IDs for that type
   - Create new subset h5ad
   - Recursively annotate children
2. Collect all levels in `HierarchicalAnnotation`
3. Export combined TSV with level indicators

**Example Hierarchy**:
```
Level 0: Epithelial Cell (5000 cells)
  ├─ Level 1: Airway Epithelial Cell (3000 cells)
  │   ├─ Level 2: Ciliated Cell (1500 cells)
  │   └─ Level 2: Secretory Cell (1500 cells)
  └─ Level 1: Alveolar Epithelial Cell (2000 cells)
      ├─ Level 2: AT1 Cell (1200 cells)
      └─ Level 2: AT2 Cell (800 cells)
```

---

## Configuration

### Nextflow Parameters

```bash
# Required
--input_h5ad FILE       # QCed h5ad with 'counts' layer
--tree_json FILE        # Cell type hierarchy JSON

# Optional
--context STRING        # Biological context (default: "single-cell RNA-seq data")
--batch_key STRING      # Batch variable name
--integration           # Enable Harmony integration
--cpus INT              # CPUs per task (default: 16)
--outdir DIR            # Output directory (default: results)
```

### Python API Parameters

#### `annotate_cell_types()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expected_cells` | Dict | Required | Cell types at this level with definitions |
| `tree_file` | str | Required | Path to hierarchy JSON |
| `general_context` | str | Required | Biological context |
| `nametag` | str | Required | Unique identifier for this run |
| `input_file` | str | Required | Path to h5ad file |
| `batch_key` | List[str] | None | Batch correction variables |
| `integration` | bool | False | Enable Harmony integration |
| `cpus_per_task` | int | 16 | Number of CPUs |
| `current_level` | int | 0 | Hierarchy level (0=root) |
| `min_cells_for_subtype` | int | 1000 | Minimum cells for recursion |

#### `prepare_subset_dataset()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nametag` | str | Required | Dataset identifier |
| `input_file` | str | Required | Path to h5ad file |
| `cell_id_file` | str | None | Cell IDs to subset (TSV) |
| `batch_key` | List[str] | None | Batch correction variables |
| `integration` | bool | False | Enable Harmony |
| `cpus` | int | 16 | Number of CPUs |

### Cell Type Hierarchy JSON Format

```json
{
  "Major Cell Type": {
    "definition": "Detailed biological definition",
    "children": {
      "Subtype 1": {
        "definition": "Subtype definition",
        "children": {}
      },
      "Subtype 2": {
        "definition": "Subtype definition",
        "children": {
          "Sub-subtype": {
            "definition": "...",
            "children": {}
          }
        }
      }
    }
  }
}
```

---

## Output Files

### Directory Structure

```
results/
├── annotation_summary.html         # Interactive HTML report
├── final_annotations_combined.tsv  # All annotations merged
├── pipeline_report.html            # Nextflow execution report
├── pipeline_timeline.html          # Timeline visualization
├── pipeline_trace.txt              # Resource usage
└── annotations/
    └── work/
        ├── responses/              # LLM response JSONs
        ├── figures/                # UMAP and cluster trees
        ├── data/                   # Processed h5ad files
        └── logs/                   # Execution logs
```

### Final Annotations Format

```tsv
cell_id              leiden_0_15  ann                  annotation_level  nametag
AAACCTGAGCGATATA_1  0            Epithelial Cell      lev0              lev0
AAACCTGAGCGATATA_1  2            AT2 Cell            lev1              lev1_Epithelial_Cell
```

### LLM Response JSON

**Cluster Selection**:
```json
{
  "nametag": "lev0",
  "level": 0,
  "timestamp": "2025-10-27T11:00:00",
  "expected_cell_types": ["Epithelial Cell", "Immune Cell"],
  "cluster_selection_response": {
    "resolution": 0.15,
    "level_color": "red",
    "justification": "..."
  }
}
```

**Cell Type Annotation**:
```json
{
  "nametag": "lev0",
  "level": 0,
  "timestamp": "2025-10-27T11:05:00",
  "annotation_response": [
    {
      "cluster_id": "0",
      "cell_type_hypotheses": "Alveolar Type 2 Cell",
      "justification": "Strong SFTPC, SFTPA1 expression...",
      "key_markers_cited": ["SFTPC", "SFTPA1", "SLC34A2"],
      "confidence": "High"
    }
  ]
}
```

---

## Advanced Usage

### Custom Cell Thresholds

```bash
# For small datasets
python bin/run_annotation.py annotate \
  --min-cells 500 \
  --input data/small_dataset.h5ad \
  --tree data/lung.json \
  --celltype "Epithelial Cell" \
  --context "lung organoid"
```

### Manual Subtype Annotation

```python
from bin.lib import annotate_cell_types

# Annotate only a specific branch
collector = annotate_cell_types(
    expected_cells={
        "AT1 Cell": {"has_children": False, "definition": "..."},
        "AT2 Cell": {"has_children": True, "definition": "..."}
    },
    tree_file="data/lung.json",
    general_context="alveolar epithelium",
    nametag="lev2_alveolar",
    input_file="work/data/lev1_Epithelial_Cell_subset.h5ad",
    batch_key=None,
    integration=False,
    cpus_per_task=8,
    current_level=2,
    min_cells_for_subtype=500
)
```

### Batch Processing Multiple Datasets

```bash
for dataset in data/*.h5ad; do
  nextflow run main.nf \
    --input_h5ad "$dataset" \
    --tree_json data/lung.json \
    --context "lung tissue" \
    --outdir "results/$(basename $dataset .h5ad)"
done
```

### Error Recovery

```python
from bin.lib import TreeValidationError

try:
    collector = annotate_cell_types(...)
except TreeValidationError as e:
    print(f"Invalid tree structure: {e}")
except Exception as e:
    print(f"Annotation failed: {e}")
    # Check logs in work/logs/
```

---

## Troubleshooting

### Common Issues

#### 1. Import Errors

**Problem**: `ModuleNotFoundError: No module named 'scanpy'`

**Solution**:
```bash
source venv/bin/activate  # Activate virtual environment
pip install -r requirements.txt
```

#### 2. LLM API Errors

**Problem**: `google.api_core.exceptions.PermissionDenied`

**Solution**:
```bash
export GOOGLE_API_KEY="your-actual-api-key"
```

#### 3. Small Dataset Warnings

**Problem**: `Only 800 cells found for cell type 'X' (minimum 1000 required)`

**Solution**:
```bash
--min-cells 500  # Lower threshold
```

#### 4. Memory Issues

**Problem**: Process killed due to memory

**Solution**:
```groovy
// In nextflow.config
process {
    withName: 'ANNOTATE_CELL_TYPE' {
        memory = '32 GB'  // Increase memory
    }
}
```

#### 5. Nextflow Pipeline Failures

**Problem**: Pipeline exits with errors

**Check**:
```bash
# View Nextflow logs
cat .nextflow.log

# Check individual task logs
ls -la work/*/*/.command.log

# Validate syntax
nextflow config main.nf
```

### Debug Mode

```bash
# Enable detailed logging
export PYTHONPATH="${PWD}/bin"
python -c "
from lib.logger import PipelineLogger
import logging
logger = PipelineLogger.get_logger('debug', level=logging.DEBUG)
"

# Run with Nextflow debug
nextflow run main.nf -with-trace -with-report -with-timeline
```

---

## API Reference

### Core Functions

#### `annotate_cell_types()`
Main hierarchical annotation function. See [Configuration](#configuration) for parameters.

#### `prepare_subset_dataset()`
Preprocesses data: normalization, HVG, PCA, UMAP, clustering.

**Returns**: `Dict[str, str]` with keys:
- `encoded_cluster_tree`: Base64 PNG
- `subset_adata`: Path to processed h5ad
- `cluster_tree_file`: Path to cluster tree PNG

#### `prepare_celltype_inputs()`
Extracts features for annotation: markers, pathways, UMAP, adjacency.

**Returns**: `Dict` with keys:
- `celltyper_top_20_genes_based_on_scores`
- `celltyper_top_20_genes_based_on_pct_diff`
- `celltyper_top_20_pws`
- `consolidator_neighbors`
- `encoded_umap`

### Utility Functions

#### `find_node_data(tree_dict, cell_name, raise_on_missing=False)`
Recursively finds a node in the cell type tree.

**Raises**: `TreeValidationError` if `raise_on_missing=True` and node not found.

#### `get_immediate_children(node_data)`
Extracts direct children of a tree node.

**Returns**: `Dict[str, Dict[str, Union[str, bool]]]`

### Classes

#### `HierarchicalAnnotation`
Manages multi-level annotations.

**Methods**:
- `add_annotation(nametag, annotation_df, level)`: Add annotation
- `export_final_annotations(output_file)`: Export to TSV
- `get_level_annotations(level)`: Get all annotations at level

#### `CellTypeAnnotationClient`
Gemini API wrapper.

**Methods**:
- `select_cluster_resolution(prompt, image_data)`: Select resolution
- `annotate_cell_types(prompt, image_data)`: Annotate clusters

**Features**:
- Automatic retry (max 3 attempts)
- Robust JSON parsing
- Response validation

---

## Testing

```bash
# Module tests
python test_modules.py

# Pipeline tests
python test_pipeline.py

# Comprehensive improvement tests
python test_improvements.py

# All tests
./run_all_tests.sh
```

---

## Changelog

### v1.0.0 (Initial Release)

- Hierarchical LLM-based annotation
- Nextflow pipeline integration
- HTML report generation
- Parallel cell type processing
