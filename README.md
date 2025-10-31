# LLM CellTyper

**Hierarchical, LLM-based cell type annotation for single-cell RNA-seq data**

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Nextflow](https://img.shields.io/badge/nextflow-21.04.0+-brightgreen.svg)](https://www.nextflow.io/)

## 🚀 Overview

LLM CellTyper is an automated pipeline that leverages Large Language Models (Gemini) to perform hierarchical cell type annotation on single-cell RNA-seq data. It combines computational biology techniques with AI-powered decision making to provide accurate, reproducible cell type annotations.

### ✨ Key Features

- **🤖 Automated Annotation**: LLM-based cell type assignment with confidence scores
- **🌳 Hierarchical Processing**: Recursive annotation from major cell types to subtypes  
- **⚡ Parallel Execution**: Process multiple cell types simultaneously via Nextflow
- **📊 Interactive Reports**: Beautiful HTML summaries with visualizations
- **📝 Complete Audit Trail**: All LLM decisions saved to JSON
- **⚙️ Resource Control**: Configurable CPU and memory usage
- **🔍 Advanced Logging**: Timestamped logs with configurable verbosity levels

## Quick Start

### Using Nextflow Pipeline (Recommended)

```bash
# Run annotation for all major cell types in parallel
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --context "lung tissue from healthy adult" \
  --cpus 16

# View results
firefox results/annotation_summary.html
```

### Using Python Directly

```python
from cell_annotator import celltype_annotor

celltype_annotor(
    expected_cells=my_cells,
    tree_file="data/lung.json",
    general_context="lung tissue",
    nametag="lev0",
    input_file="data/test.h5ad",
    batch_key=["batch"],
    cpus_per_task=16
)
```

## Features

✅ **Automated Annotation**: LLM-based cell type assignment with confidence scores  
✅ **Hierarchical Processing**: Recursive annotation of cell type subtypes  
✅ **Parallel Execution**: Process multiple cell types simultaneously  
✅ **Interactive Reports**: Beautiful HTML summaries with visualizations  
✅ **Complete Audit Trail**: All LLM decisions saved to JSON  
✅ **Resource Control**: Configurable CPU and memory usage  
✅ **Advanced Logging**: Timestamped logs with configurable verbosity levels  

## Documentation

- **[Complete Documentation](docs/README.md)**: Comprehensive guide covering all features, API reference, and troubleshooting
- **[LLM & SLURM Configuration](docs/LLM_AND_SLURM_CONFIGURATION.md)**: Guide for configuring LLM models and SLURM cluster execution
- **[Logging Usage Guide](docs/LOGGING_USAGE_GUIDE.md)**: Detailed guide for log levels, timestamps, and best practices
- **[Test Suite](tests/README.md)**: Test documentation and how to run tests
- **[Quick Start](#quick-start)**: Get started in 5 minutes

## 🔧 Installation

### Prerequisites

- Python >= 3.8
- Nextflow >= 21.04.0
- Google Gemini API access

### Setup

```bash
# Clone repository
git clone https://github.com/yourusername/llm_celltyper.git
cd llm_celltyper

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Configure Gemini API
export GOOGLE_API_KEY="your-api-key"
```

### Verify Installation

```bash
# Run module tests
python tests/test_modules.py

# Run pipeline tests  
python tests/test_pipeline.py

# Validate Nextflow configuration
nextflow config main.nf
```

## 📁 Project Structure

```
llm_celltyper/
├── main.nf                         # Nextflow pipeline
├── nextflow.config                 # Pipeline configuration
├── requirements.txt                # Python dependencies
├── README.md                       # This file
├── bin/                            # Executable scripts
│   ├── run_annotation.py           # CLI wrapper
│   ├── generate_summary.py         # HTML report generation
│   ├── analyze_logs.py             # Log analysis tool
│   └── lib/                        # Core library modules
│       ├── annotator.py            # Main annotation logic
│       ├── llm_client.py           # LLM API wrapper
│       ├── data_processing.py      # Data preprocessing
│       ├── tree_utils.py           # Tree navigation
│       ├── marker_genes.py         # Marker gene analysis
│       ├── pathway_enrichment.py   # Pathway analysis
│       ├── cluster_analysis.py     # Cluster relationships
│       ├── logger.py               # Logging system
│       └── prompts.py              # LLM prompts
├── data/                           # Example data
│   ├── lung.json                   # Cell type hierarchy example
│   └── test.h5ad                   # Example dataset
├── docs/                           # Documentation
│   ├── README.md                   # Complete documentation
│   ├── LLM_AND_SLURM_CONFIGURATION.md
│   └── LOGGING_USAGE_GUIDE.md
└── tests/                          # Test suite
    ├── test_modules.py             # Module tests
    └── test_pipeline.py            # Pipeline tests
```

## Usage Examples

### 1. List Available Cell Types

```bash
python run_annotation.py list --tree data/lung.json
```

### 2. Annotate Specific Cell Type

```bash
python run_annotation.py annotate \
  --input data/test.h5ad \
  --tree data/lung.json \
  --celltype "Endothelial Cell" \
  --context "lung tissue" \
  --cpus 16
```

### 3. Run Complete Pipeline

```bash
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --context "lung tissue from healthy adult" \
  --batch_key dataset \
  --integration \
  --cpus 16 \
  --outdir results
```

### 4. Generate Summary Report

```bash
python generate_summary.py
```

## Logging Configuration

The pipeline now includes enhanced logging features with timestamps and configurable verbosity levels.

### Available Log Levels

- **DEBUG**: Most verbose - shows all messages including detailed execution flow
- **INFO**: Default - shows informational messages, warnings, and errors
- **WARNING**: Shows only warnings and errors
- **ERROR**: Least verbose - shows only error messages

### Using Log Levels

#### With Nextflow Pipeline

```bash
# Enable debug logging for troubleshooting
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --log_level DEBUG

# Reduce verbosity (production mode)
nextflow run main.nf \
  --input_h5ad data/test.h5ad \
  --tree_json data/lung.json \
  --log_level WARNING
```

#### With Python Scripts Directly

```bash
# Enable debug logging
python bin/run_annotation.py annotate \
  --input data/test.h5ad \
  --tree data/lung.json \
  --celltype "T cell" \
  --context "lung tissue" \
  --log-level DEBUG

# Analyze logs with minimal output
python bin/analyze_logs.py \
  --log-files logs/*.log \
  --log-level ERROR
```

### Log Output Format

All console and file logs now include timestamps:

```
2025-10-29 14:30:15 | INFO | Starting annotation for: T cell
2025-10-29 14:30:16 | INFO | Configuration:
2025-10-29 14:30:16 | INFO |   Input: test.h5ad
2025-10-29 14:30:17 | DEBUG | Clustering with resolution 0.5
```

### Important Notes

- **Console logs** respect the `--log-level` parameter
- **Log files** (in `logs/` directory) always contain full DEBUG-level information
- Default log level is **INFO** (same as before, backward compatible)
- For detailed logging guide, see `docs/LOGGING_USAGE_GUIDE.md`

## Testing

All test files are located in the `tests/` directory.

```bash
# Run individual test suites
python tests/test_modules.py      # Module & utility tests (4 tests)
python tests/test_pipeline.py     # Pipeline & config tests (5 tests)

# Run all tests at once
for test in tests/test_*.py; do python $test || exit 1; done

# Quick validation
nextflow config main.nf           # Validate Nextflow syntax

# Test logging improvements
python test_logging_simple.py     # Verify timestamp and log level features
```

**Expected Results**: All tests should pass (18/18 tests passing)

See `tests/README.md` for detailed test documentation.

## 📋 Output

The pipeline generates:

1. **📊 Annotation Results** (TSV files): Cell type assignments per cell
2. **📄 Response JSONs**: Complete LLM decision audit trail  
3. **📈 Visualizations**: UMAP plots and cluster trees
4. **🌐 Interactive HTML**: Comprehensive annotation summary

## 🧪 Testing

```bash
# Run individual test suites
python tests/test_modules.py      # Module & utility tests
python tests/test_pipeline.py     # Pipeline & config tests

# Run all tests
for test in tests/test_*.py; do python $test || exit 1; done

# Validate Nextflow syntax
nextflow config main.nf
```

## 📚 Documentation

- **[Complete Documentation](docs/README.md)**: Comprehensive guide with API reference
- **[LLM & SLURM Configuration](docs/LLM_AND_SLURM_CONFIGURATION.md)**: LLM models and cluster setup
- **[Logging Guide](docs/LOGGING_USAGE_GUIDE.md)**: Log levels and best practices
- **[Test Documentation](tests/README.md)**: How to run and write tests

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request
