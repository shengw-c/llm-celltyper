# Annotation Aggregation System

## Overview

The annotation aggregation system consolidates hierarchical cell type annotations from the entire pipeline into two comprehensive outputs:

1. **Hierarchical JSON** - A nested tree structure showing parent-child relationships
2. **Cell Mapping TSV** - A flat table mapping each cell to annotations at all hierarchy levels

## Quick Start

### Automatic Integration (Recommended)

The aggregation runs automatically as part of the Nextflow pipeline:

```bash
nextflow run main.nf \
    --input_h5ad data/test.h5ad \
    --tree_json data/lung.json
```

Find outputs in:
- `results/annotation/consolidated_annotations.json`
- `results/annotation/consolidated_annotations.tsv`

### Manual Testing

Test the aggregation independently:

```bash
./test_aggregation.sh
```

### Manual Execution

```bash
python bin/lib/aggregate_annotations.py \
    --responses-dir results/annotation/responses \
    --data-dir results/annotation/data \
    --output-json consolidated_annotations.json \
    --output-tsv consolidated_annotations.tsv \
    --pretty
```

## Output Formats

### 1. Hierarchical JSON

**File**: `consolidated_annotations.json`

**Structure**:
```json
{
  "nametag": "lev0_all_types",
  "level": 0,
  "timestamp": "2025-10-30T...",
  "annotation_response": [
    {
      "cluster_id": "0",
      "unique_id": "0",
      "cell_type_hypotheses": "Immune Cell",
      "justification": "...",
      "key_markers_cited": ["PTPRC", "CD45"],
      "confidence": "High",
      "hierarchy_path": "lev0_all_types",
      "children": [
        {
          "cluster_id": "0",
          "unique_id": "0.0",
          "cell_type_hypotheses": "Myeloid Cell",
          "hierarchy_path": "lev0_Immune_Cell",
          "children": [...]
        }
      ]
    }
  ]
}
```

**Key Features**:
- Nested structure preserves parent-child relationships
- `unique_id` provides globally unique cluster identifier (e.g., "0.0.1.2")
- `hierarchy_path` shows the nametag of this level
- `children` array contains sub-annotations
- All original annotation metadata preserved

**Use Cases**:
- Visualizing the cell type hierarchy tree
- Interactive exploration of annotation decisions
- Understanding the reasoning at each level
- Programmatic navigation of the hierarchy

### 2. Cell Mapping TSV

**File**: `consolidated_annotations.tsv`

**Structure**:
```
cell_id             ann_level0      ann_level1      ann_level2      ann_level3      ann_finest              unique_id
AAACCTGAGACAAAGG    Immune Cell     Myeloid Cell    Macrophage      Alveolar M...   Alveolar Macrophage     0.0.0.3
AAACCTGAGAGATGGG    Epithelial Cell Airway Ep...    Ciliated Cell   Mature Cil...   Mature Ciliated Cell    4.0.5.8
```

**Columns**:
- `cell_id` - Cell barcode/identifier
- `ann_level0` through `ann_levelN` - Annotation at each hierarchy level
- `ann_finest` - The most specific (deepest) annotation for this cell
- `unique_id` - Unique cluster identifier (matches JSON hierarchy)

**Key Features**:
- One row per cell
- Columns for each annotation level (0 to N)
- Easy filtering and grouping by cell type
- Compatible with standard data analysis tools

**Use Cases**:
- Adding annotations to AnnData objects
- Statistical analysis of cell type distributions
- Filtering cells by annotation
- Downstream visualization (UMAP with colors)

## Unique ID System

Each cluster is assigned a unique hierarchical identifier:

- Root level: `"0"`, `"1"`, `"2"`, ...
- First subdivision: `"0.0"`, `"0.1"`, `"0.2"`, ...
- Second subdivision: `"0.0.0"`, `"0.0.1"`, `"0.0.2"`, ...
- And so on...

**Example**:
```
0                    → Immune Cell
├─ 0.0              → Myeloid Cell
│  ├─ 0.0.0         → Macrophage
│  │  └─ 0.0.0.3    → Alveolar Macrophage
│  └─ 0.0.1         → Monocyte
└─ 0.1              → Lymphoid Cell
```

**Benefits**:
- Globally unique across entire hierarchy
- Encodes parent-child relationships
- Sortable and comparable
- Easy to determine relationship between clusters

## Usage Examples

### Example 1: Load and Explore Hierarchy

```python
import json

# Load hierarchical JSON
with open('consolidated_annotations.json', 'r') as f:
    data = json.load(f)

# Print root-level cell types
for annotation in data['annotation_response']:
    cell_type = annotation['cell_type_hypotheses']
    unique_id = annotation['unique_id']
    n_children = len(annotation.get('children', []))
    print(f"{cell_type} (ID: {unique_id}) - {n_children} subtypes")
```

### Example 2: Analyze Cell Distributions

```python
import pandas as pd

# Load cell mapping
df = pd.read_csv('consolidated_annotations.tsv', sep='\t')

# Count cells by major type
major_counts = df['ann_level0'].value_counts()
print("Major cell types:")
print(major_counts)

# Get finest annotation distribution
finest_counts = df['ann_finest'].value_counts().head(20)
print("\nTop 20 finest annotations:")
print(finest_counts)
```

### Example 3: Filter Cells by Type

```python
import pandas as pd

df = pd.read_csv('consolidated_annotations.tsv', sep='\t')

# Get all immune cells
immune_cells = df[df['ann_level0'] == 'Immune Cell']
print(f"Found {len(immune_cells)} immune cells")

# Get all macrophages (at any level)
level_cols = [c for c in df.columns if c.startswith('ann_level')]
macrophages = df[df[level_cols].apply(lambda row: 'Macrophage' in row.values, axis=1)]
print(f"Found {len(macrophages)} macrophages")

# Export specific cell types
macrophages.to_csv('macrophages_only.tsv', sep='\t', index=False)
```

### Example 4: Find Hierarchy Path

```python
import json

def find_path(data, unique_id, path=[]):
    """Find the full hierarchy path for a cluster."""
    for annotation in data.get('annotation_response', []):
        current_id = annotation.get('unique_id', '')
        cell_type = annotation['cell_type_hypotheses']
        
        if current_id == unique_id:
            return path + [cell_type]
        
        if 'children' in annotation:
            child_data = {'annotation_response': annotation['children']}
            result = find_path(child_data, unique_id, path + [cell_type])
            if result:
                return result
    return None

# Load data
with open('consolidated_annotations.json', 'r') as f:
    data = json.load(f)

# Find path for a specific cluster
path = find_path(data, '0.0.0.3')
print(" -> ".join(path))  # Output: Immune Cell -> Myeloid Cell -> Macrophage -> Alveolar Macrophage
```

### Example 5: Add to AnnData

```python
import pandas as pd
import scanpy as sc

# Load annotations
annotations_df = pd.read_csv('consolidated_annotations.tsv', sep='\t', index_col='cell_id')

# Load AnnData object
adata = sc.read_h5ad('data.h5ad')

# Add annotations to obs
for col in annotations_df.columns:
    adata.obs[col] = annotations_df.loc[adata.obs.index, col]

# Now you can plot with annotations
sc.pl.umap(adata, color='ann_finest')
```

## Advanced Usage

### Complete Example Script

See `example_use_consolidated_outputs.py` for a comprehensive example showing:
- Loading both output formats
- Navigating the hierarchy
- Computing statistics
- Finding cells by cell type
- Exporting subsets

Run it with:
```bash
python example_use_consolidated_outputs.py
```

### Integration with Nextflow

The AGGREGATION process is automatically included in the workflow:

```groovy
AGGREGATION(
    all_annotation_responses,  // All *_cell_annotation.json files
    all_tsv_files              // All *_subset.tsv files
)
```

Outputs are published to `${params.outdir}/annotation/` and used by downstream processes.

## Troubleshooting

### Missing Annotations

**Problem**: Some cells have `NaN` or `Unknown` in annotation columns.

**Explanation**: 
- Cells not processed at deeper levels will have fewer annotation levels filled
- Some clusters may be marked as "Unknown" by the LLM
- This is expected and indicates the annotation hierarchy depth varies by cell type

### Column Order

**Problem**: Annotation levels seem out of order (e.g., ann_level9 before ann_level10).

**Status**: This is actually correct! The script uses numeric sorting. If you see ann_level9, it means:
- There's no ann_level8 in your data (hierarchy gap)
- Or ann_level8 exists but sorting displayed incorrectly in console

Verify with: `head -1 file.tsv | tr '\t' '\n'`

### Performance

**Problem**: Aggregation is slow for large datasets.

**Solutions**:
- Aggregation scales with number of annotation files, not cells
- Typical runtime: < 1 minute for 30 annotation files and 500K cells
- Use `--pretty` flag only for human-readable output (increases file size)

## File Formats

### JSON Schema

The consolidated JSON follows this schema:

```
Root Object:
  - nametag: string (e.g., "lev0_all_types")
  - level: integer (0)
  - timestamp: ISO datetime string
  - general_context: string
  - expected_cell_types: array of strings
  - resolution: float
  - annotation_response: array of Annotation Objects

Annotation Object:
  - cluster_id: string
  - unique_id: string (hierarchical)
  - cell_type_hypotheses: string
  - justification: string
  - key_markers_cited: array of strings
  - confidence: string (High/Medium/Low)
  - hierarchy_path: string (nametag of this level)
  - children: array of Annotation Objects (recursive)
```

### TSV Schema

```
Column         Type     Description
-----------    ------   --------------------------------------------------
cell_id        string   Cell barcode/identifier
ann_level0     string   Annotation at level 0 (major cell type)
ann_level1     string   Annotation at level 1 (sub-type)
...            ...      Additional levels as needed
ann_levelN     string   Annotation at deepest level
ann_finest     string   Most specific annotation (may skip levels)
unique_id      string   Hierarchical cluster ID (e.g., "0.0.1.2")
```

## References

- **Main Script**: `bin/lib/aggregate_annotations.py`
- **Pipeline Integration**: `main.nf` (AGGREGATION process)
- **Test Script**: `test_aggregation.sh`
- **Example Usage**: `example_use_consolidated_outputs.py`
- **Implementation Notes**: `AGGREGATION_IMPLEMENTATION.md`

## Support

For issues or questions:
1. Check this README
2. Review `AGGREGATION_IMPLEMENTATION.md` for technical details
3. Run `test_aggregation.sh` to verify your setup
4. Check Nextflow logs in `work/` directory
