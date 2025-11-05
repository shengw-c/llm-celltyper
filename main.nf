#!/usr/bin/env nextflow

/*
 * Cell Type Annotation Pipeline
 * 
 * This pipeline performs hierarchical cell type annotation on single-cell RNA-seq data
 * using LLM-based annotation for each major cell type in parallel.
 *
 * Inputs:
 *   - QCed h5ad file
 *   - Cell type hierarchy JSON file
 *
 * Outputs:
 *   - Annotation results for each cell type
 *   - Interactive HTML summary report
 */

nextflow.enable.dsl=2

/*
 * Pipeline parameters are defined in nextflow.config
 * They can be overridden via command line: --param_name value
 */

/*
 * Print help message
 */
def helpMessage() {
    log.info """
    ============================================================================
    Cell Type Annotation Pipeline
    ============================================================================
    
    Usage:
      nextflow run main.nf --input_h5ad <file> --tree_json <file> [options]
    
    Required Arguments:
      --input_h5ad FILE       QCed h5ad file with single-cell data
      --tree_json FILE        Cell type hierarchy JSON file
    
    Optional Arguments:
      --context STRING        Biological context description (default: "single-cell RNA-seq data")
      --batch_key STRING      Batch key for integration (default: null)
      --integration           Enable Harmony integration (default: false)
      --cpus INT              Number of CPUs per task (default: 16)
      --outdir DIR            Output directory (default: results)
      --min_cells INT         Minimum cells for subtype annotation (default: 1000)
      --condition_num INT     Number of expected condition-driven groups or artifacts (default: 1)
      
      LLM Configuration:
      --llm_model_general STRING      LLM model for general purpose (default: gemini-2.5-flash)
      --llm_model_complicated STRING  LLM model for complicated tasks (default: gemini-2.5.pro)
      --llm_nreplies INT              Number of replies to request from LLM (default: 1)
      --llm_max_retries INT    Maximum retries for LLM API calls (default: 3)
      
      Feature Extraction:
      --top_genes INT         Number of top marker genes per cluster (default: 20)
      --top_pathways INT      Number of top pathways per cluster (default: 20)
      
      GSEA Configuration:
      --gsea_databases STRING Comma-separated list of gene set databases (default: MSigDB_Hallmark_2020,KEGG_2021_Human)
      
      Logging Configuration:
      --log_level STRING      Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)

      --help                  Show this      help message
    
    Profiles:
      -profile standard       Local execution (default)
      -profile slurm          SLURM cluster execution
      -profile slurm_gpu      SLURM with GPU support
      -profile test           Minimal resources for testing
      -profile production     Production with enhanced monitoring
    
    Example:
      nextflow run main.nf \\
        --input_h5ad data/test.h5ad \\
        --tree_json data/lung.json \\
        --context "lung tissue from healthy adult" \\
        --batch_key donor_id \\
        --cpus 8 \\
        --llm_model gemini-2.0-flash-exp \\
        -profile slurm
    
    ============================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
 * Validate required parameters
 */
if (!params.input_h5ad) {
    log.error "ERROR: --input_h5ad is required"
    helpMessage()
    exit 1
}

if (!params.tree_json) {
    log.error "ERROR: --tree_json is required"
    helpMessage()
    exit 1
}

/*
 * Print pipeline configuration
 */
log.info """
============================================================================
Cell Type Annotation Pipeline
============================================================================
Input h5ad           : ${params.input_h5ad}
Tree JSON            : ${params.tree_json}
Context              : ${params.context}
Batch key            : ${params.batch_key ?: 'none'}
Integration          : ${params.integration}
CPUs per task        : ${params.cpus}
Output dir           : ${params.outdir}
Min cells            : ${params.min_cells}
Condition num        : ${params.condition_num}

LLM Configuration:
  Model (general).   : ${params.llm_model_general}
  Model (complicated): ${params.llm_model_complicated}
  Max retries        : ${params.llm_max_retries}
  Number of replies  : ${params.llm_nreplies}

Feature Extraction:
  Top genes          : ${params.top_genes}
  Top pathways       : ${params.top_pathways}

GSEA Configuration:
  Databases          : ${params.gsea_databases}

Logging Configuration:
  Log level          : ${params.log_level}
============================================================================
"""

/*
 * Process 1: Split dataset by all top-level cell types
 * This runs ONCE on the full dataset and creates separate .h5ad files
 */
process SPLIT_DATASET {
    tag "Split by all cell types"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "responses/*.json"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "figures/*.png"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "data/*.tsv"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "data/*.h5ad"
    publishDir "${params.outdir}/logs/", mode: 'copy', pattern: "logs/*.log", saveAs: { filename -> "SPLIT_${filename.split('/')[-1]}" }
    
    input:
    tuple path(input_h5ad), path(tree_json)
    
    output:
    path "split/*.h5ad", emit: split_files
    path "responses/*_cluster_selection.json", emit: cluster_responses, optional: true
    path "responses/*_cell_annotation.json", emit: annotation_responses, optional: true
    path "figures/*.png", emit: figures, optional: true
    path "data/*.h5ad", emit: h5ad_files, optional: true
    path "data/*.tsv", emit: tsv_files, optional: true
    path "logs/*.log", emit: logs, optional: true
    
    script:
    def integration_flag = params.integration ? '--integration' : ''
    def batch_key_arg = params.batch_key ? "--batch-key ${params.batch_key}" : ''
    """
    #!/bin/bash
    set -e
    source ${projectDir}/venv/bin/activate
    
    # Run split mode - processes ALL top-level cell types at once
    python ${projectDir}/bin/run_annotation.py split \\
        --input ${input_h5ad} \\
        --tree ${tree_json} \\
        --context "${params.context}" \\
        ${batch_key_arg} \\
        ${integration_flag} \\
        --cpus ${params.cpus} \\
        --output-dir . \\
        --min-cells ${params.min_cells} \\
        --condition-num ${params.condition_num} \\
        --llm-model-general "${params.llm_model_general}" \\
        --llm-model-complicated "${params.llm_model_complicated}" \\
        --llm-nreplies ${params.llm_nreplies} \\
        --gsea-databases "${params.gsea_databases}" \\
        --top-genes ${params.top_genes} \\
        --top-pathways ${params.top_pathways} \\
        --log-level ${params.log_level}
    
    echo "✓ Split dataset completed"
    ls -lh split/
    """
}

/*
 * Process 2: Run nested annotation for each split file
 * This runs in PARALLEL for each cell type
 */
process RUN_NESTED_ANNOTATION {
    tag "${split_file.baseName}"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "responses/*.json"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "figures/*.png"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "data/*.tsv" 
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "data/*.h5ad"
    publishDir "${params.outdir}/logs/", mode: 'copy', pattern: "logs/*.log", saveAs: { filename -> 
        def celltype = split_file.baseName
        def logname = filename.split('/')[-1]
        // If log already starts with cell type name, keep as-is, otherwise prefix it
        logname.startsWith(celltype) ? logname : "${celltype}_${logname}"
    }
    errorStrategy 'ignore'  // Continue with other cell types even if one fails
    
    input:
    tuple path(split_file), path(tree_json)
    
    output:
    path "responses/*_cluster_selection.json", emit: cluster_responses, optional: true
    path "responses/*_cell_annotation.json", emit: annotation_responses, optional: true
    path "figures/*.png", emit: figures, optional: true
    path "data/*.h5ad", emit: h5ad_files, optional: true
    path "data/*.tsv", emit: tsv_files, optional: true
    path "data/*_final_annotations.tsv", emit: final_annotations, optional: true
    path "logs/*.log", emit: logs, optional: true
    path ".exitcode", emit: exitcode, optional: true
    
    script:
    // Extract cell type name from filename (remove .h5ad extension)
    def cell_type = split_file.baseName.replace('_', ' ')
    def integration_flag = params.integration ? '--integration' : ''
    def batch_key_arg = params.batch_key ? "--batch-key ${params.batch_key}" : ''
    """
    #!/bin/bash
    set -e
    
    source ${projectDir}/venv/bin/activate
    
    # Run nested annotation on the split file
    python ${projectDir}/bin/run_annotation.py annotate \\
        --input ${split_file} \\
        --tree ${tree_json} \\
        --celltype "${cell_type}" \\
        --context "${params.context}" \\
        ${batch_key_arg} \\
        ${integration_flag} \\
        --cpus ${params.cpus} \\
        --output-dir . \\
        --min-cells ${params.min_cells} \\
        --condition-num ${params.condition_num} \\
        --llm-model-general "${params.llm_model_general}" \\
        --llm-model-complicated "${params.llm_model_complicated}" \\
        --llm-nreplies ${params.llm_nreplies} \\
        --gsea-databases "${params.gsea_databases}" \\
        --top-genes ${params.top_genes} \\
        --top-pathways ${params.top_pathways} \\
        --log-level ${params.log_level}
    
    # Save exit code for error tracking
    echo \$? > .exitcode
    
    # Log completion status
    if [ \$? -eq 0 ]; then
        echo "✓ Successfully annotated ${cell_type}" >> logs/status.log
    else
        echo "✗ Failed to annotate ${cell_type}" >> logs/status.log
    fi
    """
}

/*
* Process 3: Aggregate responses and annotations
*/
process AGGREGATION {
    tag "Aggregate annotations"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "consolidated_*.{json,tsv}"
    
    input:
    path annotation_responses, stageAs: "responses/*"
    path tsv_files, stageAs: "data/*"
    path tree_json
    
    output:
    path "consolidated_annotations.json", emit: consolidated_json
    path "consolidated_annotations.tsv", emit: consolidated_tsv
    
    script:
    """
    #!/bin/bash
    set -e

    source ${projectDir}/venv/bin/activate
    
    # Run the aggregation script
    python ${projectDir}/bin/lib/aggregate_annotations.py \\
        --responses-dir responses \\
        --data-dir data \\
        --tree-json ${tree_json} \\
        --output-json consolidated_annotations.json \\
        --output-tsv consolidated_annotations.tsv \\
        --pretty
    
    echo "✓ Aggregation completed"
    """
}


/*
 * Process 3: Analyze logs using LLM
 */
process ANALYZE_LOGS {
    tag "Analyze logs"
    
    input:
    path files, stageAs: "logs_input/?/*"  // Put all logs in a subdirectory to avoid name collisions
    
    output:
    path "log_analysis.json", emit: log_analysis, optional: true
    
    script:
    """
    #!/bin/bash
    set -e

    source ${projectDir}/venv/bin/activate
    
    python3 << 'EOF'
import os
import shutil
import glob

# Run log analysis
import sys
sys.path.insert(0, "${projectDir}/bin")

# Only run if we have log files
log_files = [f for f in glob.glob("logs_input/**/*") if f.endswith('.log')]
log_files = ",".join(log_files)
if log_files:
    print(f"Found {len(log_files)} log files for analysis")
    import subprocess
    subprocess.run([
       "python", "${projectDir}/bin/analyze_logs.py",
        "--log-files", log_files,
        "--output", "log_analysis.json",
        "--log-level", "${params.log_level}",
        "--llm-model", "${params.llm_model_general}"
    ])
else:
    print("No log files found, skipping log analysis!")
EOF
    """
}

/*
 * Process 4: Generate HTML summary report
 */
process GENERATE_SUMMARY {
    tag "Generate summary"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path cluster_responses, stageAs: "responses/*"
    path annotation_responses, stageAs: "responses/*"
    path figures, stageAs: "figures/*"
    path log_analysis
    path consolidated_annotations  // Single consolidated TSV file from AGGREGATION
    path lev0_annotations
    
    output:
    path "annotation_summary.html", emit: html_report
    path "final_annotations_combined.tsv", emit: final_annotations_combined, optional: true
    
    script:
    """
    #!/bin/bash
    set -e
    source ${projectDir}/venv/bin/activate
    python3 << 'EOF'
import os
import shutil
import pandas as pd
import numpy as np

# Process consolidated annotations from AGGREGATION
if os.path.exists("${consolidated_annotations}"):
    print("Reading consolidated annotations from AGGREGATION...")
    combined_df = pd.read_csv("${consolidated_annotations}", sep='\\t')
    
    # Merge with level 0 annotations if provided
    if os.path.exists("${lev0_annotations}"):
        print("Merging with level 0 annotations...")
        lev0_df = pd.read_csv("${lev0_annotations}", sep='\\t', index_col=0)
        lev0_df = (
            lev0_df.reset_index(names='cell_id')
            .rename(columns={"ann":"ann_major"})
            [['cell_id', 'ann_major']]
        )
        
        # Merge on cell_id
        combined_df = lev0_df.merge(combined_df, on='cell_id', how='left')
    
    # Make sure ann_major is the first annotation column, and ann_finest is the last
    ann_cols = [c for c in combined_df.columns if c.startswith('ann_')]
    
    # Ensure ann_finest exists and is the deepest annotation (not NaN)
    # Replace "Unknown" as NaN to avoid the finest being "Unknown"
    level_cols = sorted([c for c in combined_df.columns if c.startswith('ann_level')])
    if level_cols:
        combined_df['ann_finest'] = combined_df.replace("Unknown", np.nan)[level_cols].ffill(axis=1).iloc[:, -1]
    
    # Order columns: cell_id, ann_major, ann_level0, ann_level1, ..., ann_finest, unique_id
    base_cols = ['cell_id']
    if 'ann_major' in combined_df.columns:
        base_cols.append('ann_major')
    
    remaining_ann_cols = sorted([c for c in ann_cols if c not in ['ann_major', 'ann_finest'] and c.startswith('ann_level')])
    
    end_cols = []
    if 'ann_finest' in combined_df.columns:
        end_cols.append('ann_finest')
    if 'unique_id' in combined_df.columns:
        end_cols.append('unique_id')
    
    ordered_cols = base_cols + remaining_ann_cols + end_cols
    combined_df = combined_df[ordered_cols]
    
    # Save to file
    combined_df.to_csv("final_annotations_combined.tsv", sep='\\t', index=False)
    print(f"✓ Created final_annotations_combined.tsv")
    print(f"  Total cells annotated: {len(combined_df)}")
    
    # Show column summary
    ann_cols = [c for c in combined_df.columns if c.startswith('ann_')]
    print(f"  Annotation columns: {', '.join(ann_cols)}")
else:
    print("ERROR: No consolidated annotations found from AGGREGATION!")

# Generate HTML summary
import sys
sys.path.insert(0, "${projectDir}/bin")
from generate_summary import create_annotation_summary

html_file = create_annotation_summary(
    work_dir=".",
    output_file="annotation_summary.html"
)

print(f"✓ Summary report generated: {html_file}")
EOF
    """
}

/*
 * Main workflow - Two-Stage Optimized Pipeline
 */
workflow {
    // Input channels
    input_h5ad_ch = Channel.fromPath(params.input_h5ad, checkIfExists: true)
    tree_json_ch = Channel.fromPath(params.tree_json, checkIfExists: true)
    
    // STAGE 1: Split dataset (runs ONCE on full dataset)
    SPLIT_DATASET(
        input_h5ad_ch.combine(tree_json_ch)
    )
    
    // Collect stage 1 outputs
    split_cluster_responses = SPLIT_DATASET.out.cluster_responses
    split_annotation_responses = SPLIT_DATASET.out.annotation_responses
    split_figures = SPLIT_DATASET.out.figures
    split_logs = SPLIT_DATASET.out.logs
    split_annotation = SPLIT_DATASET.out.tsv_files
    
    // STAGE 2: Nested annotation for each split file (runs in PARALLEL)
    // Flatten the split files channel and combine each with tree_json
    split_files_ch = SPLIT_DATASET.out.split_files
        .flatten()
        .combine(tree_json_ch)
    
    RUN_NESTED_ANNOTATION(split_files_ch)
    
    // Collect all outputs from stage 2
    nested_cluster_responses = RUN_NESTED_ANNOTATION.out.cluster_responses.collect()
    nested_annotation_responses = RUN_NESTED_ANNOTATION.out.annotation_responses.collect()
    nested_figures = RUN_NESTED_ANNOTATION.out.figures.collect()
    nested_logs = RUN_NESTED_ANNOTATION.out.logs.collect()
    nested_final_annotations = RUN_NESTED_ANNOTATION.out.final_annotations.collect()
    
    // Combine outputs from both stages
    all_cluster_responses = split_cluster_responses.mix(nested_cluster_responses).collect()
    all_annotation_responses = split_annotation_responses.mix(nested_annotation_responses).collect()
    all_figures = split_figures.mix(nested_figures).collect()
    all_logs = split_logs.mix(nested_logs).collect()
    
    // Collect all TSV files for aggregation
    split_tsv_files = SPLIT_DATASET.out.tsv_files
    nested_tsv_files = RUN_NESTED_ANNOTATION.out.tsv_files.collect()
    all_tsv_files = split_tsv_files.mix(nested_tsv_files).collect()
    
    // Run aggregation to consolidate all annotations
    // Generates:
    // (1) hierarchical JSON with nested children and unique IDs for each cluster
    // (2) cell-to-cluster mapping TSV with annotations at each level
    AGGREGATION(
        all_annotation_responses,
        all_tsv_files,
        tree_json_ch
    )

    // Skip log analysis due to file name collisions - analyze manually if needed
    ANALYZE_LOGS(all_logs)
    
    // Generate summary report
    GENERATE_SUMMARY(
        all_cluster_responses,
        all_annotation_responses,
        all_figures,
        ANALYZE_LOGS.out.log_analysis,
        AGGREGATION.out.consolidated_tsv,  // Use consolidated TSV from aggregation
        split_annotation // pass split annotations for level 0 merging
    )
}

/*
 * Completion notification
 */
workflow.onComplete {
    // Collect status information
    def status_msg = workflow.success ? 'SUCCESS' : 'FAILED'
    def status_icon = workflow.success ? '✓' : '✗'
    
    log.info """
    ============================================================================
    Pipeline completed!
    ============================================================================
    Status              : ${status_icon} ${status_msg}
    Duration            : ${workflow.duration}
    Results dir         : ${params.outdir}
    HTML report         : ${params.outdir}/annotation_summary.html
    Final annotations   : ${params.outdir}/final_annotations_combined.tsv
    
    Consolidated Outputs:
    Hierarchical JSON   : ${params.outdir}/annotation/consolidated_annotations.json
    Cell mapping TSV    : ${params.outdir}/annotation/consolidated_annotations.tsv
    ============================================================================
    """.stripIndent()
    
    if (!workflow.success) {
        log.error """
        ⚠️  Some processes may have failed. Check the logs for details.
        Individual cell type failures are logged but don't stop the pipeline.
        Review ${params.outdir}/logs/ directory for error details.
        """
    }
}

workflow.onError {
    log.error """
    ============================================================================
    Pipeline execution error!
    ============================================================================
    Error message: ${workflow.errorMessage}
    Error report : ${workflow.errorReport}
    ============================================================================
    """
}
