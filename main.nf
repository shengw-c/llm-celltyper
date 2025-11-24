#!/usr/bin/env nextflow

/*
 * Cell Type Annotation Pipeline v2
 * 
 * This pipeline performs hierarchical cell type annotation on single-cell RNA-seq data
 * using iterative LLM-based annotation with automatic cluster refinement.
 *
 * Inputs:
 *   - QCed h5ad file
 *
 * Outputs:
 *   - Hierarchical annotation results
 *   - Cluster data and figures
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
    Cell Type Annotation Pipeline v2
    ============================================================================
    
    Usage:
      nextflow run main.nf --input_h5ad <file> [options]
    
    Required Arguments:
      --input_h5ad FILE       QCed h5ad file with single-cell data
    
    Optional Arguments:
      --context STRING        Biological context description (default: "single-cell RNA-seq data")
      --batch_key STRING      Batch key for integration (default: null)
      --integration           Enable Harmony integration (default: false)
      --cpus INT              Number of CPUs per task (default: 16)
      --outdir DIR            Output directory (default: results)
      --min_cells INT         Min cells for subtype annotation (default: 1000)
      
      LLM Configuration:
      --llm_model STRING       LLM model for annotation and harmonization (default: gemini-2.5-pro)
      --llm_max_retries INT    Maximum retries for LLM API calls (default: 3)
      --mock_llm               Enable mock LLM mode for testing (no API calls)
      
      Clustering Configuration:
      --max_resolution FLOAT  Maximum resolution for Leiden clustering (default: 1.0)
      --transition_cutoff FLOAT  Stability threshold for clustering (default: 0.1)
      
      Feature Extraction:
      --top_genes INT         Number of top marker genes per cluster (default: 20)
      --top_pathways INT      Number of top pathways per cluster (default: 20)
      
      GSEA Configuration:
      --gsea_databases STRING Comma-separated list of gene set databases 
                              (default: MSigDB_Hallmark_2020,KEGG_2021_Human)
      
      Logging Configuration:
      --log_level STRING      Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)

      --help                  Show this help message
    
    Profiles:
      -profile standard       Local execution (default)
      -profile slurm          SLURM cluster execution
      -profile test           Minimal resources for testing with mock LLM
    
    Example:
      nextflow run main.nf \\
        --input_h5ad data/test.h5ad \\
        --context "lung tissue from healthy adult" \\
        --batch_key donor_id \\
        --cpus 8 \\
        --mock_llm \\
        -profile test
    
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

/*
 * Print pipeline configuration
 */
log.info """
============================================================================
Cell Type Annotation Pipeline v2
============================================================================
Input h5ad           : ${params.input_h5ad}
Context              : ${params.context}
Batch key            : ${params.batch_key ?: 'none'}
Integration          : ${params.integration}
CPUs per task        : ${params.cpus}
Output dir           : ${params.outdir}
Min cells            : ${params.min_cells}

LLM Configuration:
  Model              : ${params.llm_model}
  Max retries        : ${params.llm_max_retries}
  Mock mode          : ${params.mock_llm}

Clustering Configuration:
  Max resolution     : ${params.max_resolution}
  Transition cutoff  : ${params.transition_cutoff}

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
 * Process: Run hierarchical cell type annotation with harmonization
 * This process performs clustering, annotation, and harmonization in one step
 */
process ANNOTATE_CELL_TYPES {
    tag "Hierarchical annotation"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "responses/*.json"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "figures/*.png"
    publishDir "${params.outdir}/annotation/", mode: 'copy', pattern: "data/*.{tsv,csv}"
    publishDir "${params.outdir}/logs/", mode: 'copy', pattern: "logs/*.log"
    
    input:
    path input_h5ad
    
    output:
    path "responses/*.json", emit: annotation_responses, optional: true
    path "figures/*.png", emit: figures, optional: true
    path "data/*.tsv", emit: tsv_files, optional: true
    path "data/*.csv", emit: csv_files, optional: true
    path "data/hierarchical_annotation_complete.csv", emit: final_annotations, optional: true
    path "logs/*.log", emit: logs, optional: true
    
    script:
    def integration_flag = params.integration ? '--integration' : ''
    def batch_key_arg = params.batch_key ? "--batch-key ${params.batch_key}" : ''
    def mock_llm_env = params.mock_llm ? 'export MOCK_LLM_MODE=1' : ''
    """
    #!/bin/bash
    set -e
    
    ${mock_llm_env}
    
    # Run hierarchical annotation with harmonization
    python ${projectDir}/bin/run_annotation.py \\
        --input ${input_h5ad} \\
        --context "${params.context}" \\
        ${batch_key_arg} \\
        ${integration_flag} \\
        --cpus ${params.cpus} \\
        --output-dir . \\
        --min-cells ${params.min_cells} \\
        --max-resolution ${params.max_resolution} \\
        --transition-cutoff ${params.transition_cutoff} \\
        --llm-model "${params.llm_model}" \\
        --llm-max-retries ${params.llm_max_retries} \\
        --gsea-databases "${params.gsea_databases}" \\
        --top-genes ${params.top_genes} \\
        --top-pathways ${params.top_pathways} \\
        --log-level ${params.log_level}
    
    echo "✓ Annotation and harmonization completed"
    """
}

/*
 * Main workflow - Single-process pipeline with integrated harmonization
 */
workflow {
    // Input channel
    input_h5ad_ch = Channel.fromPath(params.input_h5ad, checkIfExists: true)
    
    // Run annotation with harmonization
    ANNOTATE_CELL_TYPES(input_h5ad_ch)
}

/*
 * Completion notification
 */
workflow.onComplete {
    def status_msg = workflow.success ? 'SUCCESS' : 'FAILED'
    def status_icon = workflow.success ? '✓' : '✗'
    
    log.info """
    ============================================================================
    Pipeline completed!
    ============================================================================
    Status              : ${status_icon} ${status_msg}
    Duration            : ${workflow.duration}
    Results dir         : ${params.outdir}
    Final annotations   : ${params.outdir}/annotation/data/hierarchical_annotation_complete.csv
    Responses           : ${params.outdir}/annotation/responses/
    Figures             : ${params.outdir}/annotation/figures/
    ============================================================================
    """.stripIndent()
    
    if (!workflow.success) {
        log.error """
        ⚠️  Pipeline failed. Check the logs for details.
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
