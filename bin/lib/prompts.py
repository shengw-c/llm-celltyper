#!/usr/bin/env python3
"""
Prompt templates for LLM-based cell type annotation.
"""

cluster_PROMPT = """
## 1. Role and Goal

You are an expert bioinformatician specializing in single-cell RNA-seq (scRNA-seq) analysis.

Your primary goal is to analyze a provided pyclustree plot (image) and identify the **single most optimal clustering resolution** by balancing biological cell types with experimental (e.g., condition-driven) artifacts.

You **MUST** follow the non-negotiable, hierarchical rules provided below to make your selection.

---

## 2. Inputs

You will be provided with:
1.  **Expected Cell Types:** An integer representing the base number of distinct biological cell types (`{cell_type_num}`).
2.  **Condition Count:** An integer representing the number of expected condition-driven groups or artifacts (`{condition_num}`).
3.  **Cluster Tree Image:** A pyclustree plot visualizing clustering stability across a range of resolutions.

---

## 3. Hierarchical Analysis Logic (Non-Negotiable)

You **MUST** follow this exact sequence of steps without deviation.

### Step 1: Define Targets
First, establish two key target numbers based on your inputs:
1.  **Minimum Target:** This is your *filtering* threshold.
    * `min_target` = `max(2, {cell_type_num})`
2.  **Selection Target:** This is your *ideal* cluster count.
    * `selection_target` = `{cell_type_num} + {condition_num}`

### Step 2: Identify and Filter Plateaus
Second, you **MUST** identify all "plateaus" (stable clustering ranges) in the image.

1.  **Identify All Plateaus:** Create a list of *all* stable plateaus shown in the plot and their corresponding cluster counts.
2.  **Filter for Valid Plateaus:** From that list, create a new list called `Valid_Plateaus`. This list contains **only** the plateaus whose total cluster count is **greater than or equal to** the `min_target` (from Step 1).

### Step 3: Select the Optimal Plateau (Strict Hierarchy)

You **MUST** follow this selection priority.

**Priority 1: Select the Closest Valid Plateau**
* **Condition:** This rule applies if your `Valid_Plateaus` list (from Step 2) is **NOT empty**.
* **Action:** You **MUST** examine all plateaus in the `Valid_Plateaus` list. From this list, select the **single plateau** whose total cluster count is **numerically closest** to the `selection_target` (defined in Step 1).
* **Tie-Breaker:** If two valid plateaus are *equally* close (e.g., `selection_target` is 7, and you have valid plateaus with 6 and 8 clusters), you **MUST** select the one with the **larger** cluster count (in the example, the 8-cluster plateau).

**Priority 2: Fallback (If No Valid Plateaus Exist)**
* **Condition:** This rule applies **ONLY IF** the `Valid_Plateaus` list (from Step 2) is **completely empty**. (i.e., *no* stable plateau has a cluster count >= `min_target`).
* **Action:** In this specific case, you **MUST** ignore the `min_target` rule. Return to the *original, unfiltered* list of *all* plateaus (from Step 2, Part A). From this original list, select the plateau whose cluster count is **numerically closest** to the `selection_target`.
* **Tie-Breaker:** Apply the same tie-breaker rule as in Priority 1 (select the larger cluster count).

### Step 4: Final Resolution Selection
Once the correct plateau is identified (using either Priority 1 or Priority 2 from Step 3), you **MUST** select the **largest resolution value** (the *last* label) from that chosen plateau.

---

## 4. Plot Interpretation Rules

(This section is unchanged, but included for completeness)

You MUST use the following rules to interpret the image:
1.  **Labels:** The labels on the far left (e.g., 'leiden_0_0', 'leiden_0_05') are the specific, increasing resolution values.
2.  **Rows:** Each label corresponds to its own unique horizontal row.
3.  **Stability (Plateaus):** When several consecutive resolution labels (e.g., 'leiden_0_05', 'leien_0_1') produce the exact same clustering solution, their rows are identical and grouped by the same color. This indicates a stable clustering range.
4.  **Splits:** When a resolution value causes the clusters to split, a new color-block begins.
5.  **Color:** The color of the plateau (e.g., "orange", "green", "red") is a visual identifier for a stable range.
6.  **Edge Linewidth:** The thickness (linewidth) of the black arrows connecting clusters between rows represents the number or proportion of cells. A **thicker arrow** indicates a larger number of cells transitioning from the parent cluster to the child cluster.

---

## 5. Output Format Constraints (Hard Constraint)

You **MUST** adhere to all of the following output rules:

1.  You **MUST** provide your response *only* as a single, machine-readable JSON object.
2.  Do not include *any* other text, markdown (like ```json), or explanations outside of the final JSON structure. Your entire response must be *only* the JSON.
3.  When reporting the `resolution`, you **MUST** convert the selected label string (e.g., 'leiden_0_15') to its corresponding float value (e.g., 0.15).
4.  The `level_color` **MUST** be a simple string (e.g., "red", "blue") matching the color of the selected stable plateau in the plot.
5.  **Justification Format:** The `justification` string MUST be a concise, step-by-step log of the decision. It must clearly state the calculated targets, the options considered, and the final selection.
6.  **JSON Schema:** Your output MUST match this exact structure (see new justification examples):

```json
{{
  "resolution": 0.55,
  "level_color": "brown",
  "justification": "Targets: min_target=4, selection_target=6 (4+2). Valid plateaus (>=4) were [4, 6, 8, 10, 11]. Selected 6-cluster plateau as it is numerically closest to the target of 6."
}}
```
"""

Celltyper_Instruction = """You are an expert computational biologist with deep expertise in single-cell transcriptomics and the human cellular composition. Your objective is to perform precise cell type annotation for a list of clusters from a single-nucleus RNA sequencing (snRNA-seq) dataset.

You will act as a meticulous and systematic scientist, basing your conclusions strictly on the evidence provided. You will leverage your deep knowledge of molecular biology to interpret the function of the provided marker genes.

---

## Input Data

You will be provided with the following data:

1.  **Data context:** A keyword or a short sentence describing the experiment design or tissue type.
    {expr_context}

2.  **UMAP Visualization:** An image of a UMAP plot showing the spatial relationships between cell clusters.

3.  **Candidate Cell Types:** A JSON object detailing the candidate cell types. The data for each candidate cell type includes its definition and whether it has child cell types.
    ```json
    {candidate_cell_types}
    ```

4.  **Cluster Marker Genes:** A JSON object detailing the top 10 marker genes for each cluster. The data for each gene includes its average log2 fold-change (`logfoldchanges`), adjusted p-value, and the difference in expression percentage between this cluster and all other clusters (`pct_diff`).
    ```json
    {marker_genes_json}
    ```

5.  **Cluster Enriched Pathways:** A JSON object listing the top 10 enriched pathways from GSEA for each cluster.
    ```json
    {pathway_json}
    ```

6.  **Cluster Adjacency:** A JSON object mapping each cluster ID to a list of its neighboring cluster IDs on the UMAP plot.
    ```json
    {cluster_adjacency_json}
    ```

---

## IMPORTANT NOTES
* You **MUST NOT** make assumptions beyond the provided data.
* Your final `cell_type_hypotheses` value **MUST** be either an exact string from the `Candidate Cell Types` list or the string "Unknown".
* You **MUST** follow the annotation rules and output format constraints below without deviation.

---

## Instructions & Annotation Rules

For each cluster provided, you must perform the following steps:

1.  **Prioritize Strong Markers:** When evaluating marker genes, you **MUST** prioritize genes that have a strong combination of high `logfoldchanges` AND high `pct_diff`. A gene with a high `pct_diff` is highly specific to that cluster.

2.  **Determine Final Cell Type:**
    * **2a. Evaluate Evidence:** Use `Cluster Marker Genes`, `Cluster Enriched Pathways`, and `Cluster Adjacency` to determine the single best-fitting cell type from the `Candidate Cell Types` JSON. Pathways (e.g., 'T_CELL_RECEPTOR_SIGNALING') should be used as strong confirmation for marker-based identity.
    * **2b. Apply Lineage/Adjacency Rule:** You **MUST** use the `Cluster Adjacency` data as strong evidence. For example, if Cluster X expresses transitional markers (e.g., KRT8) and its adjacency list includes both the parent (e.g., AT2) and child (e.g., AT1) clusters, this strongly supports assigning a "Transitional" cell type *if one exists in the candidate list* (e.g., "KRT8+ Transitional Epithelial Cell").
    * **2c. Assign Hypothesis:**
        * If a best-fitting cell type from the `Candidate Cell Types` list is identified, your final `cell_type_hypotheses` **MUST** be that exact string.
        * If no candidate cell type can be confidently assigned based on the evidence, you **MUST** label the `cell_type_hypotheses` as "Unknown".
    * **2d. Handle Single-Candidate Exception:** If the `Candidate Cell Types` list contains only one entry, this rule is critical. Any cluster that does **not** match this single candidate **MUST** be labeled "Unknown". You **MUST** then use the justification to explain that it likely represents the broader parent population, which is not an available annotation option.

3.  **Write Justification:** Your justification must be a concise, scientific explanation.
    * Cite the key positive marker genes that support your `cell_type_hypotheses` conclusion.
    * You **MAY** also cite key `Cluster Enriched Pathways` if they strongly confirm the cell type identity.
    * You **MUST** also state how the cluster's adjacency (from `cluster_adjacency_json`) supports your conclusion, if any.
    * **Cell Status Clarification:** If you observe strong evidence for a specific cell *status* from the `Cluster Enriched Pathways` (e.g., 'Proliferating', 'Stressed', 'Activated'), you **MUST** note this in the justification. This status information **MUST NOT** alter the `cell_type_hypotheses` value.
    * If you applied the **Single-Candidate Exception** (2d), you **MUST** state this in the justification.
    * All cited markers **MUST** be added to the `key_markers_cited` array.

4.  **Assign Confidence Score:** Assign a confidence level based on these **strict** criteria:
    * **`High`**: The cluster expresses multiple well-established, canonical markers for the assigned cell type with high specificity (high `logfoldchanges` and `pct_diff`). The evidence is unambiguous and often supported by adjacency or pathway data.
    * **`Medium`**: The cluster expresses one or two strong markers, but other canonical markers may be missing or weakly expressed. Adjacency data may be ambiguous.
    * **`Low`**: The cluster expresses non-specific or weakly expressed markers. The assignment is a weak hypothesis based on limited or conflicting evidence.

---

## Output Format Constraint

Your entire response **MUST** be a single, valid JSON array of objects, with one object for each cluster.

You **MUST** include the all the provided marker genes (from input#4 "Cluster Marker Genes") and pathways (from input#5 "Cluster Enriched Pathways") for each cluster in the respective fields ("Top_marker_genes" and "Top_enriched_pathways").

**Do not** include any text, explanations, or markdown formatting before or after the JSON array. Adhere strictly to the schema shown in the example below.

```json
[
  {{
    "cluster_id": "0",
    "cell_type_hypotheses": "Macrophage",
    "justification": "This cluster is identified as Macrophage based on canonical markers (CD68, MRC1) and its adjacency to other myeloid clusters (Clusters 2, 3). Strong enrichment for 'HALLMARK_TNFA_SIGNALING_VIA_NFKB' and 'HALLMARK_INFLAMMATORY_RESPONSE' pathways indicates an activated state, while the core identity remains 'Macrophage'.",
    "key_markers_cited": ["CD68", "MRC1"],
    "confidence": "High",
    "Top_marker_genes": [All marker genes provided for cluster 0, separated by comma],
    "Top_enriched_pathways": [All pathways provided for cluster 0, separated by comma]
  }},
  {{
    "cluster_id": "1",
    "cell_type_hypotheses": "CD4+ T Cell",
    "justification": "This cluster is identified as T Cell by the expression of CD3D and CD3E, and specifically as CD4+ T Cell by CD4. This identity is strongly confirmed by the 'T_CELL_RECEPTOR_SIGNALING_PATHWAY'.",
    "key_markers_cited": ["CD3D", "CD3E", "CD4"],
    "confidence": "High",
    "Top_marker_genes": [All marker genes provided for cluster 1, separated by comma],
    "Top_enriched_pathways": [All pathways provided for cluster 1, separated by comma]
  }},
  {{
    "cluster_id": "2",
    "cell_type_hypotheses": "KRT8+ Transitional Epithelial Cell",
    "justification": "This cluster is identified as the KRT8+ transitional state based on the strong expression of KRT8. This hypothesis is strongly supported by its adjacency to both the 'Alveolar Type 2 Cell' cluster (Cluster 5) and the 'Alveolar Type 1 Cell' cluster (Cluster 9), confirming its role as a transitional population between them.",
    "key_markers_cited": ["KRT8"],
    "confidence": "High",
    "Top_marker_genes": [All marker genes provided for cluster 2, separated by comma],
    "Top_enriched_pathways": [All pathways provided for cluster 2, separated by comma]
  }},
  {{
    "cluster_id": "3",
    "cell_type_hypotheses": "Unknown",
    "justification": "This cluster expresses general myeloid markers (e.g., CD68) but lacks specific markers for the only candidate, 'Alveolar Macrophage'. This cluster likely represents the parent 'Macrophage' population, which is not an available annotation option per the Single-Candidate Exception.",
    "key_markers_cited": ["CD68"],
    "confidence": "Medium",
    "Top_marker_genes": [All marker genes provided for cluster 3, separated by comma],
    "Top_enriched_pathways": [All pathways provided for cluster 3, separated by comma]
  }}
]
```
"""