from pydantic import BaseModel, Field
from typing import List

class CellTypeHypothesis(BaseModel):
    cluster: str = Field(description="The cluster ID as provided in the input.")
    cell_type_hypothesis: str = Field(description="The hypothesized cell type for the cluster (no description included).")
    cell_type_description: str = Field(description="The concise description of the hypothesized cell type.")
    confidence: str = Field(description="Confidence level of the hypothesis: High, Medium, Low.")
    justification: str = Field(description="Detailed justification for the hypothesis.")
    split_status: str = Field(description="Status indicating if further subclustering is beneficial: Yes, No.")
    next_round_context: str = Field(description="Context to be used in next round annotation if needed")
    key_markers_cited: List[str]
    key_pathways_cited: List[str]
    note: str = Field(description="The log of the harmonization process.")

class CellTypeHypotheses(BaseModel):
    hypotheses: List[CellTypeHypothesis]

class CellTypeHarmonized(BaseModel):
    cluster: str = Field(description="The cluster ID as provided in the input.")
    harmonized_cell_type: str = Field(description="The harmonized cell type annotation for the cluster.")
    harmonization_log: str = Field(description="The log of the harmonization process.")

class CellTypeHarmonizations(BaseModel):
    harmonizations: List[CellTypeHarmonized]