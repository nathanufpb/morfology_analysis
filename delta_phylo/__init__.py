"""
delta-phylo: Morphological phylogenetic analyses using DELTA format datasets.

A modern Python library for systematists, taxonomists, and evolutionary
biologists working with morphological characters and DELTA data.

Quick start::

    from delta_phylo import DeltaReader, MatrixBuilder, MorphologicalMatrix
    from delta_phylo import NexusWriter
    from delta_phylo import ParsimonyAnalysis, DistanceAnalysis
"""

from delta_phylo._version import __version__

# Core parser
from delta_phylo.parser import DeltaReader, Character, CharacterType, Taxon, CharacterState

# Matrix
from delta_phylo.matrix import MatrixBuilder, MorphologicalMatrix, MatrixCleaner, MatrixEncoder

# IO writers
from delta_phylo.io import NexusWriter, PhylipWriter, TNTWriter, CSVWriter

# Phylogenetic analyses
from delta_phylo.phylogeny import (
    ParsimonyAnalysis,
    DistanceAnalysis,
    gower_distance,
    LikelihoodAnalysis,
    mk_likelihood,
    TreeUtils,
)

# Utilities
from delta_phylo.utils import validate_delta_directory, validate_matrix, get_logger

__all__ = [
    "__version__",
    # parser
    "DeltaReader",
    "Character",
    "CharacterType",
    "Taxon",
    "CharacterState",
    # matrix
    "MatrixBuilder",
    "MorphologicalMatrix",
    "MatrixCleaner",
    "MatrixEncoder",
    # io
    "NexusWriter",
    "PhylipWriter",
    "TNTWriter",
    "CSVWriter",
    # phylogeny
    "ParsimonyAnalysis",
    "DistanceAnalysis",
    "gower_distance",
    "LikelihoodAnalysis",
    "mk_likelihood",
    "TreeUtils",
    # utils
    "validate_delta_directory",
    "validate_matrix",
    "get_logger",
]
