# Copyright 2016 by Kevin Ha
"""Extended alphabet for RNA contextual secondary structure

"""
from Bio.Alphabet import SecondaryStructure


class ContextualSecondaryStructure(SecondaryStructure):
    """Extended alphabet for RNA contextual secondary structure"""

    letters = "EHTBLRM"
