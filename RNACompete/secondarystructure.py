"""
Extended alphabet for RNA contextual secondary structure

"""
from Bio import Alphabet


class RNAContextualSecondaryStructure(Alphabet.SingleLetterAlphabet):
    """Extended alphabet for RNA contextual secondary structure"""

    letters = "EHTBLRM"
