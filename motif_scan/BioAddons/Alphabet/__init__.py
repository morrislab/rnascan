# Copyright 2016 by Kevin Ha
"""Extended alphabet for RNA contextual secondary structure

"""
from Bio.Alphabet import SecondaryStructure


class ContextualSecondaryStructure(SecondaryStructure):
    """Extended alphabet for RNA contextual secondary structure"""

    letters = "EHTBLRM"


class ContextualSequenceSecondaryStructure(SecondaryStructure):
    """Extended alphabet for RNA contextual secondary structure and
    sequence combined"""

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ[]"

    # conversion = {}
    # _stack = list(letters)
    # _stack.reverse()
    # for i in IUPAC.IUPACUnambiguousRNA.letters:
    #     conversion[i] = {}
    #     for j in RNAContextualSecondaryStructure.letters:
    #         conversion[i][j] = _stack.pop()
    conversion = {'A':
                      {'B': 'K',
                       'E': 'H',
                       'H': 'I',
                       'L': 'L',
                       'M': 'N',
                       'R': 'M',
                       'T': 'J'},
                  'C':
                      {'B': 'Y',
                       'E': 'V',
                       'H': 'W',
                       'L': 'Z',
                       'M': ']',
                       'R': '[',
                       'T': 'X'},
                  'G':
                      {'B': 'D',
                       'E': 'A',
                       'H': 'B',
                       'L': 'E',
                       'M': 'G',
                       'R': 'F',
                       'T': 'C'},
                  'U':
                      {'B': 'R',
                       'E': 'O',
                       'H': 'P',
                       'L': 'S',
                       'M': 'U',
                       'R': 'T',
                       'T': 'Q'}
                 }

    back_conversion = {'A': ('G', 'E'),
                       'B': ('G', 'H'),
                       'C': ('G', 'T'),
                       'D': ('G', 'B'),
                       'E': ('G', 'L'),
                       'F': ('G', 'R'),
                       'G': ('G', 'M'),
                       'H': ('A', 'E'),
                       'I': ('A', 'H'),
                       'J': ('A', 'T'),
                       'K': ('A', 'B'),
                       'L': ('A', 'L'),
                       'M': ('A', 'R'),
                       'N': ('A', 'M'),
                       'O': ('U', 'E'),
                       'P': ('U', 'H'),
                       'Q': ('U', 'T'),
                       'R': ('U', 'B'),
                       'S': ('U', 'L'),
                       'T': ('U', 'R'),
                       'U': ('U', 'M'),
                       'V': ('C', 'E'),
                       'W': ('C', 'H'),
                       'X': ('C', 'T'),
                       'Y': ('C', 'B'),
                       'Z': ('C', 'L'),
                       '[': ('C', 'R'),
                       ']': ('C', 'M')}

    @classmethod
    def convert(cls, seqletter, structletter):
        """Return SeqStruct letter for given seqletter and structletter
        """
        return cls.conversion[seqletter][structletter]

    @classmethod
    def reverse_convert(cls, seqstructletter):
        """Return the RNA letter and ContextualSecondaryStructure
        letter as a tuple given a seqstructletter
        """
        return cls.back_conversion[seqstructletter]
