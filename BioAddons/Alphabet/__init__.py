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

    @classmethod
    def convert(cls, seqletter, structletter):
      """Return SeqStruct letter for given seqletter and structletter
      """
      return cls.conversion[seqletter][structletter]


