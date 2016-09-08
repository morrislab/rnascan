# Copyright 2000-2002 Brad Chapman.
# Copyright 2004-2005 by M de Hoon.
# Copyright 2007-2015 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
# Modified Copyright 2016 by Kevin Ha
"""This class inherits Bio.Seq that adds functionality for handling
ContextualSequenceSecondaryStructure alphabets.

Specifically, will take a RNA sequence and contextual secondary structure
sequence and convert it to a unified ContextualSequenceSecondaryStructure
alphabet.
"""
from Alphabet import ContextualSequenceSecondaryStructure as RNASS
from Bio.Seq import Seq


class SeqStruct(Seq):
    """A read-only Sequence object that extends Bio.Seq

    Adds extra function for converting RNA sequence and contextual secondary
    structure sequence into a ContextualSequenceSecondaryStructure sequence
    """
    def __init__(self, seq, struct):
        # Convert sequence and struct sequences
        newseq = SeqStruct.convert(seq, struct)
        super(SeqStruct, self).__init__(newseq, RNASS)


    @staticmethod
    def convert(seq, struct):
        """Convert a seq and struct SeqRecord to a new SeqRecord with
        alphabet ContextualSequenceSecondaryStructure
        """
        if len(seq) != len(struct):
            raise ValueError(('Sequence and structure records have'
                ' different lengths'))

        seqstruct_sequence = ''
        for i,j in zip(seq, struct):
            seqstruct_sequence += RNASS.convert(i, j)

        return seqstruct_sequence

if __name__ == "__main__":
    s = SeqStruct('AGC', 'BBB')
    print s
