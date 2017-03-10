# Copyright 2013 by Michiel de Hoon.  All rights reserved.
# Copyright 2016 modified by Kevin Ha
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#



from Bio.motifs import matrix
from Bio.Alphabet import NucleotideAlphabet
import platform

class ExtendedPositionSpecificScoringMatrix(matrix.PositionSpecificScoringMatrix):
    """This new class inherits Bio.motifs.matrix.PositionSpecificScoringMatrix.
    It has been modified to support any kind of Alphabet. This allows us to
    perform motif scans on RNA sequence as well as RNA secondary structure.

    The main change is the fact that the 'ACGT' hard-coding has been replaced
    with whatever letters are in the Alphabet of the matrix. This seems to be
    sufficient enough for our purposes.
    """

    def _py_calculate(self, sequence, m, n):
        """Handles the default calcuate() method in Python.
        Moved from _calculate in the except clause below.
        """

        # The C code handles mixed case so Python version must too:
        sequence = sequence.upper()
        scores = []
        for i in range(n - m + 1):
            score = 0.0
            for position in range(m):
                letter = sequence[i + position]
                try:
                    score += self[letter][position]
                except KeyError:
                    score = float("nan")
                    break
            scores.append(score)
        return scores

    # Make sure that we use C-accelerated PWM calculations if running under CPython.
    # Fall back to the slower Python implementation if Jython or IronPython.
    try:
        from . import _pwm
        def _calculate(self, sequence, m, n):

            # Only RNA and DNA is supported right now. If sequence is
            # secondary structure, then use Python implementation.
            if not isinstance(self.alphabet, NucleotideAlphabet):
                return self._py_calculate(sequence, m, n)

            letters = ''.join(sorted(self.alphabet.letters))
            logodds = [[self[letter][i] for letter in letters]
                            for i in range(m)]
            return self._pwm.calculate(sequence, logodds)
    except ImportError:
        if platform.python_implementation() == 'CPython':
            raise
        else:
            def _calculate(self, sequence, m, n):
                return self._py_calculate(sequence, m, n)


    def calculate(self, sequence):

        # TODO - Force uppercase here and optimise switch statement in C
        # by assuming upper case?
        sequence = str(sequence)
        m = self.length
        n = len(sequence)

        scores = self._calculate(sequence, m, n)

        if len(scores) == 1:
            return scores[0]
        else:
            return scores
