import unittest
import sys
import types
sys.path.append("..")
import motif_scan as ms
from BioAddons.Alphabet import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


class ParseSequencesTestCase(unittest.TestCase):

    def setUp(self):
        self.fastas = ['test.fa']

    def test_parse_sequences_1(self):
        '''Test parsing of test sequences
        '''
        target = ms.parse_sequences(self.fastas, IUPAC.IUPACUnambiguousRNA())
        self.assertIsInstance(target, types.GeneratorType)

        for item in target:
            self.assertIsInstance(item, SeqRecord)


class ComputeBackgroundTestCase(unittest.TestCase):

    def setUp(self):
        self.fastas = ['test.fa']

    def test_compute_background_1(self):
        target = ms.compute_background(self.fastas,
                                       IUPAC.IUPACUnambiguousRNA(),
                                       verbose=False)
        expected = {'A': 0.1944,
                    'C': 0.1388,
                    'U': 0.5277,
                    'G': 0.1388}

        for key,value in expected.iteritems():
            self.assertAlmostEqual(target[key], value, 3)


if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
