import unittest
import sys
sys.path.append("..")
from BioAddons.SeqStruct import SeqStruct
from BioAddons.Alphabet import *


class SeqStructTestCase(unittest.TestCase):

    def test_convert_1(self):
        '''Test converting RNA sequence and secondary structure sequence
        to ContextualSequenceSecondaryStructure sequence
        '''
        target = SeqStruct.convert('AG', 'BE')
        expected = 'KA'
        self.assertEqual(target, expected)


    def test_convert_2(self):
        '''Test converting RNA sequence and secondary structure sequence
        to ContextualSequenceSecondaryStructure sequence
        '''
        target = SeqStruct.convert('GCAUG', 'EEEEE')
        expected = 'AVHOA'
        self.assertEqual(target, expected)


    def test_create_seqstruct_1(self):
        '''Test instantiation of a SeqStruct object
        '''
        target = SeqStruct('AG', 'BE')
        self.assertIsInstance(target, SeqStruct)
        self.assertIsInstance(target.alphabet,
                              ContextualSequenceSecondaryStructure)


if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
