import unittest
import sys
sys.path.append("..")
from BioAddons import SeqStruct


class SeqStructTestCase(unittest.TestCase):

    def test_convert_1(self):
        '''Test converting RNA sequence and secondary structure sequence
        to RNAContextualSequenceSecondaryStructure sequence
        '''
        target = SeqStruct.SeqStruct.convert('AG', 'BE')
        expected = 'KA'
        self.assertEqual(target, expected)


    def test_convert_2(self):
        '''Test converting RNA sequence and secondary structure sequence
        to RNAContextualSequenceSecondaryStructure sequence
        '''
        target = SeqStruct.SeqStruct.convert('GCAUG', 'EEEEE')
        expected = 'AVHOA'
        self.assertEqual(target, expected)

if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
