import unittest
import sys
sys.path.append("..")
from motif_scan.BioAddons.Alphabet import \
    ContextualSequenceSecondaryStructure as RNASS


class RNAContextualSequenceSecondaryStructureTestCase(unittest.TestCase):

    def test_convert_1(self):
        '''Test converting RNA letter and structural letter to
        combined letter
        '''
        target = RNASS.convert('A', 'E')
        expected = 'H'
        self.assertEqual(target, expected)

    def test_convert_2(self):
        '''Test converting RNA letter and structural letter to
        combined letter
        '''
        target = RNASS.convert('C', 'M')
        expected = ']'
        self.assertEqual(target, expected)

    def test_convert_fail_1(self):
        '''Test converting unknown letter and structural letter to
        combined letter
        '''
        with self.assertRaises(KeyError):
            RNASS.convert('D', 'B')

    def test_convert_fail_2(self):
        '''Test converting RNA letter and unknown structural letter to
        combined letter
        '''
        with self.assertRaises(KeyError):
            RNASS.convert('A', 'Z')

    def test_convert_fail_3(self):
        '''Test converting unknown RNA letter and unknown structural letter to
        combined letter
        '''
        with self.assertRaises(KeyError):
            RNASS.convert('Q', 'Z')

    def test_backconvert_1(self):
        '''Test back conversion -> given SeqStruct letter, return
        RNA letter and ContextualSecondaryStructure letter
        '''
        target = RNASS.reverse_convert('A')
        expected = ('G', 'E')
        self.assertEqual(target, expected)

    def test_backconvert_2(self):
        '''Test back conversion -> given SeqStruct letter, return
        RNA letter and ContextualSecondaryStructure letter
        '''
        target = RNASS.reverse_convert('[')
        expected = ('C', 'R')
        self.assertEqual(target, expected)

    def test_backconvert_fail_1(self):
        '''Test back conversion -> given unknown letter, return
        KeyError failure
        '''
        with self.assertRaises(KeyError):
            RNASS.reverse_convert('k')

if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
