import unittest
import sys
sys.path.append("..")
from BioAddons.Alphabet import ContextualSequenceSecondaryStructure as RNASS


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
        target = RNASS.convert('U', 'B')
        expected = 'R'
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

if __name__ == '__main__':
    print sys.argv[0]
    unittest.main()
