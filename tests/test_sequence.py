"""
PepFunNN: Protocols for the analysis of peptides using cheminformatics and bioinformatics tools
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__email__ = "raoc@novonordisk.com"

########################################################################################
# Modules to import
########################################################################################

import unittest
from pepfunn.sequence import Sequence

##########################################################################
# Functions and classes
##########################################################################

def get_properties(sequence):
    """
    Calculate properties from the sequence

    :param sequence: Peptide sequence
    :return: net charge, average hydrophobicity, isoelectric point and molecular weigth
    """
    # Create object
    pep = Sequence(sequence)

    # Properties from the sequence
    netCharge=pep.netCharge
    avgHydro=pep.avg_hydro
    isoPoint=pep.isoelectric_point

    return netCharge,avgHydro,isoPoint

##########################################################################

def get_rules(sequence):
    """
    Calculate some empirical rules

    :param sequence: Peptide sequence
    :return: failed solubility rules and failed synthesis rules
    """
    # Create object
    pep = Sequence(sequence)

    # Empirical rules
    sol=pep.solubility_rules_failed
    syn=pep.synthesis_rules_failed

    return sol,syn

##########################################################################

class TestSequence(unittest.TestCase):

    def test_sequenceProperties(self):
        """
        Test the properties predicted for some peptides
        """

        self.assertEqual(get_properties('GAANDENY'), (-2.0008411239265196, -1.22, 4.0500284194946286))
        self.assertEqual(get_properties('EAPPSYAEV'), (-1.999412677161022, 1.1600000000000001, 4.0500284194946286))
        self.assertEqual(get_properties('SDVAFRGNLLD'), (-1.0055728771384773, 0.2000000000000003, 4.533500862121582))
        self.assertEqual(get_properties('RRNLKGLNLNLH'), (3.0896133167191, -4.58, 11.999967765808105))
        self.assertEqual(get_properties('GVLKEYGV'), (-0.001722146326308235, 2.2, 6.000685691833495))
        self.assertEqual(get_properties('MCLRMTAVM'), (0.9492040695067787, 2.39, 7.999833488464356))
        self.assertEqual(get_properties('EEFELLISNS'), (-2.99679181488375, 1.3299999999999996, 4.0500284194946286))
        self.assertEqual(get_properties('SQFDLSTRRLK'), (1.9935405781644697, -5.41, 10.834379386901855))
        self.assertEqual(get_properties('KLMFKTEGPDSD'), (-1.0081847499417282, -2.28, 4.780807304382324))

    ##########################################################################

    def test_sequenceRules(self):
        """
        Test similarity and some empirical rules for peptides
        """

        self.assertEqual(get_rules('GAANDENY'), (0, 0))
        self.assertEqual(get_rules('EAPPSYAEV'), (2, 2))
        self.assertEqual(get_rules('SDVAFRGNLLD'), (2, 0))
        self.assertEqual(get_rules('RRNLKGLNLNLH'), (4, 2))
        self.assertEqual(get_rules('GVLKEYGV'), (2, 0))
        self.assertEqual(get_rules('MCLRMTAVM'), (2, 2))
        self.assertEqual(get_rules('EEFELLISNS'), (3, 2))
        self.assertEqual(get_rules('SQFDLSTRRLK'), (3, 0))
        self.assertEqual(get_rules('KLMFKTEGPDSD'), (3, 1))

##########################################################################

if __name__ == "__main__":
    unittest.main()
