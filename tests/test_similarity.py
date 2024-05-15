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
from pepfunn.similarity import Alignment
from pepfunn.similarity import pepDescriptors
from pepfunn.similarity import simMonFP

##########################################################################
# Functions and classes
##########################################################################

def run_same_len(seq1, seq2):
    """
    Run alignment between two peptides and calculate score

    :param seq1: Peptide sequence 1
    :param seq2: Peptide sequence 2
    :return: score
    """
    
    
    score1 = Alignment.align_samelen_matrix(seq1, seq2)
    score2 = Alignment.align_samelen_local(seq1, seq2)
    score3 = Alignment.align_smiles(seq1, seq2)

    return score1, score2, score3

##########################################################################

def run_diff_len(seq1, seq2, mode):
    """
    Run alignment with two peptides of different size and calculate score

    :param seq1: Peptide sequence 1
    :param seq2: Peptide sequence 2
    :param mode: mode to weight the matches
    :return score
    """

    score, align1, align2 = Alignment.align_difflen_matrix(seq1, seq2, mode=mode)

    return score, align1, align2

##########################################################################

def run_similarity(seq1, seq2, mode='same'):
    """
    Get similarity between two peptide, no matter the size

    :param seq1: Peptide sequence 1
    :param seq2: Peptide sequence 2
    :param mode: same or diff, depending on the peptides lenght
    :return similarity value
    """

    sim = Alignment.similarity_pair(seq1, seq2, mode=mode)

    return sim

##########################################################################

def generate_descriptors(sequence):
    """
    Generate descriptors for a peptide (BILN format)

    :param sequence: Peptide sequence
    :return: First 5 descriptors
    """

    # Create the dictionary
    desc=pepDescriptors(sequence)
    descriptors=desc.moranCorrelation()

    desc_pairs=[]
    counter=1
    for key, val in descriptors.items():
        if counter > 5:
            break
        desc_pairs.append((key,val))
        counter+=1

    return desc_pairs

##########################################################################

def mon_similarity(seq1, seq2, add_freq=False):
    """
    Get similarity between two peptides using monomer fingerprint

    :param seq1: Peptide sequence 1
    :param seq2: Peptide sequence 2
    :param add_freq: add frequence of pairs of triplets in the fingerprints
    :return similarity value
    """

    sim = simMonFP(seq1, seq2, radius=2, nBits=1024, add_freq=add_freq)

    return sim

########################################################################################

class TestSimilarity(unittest.TestCase):

    def test_scores(self):
        """
        Run similarity analysis between 2 sequences (equal length)
        """

        self.assertEqual(run_same_len('AFTGYW','AGTGYL'),(11.815126050420169, 4, 0.4652777777777778))
        self.assertEqual(run_same_len('NPVVHFFKNIVTPRTPPPSQ', 'AAAAAFFKNIVAAAAAAAAA'),(6.087468571256227, 6, 0.34172661870503596))
        self.assertEqual(run_same_len('LLSHYTSY', 'LLSHYTSY'),(32.0, 8, 1.0))

    ########################################################################################

    def test_scores_diff(self):
        """
        Run alignment analysis between 2 sequences of different lenght and having NNAAs too
        """

        self.assertEqual(run_diff_len("W-W-S-E-V-N-R-A-E-F", "K-T-E-E-I-S-E-V-N-I-V-A-E-F", mode="unweighted"), (7.0, 'WW-----SEVNR--AEF', '--KTEEISEVN-IVAEF'))
        self.assertEqual(run_diff_len("K-Aib-M-P","S-A-Aib-P", mode="weighted"), (8.0, 'K--AibMP', '-SAAib-P'))
        self.assertEqual(run_diff_len("L-L-S-H-Y-T-S-Y", "A-F-T-G-Y-W", mode="weighted"), (9.487773487773488, 'LLSHYTS-Y-', '-A--FT-GYW'))

    ########################################################################################

    def test_similarity(self):
        """
        Run similarity between 2 sequences of different lenght and having NNAAs too
        """

        self.assertEqual(run_similarity('AFTGYW', 'AGTGYL', mode='same'), 0.49229691876750703)
        self.assertEqual(run_similarity('NPVVHFFKNIVTPRTPPPSQ', 'AAAAAFFKNIVAAAAAAAAA', mode='same'), 0.07609335714070284)
        self.assertEqual(run_similarity('LLSHYTSY', 'LLSHYTSY', mode='same'), 1.0)
        self.assertEqual(run_similarity("W-W-S-E-V-N-R-A-E-F", "K-T-E-E-I-S-E-V-N-I-V-A-E-F", mode='diff'), 0.5916079783099616)
        self.assertEqual(run_similarity("K-Aib-M-P","S-A-Aib-P", mode='diff'), 0.5)
        self.assertEqual(run_similarity("L-L-S-H-Y-T-S-Y", "A-F-T-G-Y-W", mode='diff'), 0.34236053607351363)

    ########################################################################################

    def test_descriptors(self):
        """
        Generate descriptors for a modified peptide. Monomers should be part of the dictionary
        """

        self.assertEqual(generate_descriptors("K-Aib-M-P-C-A"),[('nrot-lag1', -0.6235294117647058), ('nrot-lag2', 0.7941176470588235), ('nrot-lag3', -0.9215686274509803),('nrot-lag4', 0.47058823529411764), ('nrot-lag5', -1.2352941176470589)])
        self.assertEqual(generate_descriptors("S-Aib-A-L-L-R-E-Iva-P"),[('nrot-lag1', 0.4625000000000001), ('nrot-lag2', -0.20714285714285713), ('nrot-lag3', -0.75),('nrot-lag4', -0.7), ('nrot-lag5', -0.475)])
        self.assertEqual(generate_descriptors("A-A-L-Nva-C-W-E-Orn"),[('nrot-lag1', 0.5201109570041611), ('nrot-lag2', -0.0032362459546925555), ('nrot-lag3', -0.0058252427184466),('nrot-lag4', -0.009708737864077667), ('nrot-lag5', -0.43042071197411)])
    
    ########################################################################################

    def test_mon_similarity(self):
        """
        Run similarity between 2 sequences of different lenght and having NNAAs too
        """

        self.assertEqual(mon_similarity('AFTGYW', 'AGTGYL'), 0.3333333333333333)
        self.assertEqual(mon_similarity('NPVVHFFKNIVTPRTPPPSQ', 'AAAAAFFKNIVAAAAAAAAA', add_freq=True), 0.2037037037037037)
        self.assertEqual(mon_similarity('LLSHYTSY', 'LLSHYTSY'), 1.0)
        self.assertEqual(mon_similarity("W-W-S-E-V-N-R-A-E-F", "K-T-E-E-I-S-E-V-N-I-V-A-E-F", add_freq=True), 0.3023255813953488)
        self.assertEqual(mon_similarity("K-Aib-M-P","S-A-Aib-P"), 0.125)
        self.assertEqual(mon_similarity("L-L-S-H-Y-T-S-Y", "A-F-T-G-Y-W"), 0.06666666666666667)


########################################################################################

if __name__ == "__main__":
    unittest.main()