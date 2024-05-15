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
import warnings
from pepfunn.library import Library
from pepfunn.pairs import MatchedPairs
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

##########################################################################
# Functions and classes
##########################################################################

def gen_pairs_ref(seed, size):
    """
    Calculate pairs using a reference peptide
    
    :param seed: Sequence used to build the library
    :param size: size of the library
    :return number of pairs
    """
    
    libtest=Library(population_size=size, mode='scanning', mode_scan='all', seeds=[seed], pairs=[(2,'A'),(3,'L'),(5,'R'),(7,'M')], positions=[2,3,5,7], from_child=True, verbose=False)
    
    dict_comparison={'ID':[], 'Sequence':[], 'MW':[]}

    for i,seq in enumerate(libtest.population):
        name=f'TEST{i+1}'
        helm = ".".join(list(seq))
        mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
        mw = Descriptors.MolWt(mol)

        dict_comparison['ID'].append(name)
        dict_comparison['Sequence'].append(seq)
        dict_comparison['MW'].append(mw)

    df = pd.DataFrame(dict_comparison)
    
    id_column='ID'
    seq_column='Sequence'
    seq_ref=seed
    name_ref='TEST'
    property_columns=['MW']

    df_sequences=MatchedPairs.get_sequences(df, id_column, seq_column, seq_ref, name_ref=name_ref, threshold=0.4)
    
    df_pairs=MatchedPairs.get_pairs(df, df_sequences, property_columns, name_ref=name_ref, operation_columns='')
    
    num_rows = df_pairs.shape[0]
                
    return num_rows

##########################################################################

def gen_pairs_gap(seed, size):
    """
    Calculate pairs using gaps with or without reference peptide
    
    :param seed: Sequence used to build the library
    :param size: size of the library
    :return number of pairs with and without a reference
    """
    
    libtest=Library(population_size=size, mode='scanning', mode_scan='all', seeds=[seed], pairs=[(1,'A'),(4,'L'),(6,'R'),(8,'M')], positions=[2,3,5,7], from_child=True, verbose=False)
    
    dict_comparison={'ID':[], 'Sequence':[], 'MW':[]}

    for i,seq in enumerate(libtest.population):
        name=f'TEST{i+1}'
        helm = ".".join(list(seq))
        mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
        mw = Descriptors.MolWt(mol)

        dict_comparison['ID'].append(name)
        dict_comparison['Sequence'].append(seq)
        dict_comparison['MW'].append(mw)

    df = pd.DataFrame(dict_comparison)
    
    id_column='ID'
    seq_column='Sequence'
    seq_ref=seed
    property_columns=['MW']
    
    warnings.filterwarnings('ignore')
    
    df_sequences = MatchedPairs.get_sequences_simple(df, id_column, seq_column)
    df_pairs = MatchedPairs.get_pairs_simple(df, df_sequences, property_columns, operation_columns='')
    
    num_rows = df_pairs.shape[0]
                
    return num_rows

########################################################################################

class TestPairs(unittest.TestCase):

    def test_pairs_ref(self):
        """
        Analyze pairs using a reference peptide
        """

        self.assertEqual(gen_pairs_ref('FNCREWCWN', 100),9900)
        self.assertEqual(gen_pairs_ref('NPVVHFFKNIVTP', 75),5550)
        self.assertEqual(gen_pairs_ref('LLSHYTSYAG', 50),2450)
        
    ########################################################################################
    
    def test_pairs_gap(self):
        """
        Analyze pairs using simple alignments with gaps
        """

        self.assertEqual(gen_pairs_gap('FNCREWCWN', 50),2450)
        self.assertEqual(gen_pairs_gap('NPVVHFFKNIVTP', 35),1190)
        self.assertEqual(gen_pairs_gap('LLSHYTSYAG', 25),600)
        
########################################################################################

if __name__ == "__main__":
    unittest.main()