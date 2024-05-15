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
from pepfunn.library import Library
from pepfunn.clustering import simClustering
from pepfunn.clustering import propClustering
from pepfunn.clustering import sequence_clustering

##########################################################################
# Functions and classes
##########################################################################

def gen_sim_cluster(seed, size):
    """
    Calculate clusters based on peptide similarity
    
    :param seed: Sequence used to build the library
    :param size: size of the library
    :return statistics from the clusters
    """
    libtest=Library(population_size=size, mode='scanning', mode_scan='all', seeds=[seed], pairs=[(1,'A'),(4,'L'),(6,'R'),(8,'M')], positions=[2,3,5,7], from_child=True, verbose=False)
    sequences=libtest.population

    # Create cluster object
    clust = simClustering(sequences=sequences)
    clust.run_clustering()
    
    # Calculate neighbors
    neighbors=clust.get_sim_reference('mol1')
    num_neigh=len(neighbors)
                    
    return num_neigh

##########################################################################

def gen_prop_cluster(seed, size):
    """
    Calculate clusters based on peptide properties
    
    :param seed: Sequence used to build the library
    :param size: size of the library
    :return statistics from the clusters
    """
    libtest=Library(population_size=size, mode='scanning', mode_scan='all', seeds=[seed], pairs=[(1,'A'),(4,'L'),(6,'R'),(8,'M')], positions=[2,3,5,7], from_child=True, verbose=False)
    sequences=libtest.population

    # Create cluster object
    pc = propClustering(sequences=sequences)
    df=pc.calc_distances(['mol1'])
    num_rows=df.shape[0]
                    
    return num_rows

##########################################################################

def gen_seq_cluster(seed, size):
    """
    Calculate clusters based on sequences
    
    :param seed: Sequence used to build the library
    :param size: size of the library
    :return statistics from the clusters
    """
    libtest=Library(population_size=size, mode='scanning', mode_scan='all', seeds=[seed], pairs=[(1,'A'),(4,'L'),(6,'R'),(8,'M')], positions=[2,3,5,7], from_child=True, verbose=False)
    sequences=libtest.population

    # Create cluster object
    clusters, clus_ids, isolated, iso_ids, df_clus = sequence_clustering(sequences=sequences, threshold=0.8)
    flag=False
    if clusters:
        flag=True
                    
    return flag

########################################################################################

class TestClustering(unittest.TestCase):
    
    ########################################################################################
    def test_sim_clustering(self):
        """
        Clustering using similarity
        """

        self.assertEqual(gen_sim_cluster('FNCREWCWN', 100),100)
        self.assertEqual(gen_sim_cluster('NPVVHFFKNIVTP', 75),75)
        self.assertEqual(gen_sim_cluster('LLSHYTSYAG', 50),50)
    
    ########################################################################################
    def test_prop_clustering(self):
        """
        Clustering using properties
        """

        self.assertEqual(gen_prop_cluster('FNCREWCWN', 100),100)
        self.assertEqual(gen_prop_cluster('NPVVHFFKNIVTP', 75),75)
        self.assertEqual(gen_prop_cluster('LLSHYTSYAG', 50),50)
        
    ########################################################################################
    
    def test_seq_clustering(self):
        """
        Clustering using sequence
        """

        self.assertEqual(gen_seq_cluster('FNCREWCWN', 100),True)
        self.assertEqual(gen_seq_cluster('NPVVHFFKNIVTP', 75),True)
        self.assertEqual(gen_seq_cluster('LLSHYTSYAG', 50),True)
        
########################################################################################

if __name__ == "__main__":
    unittest.main()