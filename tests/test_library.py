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

##########################################################################
# Functions and classes
##########################################################################

def gen_lib_pattern(size, pattern):
    """
    Generate library based on a given pattern

    :param size: size of the library
    :param pattern: Pattern to create a library (i.e. XCXXXCX)
    
    :return statistics from the generated population
    """
    
    libtest=Library(population_size=size, mode='exploration', pattern=pattern, add_phys_chem=False, mw_neigh=2, verbose=False)
    
    len_pop= len(libtest.population)
    
    list_pat = list(pattern)
    counter_hits_1=0
    counter_hits_2=0
    
    for i,char in enumerate(list_pat):
        if char != 'X':
            if char == libtest.population[0][i]:
                counter_hits_1+=1
            if char == libtest.population[-1][i]:
                counter_hits_2+=1
                
    return len_pop, counter_hits_1, counter_hits_2

##########################################################################

def gen_lib_random(size_pop, size_pep):
    """
    Generate library using random content

    :param size_pop: size of the library
    :param size_pep: peptide length. In the function multiple sizes are contemplated
    
    :return statistics from the generated population
    """

    libtest=Library(population_size=size_pop, mode='exploration', add_phys_chem=False, mw_neigh=2, min_pep_size=size_pep, max_pep_size=size_pep, verbose=False)
    
    len_pop= len(libtest.population)
    
    size1=len(libtest.population[0])
    size2=len(libtest.population[-1])

    return len_pop, size1, size2

##########################################################################

def gen_lib_seed_scan(size, seed):
    """
    Generate library bsaed on a seed. Some internal positions are changed

    :param size: size of the library
    :param seed: sequence used as a template
    
    :return statistics from the generated population
    """

    libtest=Library(population_size=size, mode='scanning', mode_scan='all', seeds=[seed], pairs=[(2,'A'),(3,'L'),(5,'R'),(7,'M')], positions=[2,3,5,7], from_child=True, verbose=False)
    
    len_pop= len(libtest.population)
    
    pos=[0,3,5,7]
    selected_aa1=[]
    selected_aa2=[]
    for p in pos:
        selected_aa1.append(libtest.population[0][p])
        selected_aa2.append(libtest.population[-1][p])
        
    return len_pop, selected_aa1, selected_aa2

##########################################################################

def gen_lib_seed_mult(size, seed):
    """
    Generate library bsaed on a seed. Some internal positions are changed

    :param size: size of the library
    :param seed: sequence used as a template
    
    :return statistics from the generated population
    """

    libtest=Library(population_size=size, mode='scanning', mode_scan='all', seeds=[seed], pairs=[(2,'A'),(3,'L'),(5,'R'),(7,'M')], positions=[2,3,5,7], single_mod=False, mw_neigh=4, no_priority=[3,5], verbose=False)
    
    len_pop= len(libtest.population)
    
    pos=[0,3,5,7]
    selected_aa1=[]
    selected_aa2=[]
    for p in pos:
        selected_aa1.append(libtest.population[0][p])
        selected_aa2.append(libtest.population[-1][p])
        
    return len_pop, selected_aa1, selected_aa2

########################################################################################

class TestLibrary(unittest.TestCase):

    def test_pattern(self):
        """
        Generate library based on a pattern
        """

        self.assertEqual(gen_lib_pattern(10, 'XCXXXCX'),(10, 2, 2))
        self.assertEqual(gen_lib_pattern(25, 'XXAXXDXXXAXX'),(25, 3, 3))
        self.assertEqual(gen_lib_pattern(100, 'ACDXXXXXXXEFG'),(100, 6, 6))
        
    ########################################################################################

    def test_random(self):
        """
        Generate a library randomly based on a desired size
        """

        self.assertEqual(gen_lib_random(10, 15), (10, 15, 15))
        self.assertEqual(gen_lib_random(25, 30), (25, 30, 30))
        self.assertEqual(gen_lib_random(100, 5), (100, 5, 5))

    ########################################################################################

    def test_seed_scan(self):
        """
        Generate a library using seeds and scanning methodologies
        """

        self.assertEqual(gen_lib_seed_scan(10, 'FNCREWCWN'), (10, ['F','R','W','W'], ['F','R','W','W']))
        self.assertEqual(gen_lib_seed_scan(25, 'NPVVHFFKNIVTP'), (25, ['N','V','F','K'], ['N','V','F','K']))
        self.assertEqual(gen_lib_seed_scan(100, 'LLSHYTSYAG'), (100, ['L','H','T','Y'], ['L','H','T','Y']))

    ########################################################################################

    def test_seed_mult(self):
        """
        Generate a library using seeds and scanning methodologies
        """

        self.assertEqual(gen_lib_seed_mult(10, 'FNCREWCWN'), (10, ['F','R','W','W'], ['F','R','W','W']))
        self.assertEqual(gen_lib_seed_mult(25, 'NPVVHFFKNIVTP'), (25, ['N','V','F','K'], ['N','V','F','K']))
        self.assertEqual(gen_lib_seed_mult(100, 'LLSHYTSYAG'), (100, ['L','H','T','Y'], ['L','H','T','Y']))
        
########################################################################################

if __name__ == "__main__":
    unittest.main()