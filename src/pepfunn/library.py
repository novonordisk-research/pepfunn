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

# System
import random
import math
import numpy as np
import multiprocessing as mp
from random import randint
import itertools
import enum
import sys
import os
import pickle
import pandas as pd
from statistics import mean

# RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors

# Internal
from pepfunn.sequence import Sequence

########################################################################################
# Classes and Functions
########################################################################################

class LibraryConstants(enum.auto):
    """
    A class to hold defaults values
    """
    AA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    CLASSES = ['CHARGE', 'POLARITY', 'SMALL']
    PROPERTIES={'A': ['NEUTRAL', 'HYDROPHOBIC', 'SMALL'],
                'R': ['POSITIVE', 'CHARGED', 'LARGE'],
                'N': ['NEUTRAL', 'POLAR', 'MEDIUM'],
                'D': ['NEGATIVE', 'CHARGED', 'MEDIUM'],
                'C': ['NEUTRAL', 'HYDROPHOBIC', 'MEDIUM'],
                'Q': ['NEUTRAL', 'POLAR', 'LARGE'],
                'E': ['NEGATIVE', 'CHARGED', 'LARGE'],
                'G': ['NEUTRAL', 'HYDROPHOBIC', 'SMALL'],
                'H': ['POSITIVE', 'CHARGED', 'LARGE'],
                'I': ['NEUTRAL', 'HYDROPHOBIC', 'MEDIUM'],
                'L': ['NEUTRAL', 'HYDROPHOBIC', 'MEDIUM'],
                'K': ['POSITIVE', 'CHARGED', 'LARGE'],
                'M': ['NEUTRAL', 'HYDROPHOBIC', 'LARGE'],
                'F': ['NEUTRAL', 'HYDROPHOBIC', 'LARGE'],
                'P': ['NEUTRAL', 'HYDROPHOBIC', 'MEDIUM'],
                'S': ['NEUTRAL', 'POLAR', 'SMALL'],
                'T': ['NEUTRAL', 'POLAR', 'MEDIUM'],
                'W': ['NEUTRAL', 'HYDROPHOBIC', 'LARGE'],
                'Y': ['NEUTRAL', 'POLAR', 'LARGE'],
                'V': ['NEUTRAL', 'HYDROPHOBIC', 'MEDIUM']}

# End of Library class-related constants definition.
############################################################

class Library:

    """
    Class to build peptide libraries based on some requirements
    """
    ########################################################################################
    def __init__(self, seeds=[], population_size=100, 
            mode='exploration', pattern='', single_mod=True,
            list_aa=[], min_pep_size=3, max_pep_size=10, add_phys_chem=False, mw_neigh=10, nb_number=2,
            mode_scan='random', positions=[], pairs=[], add_scan=False, ref_scan=['A'], no_priority=[],
            add_prop_analog=False, perc_limit=0.1, from_child=False, verbose=True):

        """
        Definition of the parameters
        """
        self.population_size = population_size
        if not list_aa:
            self.list_aa = LibraryConstants.AA
        else:
            self.list_aa = list_aa
        self.total_aa = len(self.list_aa)
        self.max_pep_size=max_pep_size
        self.min_pep_size=min_pep_size
        self.positions=positions
        self.pairs=pairs
        self.mode=mode
        self.pattern=pattern
        self.seeds = seeds
        self.add_phys_chem = add_phys_chem
        self.mw_neigh=mw_neigh
        self.add_scan=add_scan
        self.ref_scan=ref_scan
        self.add_prop_analog=add_prop_analog
        self.from_child=from_child
        self.mode_scan=mode_scan
        self.nb_number=nb_number
        self.no_priority=no_priority
        self.perc_limit=perc_limit
        self.verbose=verbose

        if self.mode not in ['exploration', 'scanning']:
            raise ValueError("The mode should be exploration or scanning. Please correct")
            sys.exit(1)

        if self.mode in ['scanning'] and not self.seeds:
            raise ValueError("For scanning mode it is required a list of seeds. Please correct")
            sys.exit(1)
        
        if self.mode=='exploration':
            if not self.pattern:
                self.population = self.generate_random_population()
            else:
                self.population = self.generate_pattern_population()

        if self.mode=='scanning':
            if single_mod:
                self.population=self.generate_single_population()
            else:
                self.population=self.generate_multiple_population()
    
    ########################################################################################
    def check_phys_chem(self, gene):
        """
        Function to check some phys chem properties to filter the library

        :param gene: peptide proposed for the library
        """
        # Create PepFuNN object
        pep = Sequence(gene)
        flag_enter=1

        # Check empirical solubility rules
        if pep.solubility_rules_failed > 1:
            flag_enter=0
        
        # Check empirical synthesis rules
        if pep.synthesis_rules_failed > 1:
            flag_enter=0
        
        # Check a range of logP values. To define depending on the project
        if pep.mol_logp > 3 or pep.mol_logp < -3:
            flag_enter=0
        
        return flag_enter
    
    ########################################################################################
    def check_mw(self, population, mw, weights, pop_counts, dist_threshold=0.2):
        """
        Function to have control on the number of neighbors with similar molecular weight
        """
        flag_enter=1
        counter=0
        # Compare against the population of peptides
        for peptide in population:
            distance=abs(mw-weights[peptide])

            # Check if the difference in MW is below the threshold
            if distance < dist_threshold:
                counter+=1
                if peptide not in pop_counts:
                    pop_counts[peptide]=1
                else:
                    pop_counts[peptide]+=1
                if pop_counts[peptide] > self.mw_neigh:
                    flag_enter=0
                    break

            # If the number of neighbors is equal to the limit the gene will not be included
            if counter==self.mw_neigh:
                flag_enter=0
                break
        
        return flag_enter,counter,pop_counts
    
    ########################################################################################
    def find_prop_analog(self, new_aa):
        """
        Function to select randomly an amino acid analog from a list
        """
        class_aa=random.choice(LibraryConstants.CLASSES)
        index_aa=LibraryConstants.CLASSES.index(class_aa)
        property=LibraryConstants.PROPERTIES[new_aa][index_aa]

        temp_property=[]
        for key_aa in LibraryConstants.PROPERTIES:
            if LibraryConstants.PROPERTIES[key_aa][index_aa] == property:
                temp_property.append(key_aa)
        
        new_aa=random.choice(temp_property)
        while new_aa not in self.list_aa:
            new_aa=random.choice(temp_property)
        
        return new_aa
    
    ########################################################################################
    def generate_pattern_population(self):
        """"
        Function to generate a population based on a pattern
        """
        population = []
        pop_counts = {}
        weights = {}

        # Iterate until the population is completed
        while len(population) < self.population_size:
            if self.verbose:
                print(f'Population length: {len(population)}')
            gene = ""
            for aa in self.pattern:
                if aa == 'X':
                    gene += random.choice(self.list_aa)
                else:
                    gene += aa

            # Start the verification
            if gene not in population:
                flag_enter=1
                helm = ".".join(list(gene))
                mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
                mw = Descriptors.MolWt(mol)
                weights[gene]=mw

                # Include phys chem filters 
                if self.add_phys_chem:
                    flag_enter=self.check_phys_chem(gene)

                # Part to check how man
                if flag_enter==1:
                    flag_enter,counter,pop_counts=self.check_mw(population, mw, weights, pop_counts)
                
                # If all the conditions are achieved, please enter the peptide
                if flag_enter==1:
                    population.append(gene)
                    pop_counts[gene]=counter
        
        return population
    
    ########################################################################################
    def generate_random_population(self):
        """"
        Function to generate a population based on a pattern
        """
        population = []
        pop_counts = {}
        weights = {}

        # Iterate until the population is completed
        while len(population) < self.population_size:
            if self.verbose:
                print(f'Population length: {len(population)}')
            gene_size = round(random.uniform(self.min_pep_size, self.max_pep_size))
            gene = ""
            for j in range(gene_size):
                gene += random.choice(self.list_aa)
            
            # Start the verification
            if gene not in population:
                flag_enter=1
                helm = ".".join(list(gene))
                mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
                mw = Descriptors.MolWt(mol)
                weights[gene]=mw

                # Include phys chem filters 
                if self.add_phys_chem:
                    flag_enter=self.check_phys_chem(gene)

                # Part to check how man
                if flag_enter==1:
                    flag_enter,counter,pop_counts=self.check_mw(population, mw, weights, pop_counts)
                
                # If all the conditions are achieved, please enter the peptide
                if flag_enter==1:
                    population.append(gene)
                    pop_counts[gene]=counter
        
        return population

    ########################################################################################
    def generate_single_population(self): 
        """
        Function to generate a library based single mutations. It can be random, basic scanning or by matched pairs, and can use child peptides as seeds
        """

        if self.mode_scan not in ['random', 'pairs', 'all']:
            raise ValueError("The mode should be all, random, pairs or scan. Please correct")
            sys.exit(1)

        if self.mode_scan in ['all', 'random'] and not self.positions:
            raise ValueError("For random or full exploration you need to provide a list with the positions. Please correct")
            sys.exit(1)
        
        if self.mode_scan in ['all', 'pairs'] and not self.pairs:
            raise ValueError("For pairs or full exploration you need to provide a list with the pairs, including the position and the aa to mutate to. Please correct")
            sys.exit(1)
        
        population = []
        pop_counts = {}
        weights = {}
        
        temp_pop=self.seeds

        # Add scan versions of the seeds
        if self.add_scan:
            scanned=[]
            for ref_aa in self.ref_scan:
                for tseq in self.seeds:
                    list_amino = list(tseq)
                    for z,amino in enumerate(tseq):
                        list_amino[z]=ref_aa
                        scanned.append(''.join(list_amino))
                        list_amino[z]=amino
            temp_pop+=scanned

        # Iterate over the population
        while len(population) < self.population_size:
            if self.verbose:
                print(f'Population length: {len(population)}')

            seed = random.choice(temp_pop)
            all_genes=[]
            if self.mode_scan in ['all', 'random']:
                new_aa=random.choice(self.list_aa)
                pos=random.choice(self.positions)
                list_seed=list(seed)
                list_seed[pos-1]=new_aa
                gene=''.join(list_seed)
                if self.mode_scan == 'all':
                    all_genes.append(gene)
            
            if self.mode_scan in ['all', 'pairs']:
                pair=random.choice(self.pairs)
                new_aa=pair[1]
                pos=pair[0]

                # Select property analogs
                if self.add_prop_analog:
                    new_aa=self.find_prop_analog(new_aa)

                list_seed=list(seed)
                list_seed[pos-1]=new_aa
                gene=''.join(list_seed)

                if self.mode_scan == 'all':
                    all_genes.append(gene)
            
            if self.mode_scan == 'all':
                gene=random.choice(all_genes)

            # Check if the gene is new
            if gene not in population:
                flag_enter=1
                helm = ".".join(list(gene))
                mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
                mw = Descriptors.MolWt(mol)
                weights[gene]=mw
                
                # Check phys chem
                if self.add_phys_chem:
                    flag_enter=self.check_phys_chem(gene)

                if flag_enter==1:
                    flag_enter,counter,pop_counts=self.check_mw(population, mw, weights, pop_counts)

                if flag_enter==1:
                    population.append(gene)
                    pop_counts[gene]=counter

                    # Activate tree-based
                    if self.from_child:
                        temp_pop.append(gene)

        return population

    ########################################################################################
    def generate_multiple_population(self): 
        """
        Function to generate a library based on initial seeds. It can be random, basic scanning or by matched pairs
        """

        if self.mode_scan not in ['random', 'pairs', 'all']:
            raise ValueError("The mode should be all, random, pairs or scan. Please correct")
            sys.exit(1)

        if self.mode_scan in ['all', 'random'] and not self.positions:
            raise ValueError("For random or full exploration you need to provide a list with the positions. Please correct")
            sys.exit(1)
        
        if self.mode_scan in ['all', 'pairs'] and not self.pairs:
            raise ValueError("For pairs or full exploration you need to provide a list with the pairs, including the position and the aa to mutate to. Please correct")
            sys.exit(1)
        

        population = []
        pos_changed={}
        pop_counts = {}
        weights = {}
        
        # Iterate over the population
        while len(population) < self.population_size:
            if self.verbose:
                print(f'Population length: {len(population)}')
            seed = random.choice(self.seeds)

            #num = random.randint(1, nb_number)
            num=self.nb_number
            tried_mutations=[]
            for attempt in range(0,num):

                all_genes=[]
                if self.mode_scan in ['all', 'random']:
                    new_aa=random.choice(self.list_aa)
                    pos=random.choice(self.positions)
                    list_seed=list(seed)
                    list_seed[pos-1]=new_aa
                    gene=''.join(list_seed)
                    if self.mode_scan == 'all':
                        all_genes.append(gene)
                
                if self.mode_scan in ['all', 'pairs']:

                    flag_pair=0
                    while flag_pair==0:
                        pair=random.choice(self.pairs)
                        lenPair=len(pair)
                        pos_pair=random.randint(1, lenPair-1)
                        new_aa=pair[pos_pair]
                        pos=pair[0]

                        if pos not in tried_mutations:
                            if pos in self.no_priority:
                                if pos in pos_changed:
                                    if pos_changed[pos]<=self.population_size * self.perc_limit:
                                        flag_pair=1
                                        pos_changed[pos]+=1
                                        tried_mutations.append(pos)
                                else:
                                    flag_pair=1
                                    pos_changed[pos]=1
                                    tried_mutations.append(pos)    
                            else:
                                flag_pair=1
                                if pos not in pos_changed:
                                    pos_changed[pos]=1
                                else:    
                                    pos_changed[pos]+=1
                                tried_mutations.append(pos)

                    # To select analogs with similar properties. Pending to select one class in particular
                    if self.add_prop_analog:
                        new_aa=self.find_prop_analog(new_aa)

                    list_seed=list(seed)
                    list_seed[pos-1]=new_aa
                    gene=''.join(list_seed)

                    if self.mode_scan == 'all':
                        all_genes.append(gene)
                
                if self.mode_scan == 'all':
                    gene=random.choice(all_genes)
                
                seed=gene

            if gene not in population:
                flag_enter=1
                helm = ".".join(list(gene))
                mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
                mw = Descriptors.MolWt(mol)
                weights[gene]=mw
                
                if self.add_phys_chem:
                    flag_enter=self.check_phys_chem(gene)

                if flag_enter==1:
                    flag_enter,counter,pop_counts=self.check_mw(population, mw, weights, pop_counts)

                if flag_enter==1:
                    population.append(gene)
                    pop_counts[gene]=counter

        return population
    
    ########################################################################################
    def filter_population(self, max_single_double=20, final_number=80):
        """"
        Function to filter the population to control the number of neighbors between the peptides

        :param max_single_double: Max number of single and double mutations between the peptides
        :param final_number: final number of peptides in the filtered list
        """
        population=self.population
        new_pop={}
        amount_pairs=[]
        for i,seq1 in enumerate(population):
            total_dismatches=[]
            for j,seq2 in enumerate(population):
                if i!=j:
                    dismatches=0
                    for k,aa in enumerate(seq1):
                        if aa!=seq2[k]:
                            dismatches+=1
                    total_dismatches.append(dismatches)
            
            # Control the number of mutation with max two changes
            if total_dismatches.count(2) + total_dismatches.count(1) > max_single_double:
                new_pop[seq1]=total_dismatches.count(2) + total_dismatches.count(1)
        
        # Sort and select the first number of elements
        sorted_dict = sorted(new_pop.items(), key=lambda item: item[1], reverse=True)
        final_pop=[]
        for i in range(0,final_number):
            try:
                final_pop.append(sorted_dict[i][0])
            except:
                raise ValueError("The required number is higher that the number of filtered peptide. Please reduce the number of candidates or decrease the max_single_double variable")
                sys.exit(1)
        
        return final_pop
   
############################################################
# End of library.py
############################################################