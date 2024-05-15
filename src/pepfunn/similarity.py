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

# External modules
import math
import pandas as pd
import re
import os
import pickle
import numpy as np
import warnings
import sys
from itertools import combinations

# BioPython
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors

# Internal
from pepfunn.sequence import Sequence
from pepfunn.sequence import SequenceConstants
from pepfunn.sequence import get_monomer_info
from pepfunn.sequence import split_outside

# Igraph
from igraph import Graph

########################################################################################
# Classes and functions
########################################################################################

class Alignment:

    """Class to allow the alignment of two sequences, including those with natural or unnatural amino acids.
       
       Note: For modified peptides the sequences should separate the residues by dashes, 
       and the unnatural should be part of the monomer dictionary
    """

    ############################################################################
    @staticmethod
    def get_matrix(path):
        """
        Read the matrix file to run the similarity analysis
        :param path: text stream with the matrix

        :return: matrix dictionary
        """
        matrix_temp = {}
        data = path.read().split('\n')
        for line in data:
            if line:
                fields = line.split()
                key = (fields[0], fields[1])
                matrix_temp[key] = float(fields[2])

        return matrix_temp

    ############################################################################  
    @staticmethod
    def align_samelen_matrix(pep_reference, pep_match, matrixInput=False):
        """
        Align position by position a peptide of reference with another one of the same length

        :param pep_reference: String with the sequence of a reference peptide
        :param pep_match: String with the sequence of a peptide that will be compared
        :param matrixInput: A similarity matrix calculated internally. This should follow
                            the format of the matrix available in data/matrix.txt

        :return score: A numerical value to rank the alignment
        """

        # Check len peptides
        if len(pep_reference)!=len(pep_match):
            warnings.warn(f"The peptides to compare should have the same length")
            sys.exit(1)

        if matrixInput:
            matrix = matrixInput
        else:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, SequenceConstants.def_matrix)
            with open(file_path, 'r') as handle:
                matrix = Alignment.get_matrix(handle)

        score_matrix = 0

        for i, p in enumerate(pep_reference):
            # Generate the tuples with the pair of amino acids to compare
            pair1 = (p, pep_match[i])
            pair2 = (pep_match[i], p)

            # Obtain the score from the matrix and sum up the values
            if pair1 in matrix:
                value = matrix[pair1]
            else:
                value = matrix[pair2]
            score_matrix += value
        
        return score_matrix

    ############################################################################
    @staticmethod
    def align_samelen_local(pep_reference, pep_match):
        """
        Align position by position a peptide of reference with another one

        :param pep_reference: String with the sequence of a reference peptide
        :param pep_match: String with the sequence of a peptide that will be compared

        :return match: Number of matches
        """
        # Check len peptides
        if len(pep_reference)!=len(pep_match):
            warnings.warn(f"The peptides to compare should have the same length")
            sys.exit(1)

        dismatch = 0

        for i, p in enumerate(pep_reference):
            # Generate the tuples with the pair of amino acids to compare
            if p != pep_match[i]:
                dismatch += 1

        match = len(pep_reference) - dismatch
        return match

    ############################################################################
    @staticmethod
    def align_smiles(pep_reference, pep_match):
        """
        Calculate similarity but using SMILES representations of the peptides

        :param pep_reference: String with the sequence of a reference peptide
        :param pep_match: String with the sequence of a peptide that will be compared

        :return smiles_similarity: based on Morgan Fingerprints and Tanimoto coefficient
        """

        # Generate molecule from sequence
        seq1=Sequence(pep_reference)

        mol1 = Chem.MolFromSmiles(seq1.smiles)
        mol1.SetProp("_Name", pep_reference)

        seq2=Sequence(pep_match)

        mol2 = Chem.MolFromSmiles(seq2.smiles)
        mol2.SetProp("_Name", pep_match)

        # Calculate the fingerprints and the similarity
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 4, 2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 4, 2048)

        smiles_similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        
        return smiles_similarity
    
    ############################################################################
    @staticmethod
    def align_smiles_mol(smiles1, smiles2):
        """
        Calculate similarity but using SMILES representations of the peptides

        :param smiles1: SMILES of molecule 1
        :param smiles2: SMILES of molecule 2

        :return smiles_similarity: based on Morgan Fingerprints and Tanimoto coefficient
        """

        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        # Calculate the fingerprints and the similarity
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 4, 2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 4, 2048)

        smiles_similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        
        return smiles_similarity

    ############################################################################
    @staticmethod
    def generate_matrix(threshold=60):
        """
        Class to generate a similarity matrix for a set of non-natural monomers

        :param threshold: Value to define a similarity threshold for the matrix
        :return: A pickle file with the generated matrix. This matrix is provided by default in the package
        """

        # Read the SDF file with the monomer information
        module_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(module_dir, SequenceConstants.def_path)
        file_path = os.path.join(data_dir, SequenceConstants.def_lib_filename)
        with open(file_path, 'rb') as handle:
            df = get_monomer_info(handle)

        # Extract the smiles from the dataframe
        totalNames=[]
        totalSmiles=[]
        for i, idx in enumerate(df.index):
            name = df.at[idx, 'm_abbr']
            smiles = df.at[idx, 'm_smiles']
            if smiles!='None':
                totalNames.append(name)
                totalSmiles.append(smiles)

        # Run pair alignments and create the matrix
        matrixP = {}
        print("Generating the matrix ...")
        for i,name1 in enumerate(totalNames):
            smiles1=totalSmiles[i]
            totalValues=[]
            for j,name2 in enumerate(totalNames):
                smiles2=totalSmiles[j]
                try:
                    mol1= Chem.MolFromSmiles(smiles1)
                    mol2 = Chem.MolFromSmiles(smiles2)

                    fp1 = rdMolDescriptors.GetMorganFingerprint(mol1,2,useChirality=True)
                    fp2 = rdMolDescriptors.GetMorganFingerprint(mol2,2,useChirality=True)
                    smiles_similarity=DataStructs.DiceSimilarity(fp1,fp2)

                    totalValues.append(smiles_similarity)
                except:
                    totalValues.append(0.0)

            for i,ele in enumerate(totalValues):
                name2=totalNames[i]
                if (name1,name2) not in matrixP and (name2,name1) not in matrixP:
                    matrixP[(name1,name2)]=10*ele-(threshold/10)


        # Output file in a local data folder
        print("The matrix with a threshold of {}% has been generated".format(threshold))
        data_dir = os.path.join(module_dir, SequenceConstants.def_path)
        file_path = os.path.join(data_dir, SequenceConstants.def_matrix)
        with open(file_path, 'w') as handle:
            for pair in matrixP:
                handle.write(f'{pair[0]} {pair[1]} {matrixP[pair]}\n')

    ########################################################################################
    @staticmethod
    def __convert_list(peptide):
        """
        Function to convert a peptide separated by dashes in a list

        :param peptide: Peptide sequence separated by dashes. This can include unnatural monomers
        :return: List with the peptide monomers
        """
        peptideList = []
        monomers = split_outside(peptide, SequenceConstants.monomer_join, '[]')
        for res in monomers:
            mon = re.sub("\(\d+,\d+\)", "", res)
            peptideList.append(str(mon))
        return peptideList

    ########################################################################################
    @staticmethod
    def align_difflen_matrix(peptide1, peptide2, mode='weighted', matrixInput=False):
        """
        Method to compare two sequences with non-natural amino acids. The input sequences should have the monomers separated by dashes.
        and the monomers should be part of the provided monomer dictionary.

        :param peptide1: Peptide input 1 with monomers separated by dashes.
        :param peptide2: Peptide input 2 with monomers separated by dashes)
        :param mode: Use the matrix or not in the alignment. Options between unweighted or weighted (default)
        :param matrix: Name of the matrix file that will be used to run the alignment

        :return: score of the alignment
        """

        # Check if the peptides are in BILN format.
        if "-" not in peptide1 or "-" not in peptide2:
            peptide1='-'.join(list(peptide1))
            peptide2='-'.join(list(peptide2))
        if mode not in ['unweighted','weighted']:
            warnings.warn("Please select one mode between weighted or unweighted (default)")
            sys.exit(1)

        # Convert the sequences to list
        peptideList1 = Alignment.__convert_list(peptide1)
        peptideList2 = Alignment.__convert_list(peptide2)

        # Without matrix
        if mode == "unweighted":
            #print("Only matching alignment")
            for a in pairwise2.align.globalxx(peptideList1, peptideList2, gap_char=["-"]):
                score=a.score
                start=a.start
                end=a.end
                
                formatted_alignment = format_alignment(*a).splitlines()
                # Prepare the alignment string
                
                # Extract only the aligned sequences with gaps, removing the leading labels
                seq1_aligned_with_gaps = formatted_alignment[0].replace(" ", "")
                seq2_aligned_with_gaps = formatted_alignment[2].replace(" ", "")
                break

        # With matrix
        elif mode == "weighted":

            # Read the matrix information
            if matrixInput:
                simMatrix = matrixInput
            else:
                module_dir = os.path.dirname(os.path.abspath(__file__))
                data_dir = os.path.join(module_dir, SequenceConstants.def_path)
                file_path = os.path.join(data_dir, SequenceConstants.def_matrix)
                with open(file_path, 'r') as handle:
                    simMatrix = Alignment.get_matrix(handle)

            for a in pairwise2.align.globaldx(peptideList1, peptideList2, simMatrix, gap_char=["-"]):
                score = a.score
                start = a.start
                end = a.end
                
                # Prepare the alignment string
                formatted_alignment = format_alignment(*a).splitlines()
                
                # Extract only the aligned sequences with gaps, removing the leading labels
                seq1_aligned_with_gaps = formatted_alignment[0].replace(" ", "")
                seq2_aligned_with_gaps = formatted_alignment[2].replace(" ", "")
                break
        
        return score, seq1_aligned_with_gaps, seq2_aligned_with_gaps

    ############################################################################
    @staticmethod
    def similarity_pair(peptide1, peptide2, mode='same', weight=True, matrixInput=False):
        """
        Function to calculate similarity between two peptide sequences

        :param peptide1: Sequence of one of the input peptides
        :param peptide2: Sequence of the second input peptide
        :param mode: Compare peptide of the same length or different length 

        :return sim_val: similarity based on the alignments between the peptides and themselves - Max value is 1
        """

        if mode not in ['same','diff']:
            warnings.warn("Please select one mode between same (default) or diff")
            sys.exit(1)

        if matrixInput:
            matrix = matrixInput
        else:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, SequenceConstants.def_matrix)
            with open(file_path, 'r') as handle:
                matrix = Alignment.get_matrix(handle)

        if mode=='same':
            # Alignment between peptide1 and peptide2
            score1_2 = Alignment.align_samelen_matrix(peptide1, peptide2, matrixInput=matrix)
        
            # Alignment between peptide1 with itself
            score1_1 = Alignment.align_samelen_matrix(peptide1, peptide1, matrixInput=matrix)
        
            # Alignment between peptide2 with itself
            score2_2 = Alignment.align_samelen_matrix(peptide2, peptide2, matrixInput=matrix)

            # Calculate similarity value
            sim_val = float(score1_2) / math.sqrt(float(score1_1 * score2_2))
        
        if mode=='diff':
            if weight:
                # Alignment between peptide1 and peptide2
                score1_2, seq1, seq2 = Alignment.align_difflen_matrix(peptide1, peptide2, matrixInput=matrix)
            
                # Alignment between peptide1 with itself
                score1_1, seq1, seq2 = Alignment.align_difflen_matrix(peptide1, peptide1, matrixInput=matrix)
            
                # Alignment between peptide2 with itself
                score2_2, seq1, seq2 = Alignment.align_difflen_matrix(peptide2, peptide2, matrixInput=matrix)

                # Calculate similarity value
                sim_val = float(score1_2) / math.sqrt(float(score1_1 * score2_2))
            else:
                score1_2, seq1, seq2 = Alignment.align_difflen_matrix(peptide1, peptide2, mode='unweighted', matrixInput=matrix)
                score1_1, seq1, seq2 = Alignment.align_difflen_matrix(peptide1, peptide1, mode='unweighted', matrixInput=matrix)
                score2_2, seq1, seq2 = Alignment.align_difflen_matrix(peptide2, peptide2, mode='unweighted', matrixInput=matrix)
                sim_val = float(score1_2) / math.sqrt(float(score1_1 * score2_2))


        # Return similarity
        return sim_val

    # End of Alignment class
    ############################################################

##########################################################################
# Additional function
########################################################################## 
def monomerFP(peptide, radius=2, nBits=1024, add_freq=False, prop_list=['heavy', 'nrot', 'hacc', 'hdon', 'nhet', 'tpsa', 'mw']):
    """
    Function to calculate a monomer-based fingerprint based on a list of monomer properties and a radius

    :param peptide: BILN representation of the molecule
    :param radius: Number of neighbors around each residue. Maximum 2
    :param nBits: Number of bits in the fingerprint
    :param add_freq: Flag to add the frequency of the fragment in the tag code
    :param prop_list: List of properties used to create the fingerprint bits
    """

    # Check if the peptides are in BILN format.
    if "-" not in peptide:
        peptide='-'.join(list(peptide))

    # Create the graph
    graph = graphPep(peptide)
    graph.createGraph()
    
    subfrags={}
    for i,v in enumerate(graph.g.vs):
        name=v["name"]
        if '-' in name:
            name='['+name+']'
        if name not in subfrags:
            subfrags[name]=1
        else:
            subfrags[name]+=1

        # get the neighbors of the current vertex
        seed=name
        for r in range(0, radius):
            if r==0:
                neighbors = graph.g.neighbors(v)
                seed_list=[]
                index_list=[i]
                for n in neighbors:
                    subname=graph.g.vs[n]["name"]
                    if '-' in subname:
                        subname='['+subname+']'
                    frag=seed+'-'+subname
                    if frag not in subfrags:
                        subfrags[frag]=1
                    else:
                        subfrags[frag]+=1
                    seed_list.append(frag)
                    index_list.append(n)
            if r==1:
                for j,n1 in enumerate(neighbors):
                    neighbors_2 = graph.g.neighbors(graph.g.vs[n1])
                    for n2 in neighbors_2:
                        if n2 not in index_list:
                            subsubname=graph.g.vs[n2]["name"]
                            if '-' in subsubname:
                                subsubname='['+subsubname+']'
                            frag=seed_list[j]+'-'+subsubname
                            if frag not in subfrags:
                                subfrags[frag]=1
                            else:
                                subfrags[frag]+=1

    # Iterate over the subfrags based on the radius
    frags={}
    for fragment in subfrags:
        monomers = split_outside(fragment, SequenceConstants.monomer_join, '[]')
        if len(monomers)==1:
            frags[fragment]=subfrags[fragment]
        if len(monomers)==2:
            if monomers[0]==monomers[1]:
                frags[fragment]=subfrags[fragment]
            else:
                fragment_rev = monomers[1]+'-'+monomers[0]
                if fragment not in frags:
                    if fragment_rev not in frags:
                        frags[fragment]=subfrags[fragment]
                    else:
                        frags[fragment_rev]+=subfrags[fragment]
        if len(monomers)==3:
            if monomers[0]==monomers[1] and monomers[1]==monomers[2]:
                frags[fragment]=subfrags[fragment]
            elif monomers[0]==monomers[2]:
                frags[fragment]=subfrags[fragment]
            else:
                fragment_rev = monomers[2]+'-'+monomers[1]+'-'+monomers[0]
                if fragment not in frags:
                    if fragment_rev not in frags:
                        frags[fragment]=subfrags[fragment]
                    else:
                        frags[fragment_rev]+=subfrags[fragment]
    
    # Store the final fragments
    final_frags={}
    for fragment in frags:
        monomers = split_outside(fragment, SequenceConstants.monomer_join, '[]')
        if len(monomers)!=1:
            final_frags[fragment]=int(frags[fragment]/2)
        else:
            final_frags[fragment]=frags[fragment]

    
    # Read the pre-calculated properties of the monomers    
    module_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(module_dir, SequenceConstants.def_path)
    file_path = os.path.join(data_dir, SequenceConstants.def_property)
    with open(file_path, 'r') as handle:
        monomers_prop = pepDescriptors.get_properties(handle) 

    # Define the fingerprint
    fp = [0]*nBits
    
    dict_fp={}
    for biln_frag in final_frags:
        # Iterate over the monomers
        monomers = split_outside(biln_frag, SequenceConstants.monomer_join, '[]')
        code=[]
        for mon in monomers:
            int_code=''
            mon = re.sub(r'\[', '', mon)
            mon = re.sub(r'\]', '', mon)
            for prop in prop_list:
                value=monomers_prop.loc[monomers_prop['name'] == mon, prop].item()
                int_code+=str(int(float(value)))
            code.append(int(int_code))
        
        tag=sum(code)
        # Frequency
        if add_freq:
            tag=tag*final_frags[biln_frag]

        # Generate the hash code and the bit index
        hash_code = hash(tag)
        bit_index = hash_code % nBits  
        
        # Store the results
        dict_fp[biln_frag]=(bit_index, tag, final_frags[biln_frag])
        fp[bit_index]=1
    
    return fp, dict_fp

########################################################################## 
def simMonFP(peptide1, peptide2, radius=2, nBits=1024, add_freq=False, prop_list=['heavy', 'nrot', 'hacc', 'hdon', 'nhet', 'tpsa', 'mw']):
    """
    Function to calculate the similarity between two biln representations

    :param peptide1: BILN representation of peptide 1
    :param peptide2: BILN representation of peptide 2
    :param radius: Number of neighbors around each residue. Maximum 2
    :param nBits: Number of bits in the fingerprint
    :param add_freq: Flag to add the frequency of the fragment in the tag code
    :param prop_list: List of properties used to create the fingerprint bits
    """

    # Create the fingerprints
    fp1, dict_fp1=monomerFP(peptide1, radius=radius, nBits=nBits, add_freq=add_freq, prop_list=prop_list)
    fp2, dict_fp2=monomerFP(peptide2, radius=radius, nBits=nBits, add_freq=add_freq, prop_list=prop_list)

    # Calculate the similarity
    intersection = sum([1 for i, j in zip(fp1, fp2) if i == j and i == 1])
    union = sum([1 for i, j in zip(fp1, fp2) if i == 1 or j == 1])
    similarity = intersection / union
    
    return similarity

########################################################################## 
def only_simMonFP(fp1, fp2):
    """
    Function to calculate the similarity between two already calculated fingerprints

    :param fp1: FP representation of peptide 1
    :param fp2: FP representation of peptide 2
    """

    # Calculate the similarity
    intersection = sum([1 for i, j in zip(fp1, fp2) if i == j and i == 1])
    union = sum([1 for i, j in zip(fp1, fp2) if i == 1 or j == 1])
    similarity = intersection / union
    
    return similarity

##########################################################################
##########################################################################
class pepDescriptors:
    """
    Class to calculate descriptors for peptide containing non-natural amino acids
    """
    ############################################################
    @staticmethod
    def generate_properties(monomer_lib=None, property_lib=None):
        """
        Function to generate the properties file

        :param monomer_lib: Flag to use a monomer library different to the internal one
        :param property_lib: Flag to generate a different property file than the internal
        """
        print("Calculating the properties of the monomers ...")

        # Dictionary with the properties
        dict_df = {'name':[], 'smiles':[], 'mw':[], 'logp':[], 'nrot':[], 'tpsa':[],
                   'hacc':[], 'hdon':[], 'nhet':[], 'heavy':[]}
        
        # Read the SDF file with the monomer information
        if monomer_lib is None:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, SequenceConstants.def_lib_filename)
            with open(file_path, 'rb') as handle:
                df = get_monomer_info(handle)
        else:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, monomer_lib)
            with open(file_path, 'rb') as handle:
                df = get_monomer_info(handle)

        # Extract the smiles from the dataframe
        totalNames=[]
        totalSmiles=[]
        for i, idx in enumerate(df.index):
            name = df.at[idx, 'm_abbr']
            smiles = df.at[idx, 'm_smiles']
            if smiles != 'None':
                mol = Chem.MolFromSmiles(smiles)

                # Calculate the descriptors
                mol_weight = Descriptors.MolWt(mol)
                mol_logp = Descriptors.MolLogP(mol)
                mol_tpsa = Descriptors.TPSA(mol)
                mol_nrot = Descriptors.NumRotatableBonds(mol)

                mol_hacc=Descriptors.NumHAcceptors(mol)
                mol_hdon=Descriptors.NumHDonors(mol)
                mol_nhet=Descriptors.NumHeteroatoms(mol)
                mol_heavy=Descriptors.HeavyAtomCount(mol)

                dict_df['name'].append(name)
                dict_df['smiles'].append(smiles)
                dict_df['mw'].append(mol_weight)
                dict_df['logp'].append(mol_logp)
                dict_df['nrot'].append(mol_nrot)
                dict_df['tpsa'].append(mol_tpsa)
                dict_df['hacc'].append(mol_hacc)
                dict_df['hdon'].append(mol_hdon)
                dict_df['nhet'].append(mol_nhet)
                dict_df['heavy'].append(mol_heavy)

        # Output file in a local data folder
        print("The properties of the monomers have been calculated")
        if property_lib is None:
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, SequenceConstants.def_property)
        else:
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, property_lib)
        with open(file_path, 'w') as handle:
            for i,name in enumerate(dict_df['name']):
                handle.write(f"{name} {dict_df['smiles'][i]} {dict_df['mw'][i]} {dict_df['logp'][i]} {dict_df['nrot'][i]} {dict_df['tpsa'][i]} {dict_df['hacc'][i]} {dict_df['hdon'][i]} {dict_df['nhet'][i]} {dict_df['heavy'][i]}\n")

    ############################################################################
    @staticmethod
    def get_properties(path):
        """
        Read the property file to calculate the descriptors
        :param path: text stream with the property

        :return: property dictionary
        """
        # Dictionary with the properties
        dict_df = {'name':[], 'smiles':[], 'mw':[], 'logp':[], 'nrot':[], 'tpsa':[], 'hacc':[], 'hdon':[], 'nhet':[], 'heavy':[]}

        data = path.read().split('\n')
        for line in data:
            if line:
                fields = line.split()
                dict_df['name'].append(fields[0])
                dict_df['smiles'].append(fields[1])
                dict_df['mw'].append(fields[2])
                dict_df['logp'].append(fields[3])
                dict_df['nrot'].append(fields[4])
                dict_df['tpsa'].append(fields[5])
                dict_df['hacc'].append(fields[6])
                dict_df['hdon'].append(fields[7])
                dict_df['nhet'].append(fields[8])
                dict_df['heavy'].append(fields[9])

        df = pd.DataFrame(dict_df)

        return df
    
    ##########################################################################
    @staticmethod
    def get_alternate_path(biln: str):
        """
        Function to return alternate paths of the BILN in case of branching
        
        :param biln: BILN representation of the peptide
        :type biln: str
        """
        chains = biln.split(SequenceConstants.chain_separator)

        newBilns = []
        mainBilns = {}
        positions = []
        for mNum, m in enumerate(chains):
            monomers = split_outside(m, SequenceConstants.monomer_join, '[]')
            mainBilns[mNum + 1] = []
            for pos,res in enumerate(monomers):
                if SequenceConstants.linker_symbol in res:
                    positions.append((mNum+1,pos+1,int(res.split(SequenceConstants.linker_symbol)[1][0])))
                mon = re.sub("\(\d+,\d+\)", "", res)
                mainBilns[mNum + 1].append(mon)
            
        connectors = {}
        for pair in positions:
            chainInt = pair[0]
            numInt = pair[1]
            posLink = pair[2]
            if posLink not in connectors:
                connectors[posLink] = [(chainInt, numInt)]
            else:
                connectors[posLink].append((chainInt, numInt))

        if connectors:
            for con in connectors:
                p1 = connectors[con][0]
                p2 = connectors[con][1]
                
                temporal_chain = mainBilns[p1[0]][:p1[1]]
                new_chain = temporal_chain + list(reversed(mainBilns[p2[0]][:p2[1]]))
                new_biln = SequenceConstants.monomer_join.join(new_chain)
                newBilns.append(new_biln)

                temporal_chain = mainBilns[p1[0]][p1[1]:]
                new_chain = mainBilns[p2[0]]+temporal_chain
                new_biln = SequenceConstants.monomer_join.join(new_chain)
                newBilns.append(new_biln)

        return newBilns, chains

    ############################################################
    def __init__(self, sequence, property_lib=None):
        """
        Start class by assigning general values

        :param sequence: Sequence of the peptide in BILN format
        :param generate_properties:Generate or not a pickle file with the pre-calculated properties. By default the file is provided
        """

        # Check if the peptide is in BILN format.
        if "-" not in sequence:
            sequence='-'.join(list(sequence))
            
        self.sequence = sequence
        
        if property_lib is None:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, SequenceConstants.def_property)
        else:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            data_dir = os.path.join(module_dir, SequenceConstants.def_path)
            file_path = os.path.join(data_dir, property_lib)
            
        with open(file_path, 'r') as handle:
            self.monomers_prop = pepDescriptors.get_properties(handle) 

    ##########################################################################
    def moranCorrelation(self, nlag=5, label=False, AAidxName = ['nrot', 'logp', 'tpsa', 'mw']) -> dict:
        """
        Calculate amino-acid based moran correlation descriptors

        :param nlag: number of neighbors accounted for the descriptors. 
                     The values should be smaller than the full peptide length
        :param label: label to differentiate the chains
        :param AAidxName: list with the properties taken into account for the descriptors

        :return descriptors: a dictionary with the numerical descriptors
        """
        fullMonomers = split_outside(self.sequence, SequenceConstants.monomer_join, '[]')

        monomers = []
        for mon in fullMonomers:
            mon = re.sub(r'\[', '', mon)
            mon = re.sub(r'\]', '', mon)
            mon = re.sub("\(\d+,\d+\)", "", mon)
            monomers.append(mon)

        totalMon = len(monomers)

        if totalMon == 1:
            print('Warning: The sequence only have one amino acid and no Moran correlation descriptors will be calculated')

        if totalMon < nlag + 1:
            print('Warning: the sequence should be larger than nlag+1: ' + str(nlag + 1) + '. The nlag was refactored based on the chain length')
            nlag = totalMon -1

        AAidx = []
        for j in AAidxName:
            property_values = []
            for i, res in enumerate(fullMonomers):
                mon = re.sub(r'\[', '', res)
                mon = re.sub(r'\]', '', mon)
                mon = re.sub("\(\d+,\d+\)", "", mon)
                property_values.append(self.monomers_prop.loc[self.monomers_prop['name'] == mon, j].item())
            AAidx.append(property_values)

        AAidx1 = np.array([float(j) for i in AAidx for j in i])
        AAidx = AAidx1.reshape((len(AAidx), len(fullMonomers)))

        propMean = np.mean(AAidx, axis=1)
        propStd = np.std(AAidx, axis=1)

        for i in range(len(AAidx)):
            for j in range(len(AAidx[i])):
                AAidx[i][j] = (AAidx[i][j] - propMean[i]) / propStd[i]

        index = {}
        for i in range(len(fullMonomers)):
            mon = re.sub(r'\[', '', fullMonomers[i])
            mon = re.sub(r'\]', '', mon)
            mon = re.sub("\(\d+,\d+\)", "", mon)
            index[mon] = i

        # Calculate the descriptors and include them in the dictionary
        descriptors = {}
        for prop in range(len(AAidxName)):
            xmean = sum([AAidx[prop][index[aa]] for aa in monomers]) / totalMon
            for n in range(1, nlag + 1):
                if totalMon > nlag:
                    num = sum([(AAidx[prop][index.get(monomers[j], 0)] - xmean) *
                                    (AAidx[prop][index.get(monomers[j + n], 0)] - xmean)
                                    for j in range(totalMon - n)]) / (totalMon - n)
                    den = sum(
                        [(AAidx[prop][index.get(monomers[j], 0)] - xmean) ** 2 for j in range(totalMon)]) / totalMon
                    rn = num / den
                else:
                    rn = 'NA'
                if label:
                    descriptors[AAidxName[prop] + '-lag' + str(n) + '-' + label] = rn
                else:
                    descriptors[AAidxName[prop] + '-lag' + str(n)] = rn

        return descriptors

    # End of pepDescriptors class
    ############################################################

##########################################################################
##########################################################################
class graphPep:

    """Class to generate a graph version of the peptide"""

    def __init__(self, peptide):
        """
        Class to generate a graph-based representation of the peptide

        :param peptide: peptide sequence in biln format
        """

        # Main variables
        self.peptide = peptide

        # Internal variables
        self.aa_list = {}
        self.linkers = {}
        self.pep_interactions = []
        self.names = []

        # Get chains
        self.chains = []
        self.pairs = {}

        if SequenceConstants.chain_separator in self.peptide:
            chains=self.peptide.split(SequenceConstants.chain_separator)
            counter_aa={}
            for i,chain in enumerate(chains):
                counter_aa[i]=0
                monomers = split_outside(chain, SequenceConstants.monomer_join, '[]')
                for mon in monomers:
                    if mon in SequenceConstants.aminoacids:
                        counter_aa[i]+=1
            chain_order = sorted(counter_aa.items(), key=lambda x: x[1], reverse=True)
            for pair in chain_order:
                self.chains.append(chains[pair[0]])
        else:
            self.chains.append(peptide)
        
        for i,chain in enumerate(self.chains):
            monomers = split_outside(chain, SequenceConstants.monomer_join, '[]')
            for j,mon in enumerate(monomers):
                
                pat=re.findall("\(\d+,\d+\)",mon)
                if len(pat)>=1:
                    num_con=list(mon).count('(')
                    
                    mon_mod = re.sub(r'\(', '+', mon)
                    mon_orig = re.sub("\(\d+,\d+\)", "", mon)
                    for x in range(0,num_con):
                        frac=mon_mod.split('+')[x+1]
                        link_id = frac.split(',')[0]
                        if link_id not in self.pairs:
                            self.pairs[link_id]=[(mon_orig+'('+frac,j+1,i)]
                        else:
                            self.pairs[link_id].append((mon_orig+'('+frac,j+1,i))

    ########################################################################################
    def __getLinkers(self):
        """
        Extract all posible inter and intra connections in the peptide
        """

        # Capture intra and inter chains
        if self.pairs:
            for pair in self.pairs:
                chTemp = []
                for m in self.pairs[pair]:
                    chainInt = m[2]+1
                    chTemp.append(chainInt)
                    numInt = m[1]
                    monInt = m[0]
                    posLink = int(pair)
                    if posLink not in self.linkers:
                        self.linkers[posLink] = [(chainInt, numInt, monInt)]
                    else:
                        self.linkers[posLink].append((chainInt, numInt, monInt))
                    if len(chTemp)==2:
                        if chTemp[0]==chTemp[1]:
                            self.linkers[posLink].append('INTRA')
                        else:
                            self.linkers[posLink].append('INTER')

    ########################################################################################
    def __assignInfoPerChain(self):
        """
        Read the peptide by chains and annotate information about the nodes

        :return: More lists with network parameters
        """

        # Iterate over the chains
        for j,ch in enumerate(self.chains):
            self.aa_list[j+1]=[]
            tempSeq = ch
            monomers = split_outside(tempSeq, SequenceConstants.monomer_join, '[]')

            # Iterate over monomers
            for i,res in enumerate(monomers):
                mon1 = re.sub("\(\d+,\d+\)", "", res)
                mon2 = re.sub(r"[\[]", "", mon1)
                mon = re.sub(r"[\]]", "", mon2)

                self.aa_list[j+1].append(mon)

    ########################################################################################
    def __graphProperties(self):
        """
        Add properties related with node label size and positioning, as well as shape. The regular interactions are also included

        :return: Lists with the assigned attributes
        """

        # Get graph order and properties of the nodes and edges
        for j,ch in enumerate(self.chains):
            names_temp=[0]*len(self.aa_list[j+1])
            for i,amino in enumerate(self.aa_list[j+1]):

                name_part=''
                pos_part = ''
                chain_part = ''
                name_part =amino
                names_temp[i] = name_part

                if i != len(self.aa_list[j+1])-1:
                    self.pep_interactions.append((len(self.names)+i,len(self.names)+i+1))

            self.names += names_temp

    ########################################################################################
    def __addInfoLinkers(self):
        """
        Internal function to add information about the intra and inter connections.
        """

        if self.linkers:
            
            for bridge in self.linkers:
                nodes = []
                chainNodes = []
                for count in range(0,2):
                    chainLink = self.linkers[bridge][count][0]
                    posLink = self.linkers[bridge][count][1]

                    if chainLink == 1:
                        node = posLink - 1
                    elif chainLink == 2:
                        node = len(self.aa_list[chainLink-1]) + posLink - 1
                    elif chainLink == 3:
                        node = len(self.aa_list[chainLink-1]) + len(self.aa_list[chainLink-2]) + posLink - 1
                    elif chainLink == 4:
                        node = len(self.aa_list[chainLink-1]) + len(self.aa_list[chainLink-2]) + len(self.aa_list[chainLink-3]) + posLink - 1
                    else:
                        for i in range(0, chainLink-1):
                            pos=i+1
                            print(chainLink-pos)
                            if i ==0:
                                node=len(self.aa_list[chainLink-pos])
                            else:
                                node+=len(self.aa_list[chainLink-pos])
                        node+=posLink-1
                        
                    nodes.append(node)
                    chainNodes.append(chainLink)

    ########################################################################################
    def createGraph(self):
        """
        Creates the iGraph object with the initial set of parameters
        """

        # Call internal functions
        self.__getLinkers()
        self.__assignInfoPerChain()
        self.__graphProperties()
        self.__addInfoLinkers()

        # Generate ipgraph object
        self.g = Graph(self.pep_interactions)
        self.g.vs["name"] = self.names
        
############################################################
## End of graphPep class
############################################################

############################################################
## End of similarity.py
############################################################