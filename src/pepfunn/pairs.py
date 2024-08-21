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

# BioPython
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

# System
import pandas as pd
import os
import itertools
import warnings
import sys
import matplotlib.pyplot as plt
from statistics import mean
import ast

# Internal
from pepfunn.sequence import SequenceConstants
from pepfunn.similarity import Alignment

########################################################################################
# Classes and functions
########################################################################################

class MatchedPairs:

    '''Class with functions to calculate matched pairs based on references and between each pair'''

    ########################################################################################
    @staticmethod
    def get_mutations(peptide1, peptide2):
        """
        Function to extract the mutations and distance for two sequences having the same size

        :param peptide1: First sequence
        :param peptide2: Second sequence

        :return mutation string and the distance (close to zero more similar)
        """
        # Generate the list of mutations
        mutations = []
        matches = 0    
        for pos, p in enumerate(peptide1):
            
            # Generate the tuples with the pair of amino acids to compare
            if p != peptide2[pos]:
                mutations.append(f"({pos+1}){p}/{peptide2[pos]}")
            else:
                matches+=1

        mutations_str = " | ".join(mutations)
        distance = len(peptide1)-matches

        return mutations_str, distance
    
    ########################################################################################
    @staticmethod
    def get_mutations_simple(peptide1, peptide2, ref1, ref2):
        """
        Function to extract the mutations and distance for two sequences of different size and with reference peptides

        :param peptide1: First sequence
        :param peptide2: Second sequence
        :param ref1: Reference sequence for first peptide
        :param ref2: Reference sequence for second peptide

        :return mutation string and the distance (close to zero more similar)
        """
        # Generate the list of mutations
        mutations = []
        matches = 0

        positions_ref1=[]
        for i,ele in enumerate(ref1):
            if ele != '-':
                positions_ref1.append(i+1)
            else:
                positions_ref1.append(0)
        
        positions_ref2=[]
        for i,ele in enumerate(ref2):
            if ele != '-':
                positions_ref2.append(i+1)
            else:
                positions_ref2.append(0)

        for pos, p in enumerate(peptide1):
            
            # Generate the tuples with the pair of amino acids to compare
            if p != peptide2[pos]:
                if ref1[pos] != '-' and ref2[pos] != '-':
                    mutations.append(f"({positions_ref1[pos]}){p}/{peptide2[pos]}")
                else:
                    matches+=1
            else:
                matches+=1
        
        mutations_str = " | ".join(mutations)
        distance = len(positions_ref1)-matches

        return mutations_str, distance
    
    ########################################################################################
    @staticmethod
    def diff_property(value_a, value_b, operation='divide'):
        """
        Method to extract the difference between two metric values

        :param value_A: First value to compare
        :param value_B: Second value to compare
        :param operation: Operation to compare the values (substract or divide)
        """

        if operation not in ['substract', 'divide']:
            raise ValueError("The operation should be substract or divide. Please correct")
            sys.exit(1)

        if operation == 'substract':
            # Subtract the smaller from the bigger number
            diff_prop = value_a - value_b
        elif operation == 'divide':
            # Divide the bigger by the smaller number, avoiding division by zero
            if value_b != 0:
                diff_prop = value_a / value_b
            else:
                # This handles the case where one or both values are zero
                diff_prop = None  

        min_prop = min(value_a, value_b)
        max_prop = max(value_a, value_b)

        return diff_prop, min_prop, max_prop
    
    ########################################################################################
    @staticmethod
    def alignment_pairs(peptide1, peptide2, matrix="BLOSUM62"):
        """
        Function to run a simple pairwise alignment and extract the aligned sequences and BLOSUM score

        :param peptide1: Sequence of first peptide to compare
        :param peptide2: Sequence of second peptide to compare
        :param matrix: Name of the matrix used to compare. Based on what is available in Biopython
        """
                
        # Create a PairwiseAligner object
        aligner = PairwiseAligner()
        sub_mat = substitution_matrices.load(matrix)
        aligner.substitution_matrix = sub_mat

        # Penalty of gaps
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.1

         # Perform the alignment using the new aligner
        alignments = aligner.align(peptide1, peptide2)

        # Take the first (best) alignment
        alignment = alignments[0]

        # Calculate the BLOSUM score
        blosum_score = alignment.score

        # Prepare the alignment string
        formatted_alignment = alignment.format().splitlines()

        # Extract only the aligned sequences with gaps, removing the leading labels
        seq1_aligned_with_gaps = formatted_alignment[0].split()[-1]
        seq2_aligned_with_gaps = formatted_alignment[2].split()[-1]

        return seq1_aligned_with_gaps, seq2_aligned_with_gaps, blosum_score

    ########################################################################################
    @staticmethod
    def get_sequences(df, id_column, seq_column, seq_ref, name_ref='', threshold=0.4):
        """
        Extract sequences from a dataframe containing the ID and sequence information

        :param df: Input dataframe with all the required data
        :param id_column: Name of the column containing the molecules ID
        :param seq_column: Name of the column containing the molecules sequences
        :param seq_ref: Reference peptide to compare. The sequences should be equal or larger than the reference
        :param name_ref: Name of the reference peptide (optional). If not provided, the REF code is assigned
        :param threshold: Matched threshold to assign the reference or not to each peptide

        :return dataframe containing the mapped sequences to the reference, and the remaining sections
        """
        if name_ref:
            sequences={'ID':[], 'Begin':[], 'End':[], 'Peptide':[], name_ref:[]}
        else:
            sequences={'ID':[], 'Begin':[], 'End':[], 'Peptide':[], 'REF':[]}
            name_ref='REF'

        # Iterate over the input dataframe
        for index, row in df.iterrows():
            id_val=row[id_column]
            if id_val not in sequences['ID']:
                try:
                    peptide = row[seq_column]
                    fragment=''
                    if len(peptide) >= len(seq_ref):
                        frags={}
                        for i in range(0, len(peptide) - len(seq_ref) + 1):
                            frag = peptide[i:i + len(seq_ref)]
                            m=Alignment.align_samelen_local(frag,seq_ref)
                            if m > len(seq_ref) * threshold:
                                frags[frag]=m
                                
                        sorted_frags = sorted(frags.items(), key=lambda item: item[1])
                        if sorted_frags:
                            fragment=sorted_frags[0][0]
                    else:
                        warnings.warn(f"The sequence is smaller than the substrate. No fragment has been obtained for {id_val}")    

                    index = peptide.find(fragment)
                    begin=peptide[0:index]
                    end=peptide[index+len(list(fragment)):]    
                
                    sequences['ID'].append(id_val)
                    sequences[name_ref].append(fragment)
                    sequences['Begin'].append(begin)
                    sequences['End'].append(end)
                    sequences['Peptide'].append(peptide)
                except:
                    pass
        
        # Generate dataframe
        df_sequences = pd.DataFrame(sequences)

        return df_sequences
    
    ########################################################################################
    @staticmethod
    def get_pairs(df, df_sequences, property_columns, name_ref='', operation_columns=''):
        """
        Function to extract all the pairs based on sequence distance and difference in their properties

        :param df: Input dataframe with the sequences and the properties
        :param df_sequence: Dataframe obtained with the mapped reference peptide on their sequences
        :param property_columns: Columns with the properties to compare in the pairs
        :param name_ref: Name of the reference peptide (optional)
        :param operation_columns: Operations used to compare the numerical values of the dataframes. 
                                  Substract and Divide are accepted
        
        :return dataframe with all the pairs calculated containing the mutations and the distances
        """

        if len(df.index) != len(df_sequences.index):
            raise ValueError("The number of sequences in both dataframes should be the same. Please correct")
            sys.exit(1)

        if not name_ref:
            name_ref='REF'

        results={'ID_A':[], 'ID_B':[], 'Mutations': [], 'Distance': []}

        # generate all possible combinations of row indices
        index_pairs = list(itertools.combinations(df.index, 2))
        # filter out pairs where both elements are the same
        index_pairs = [(a, b) for (a, b) in index_pairs if a != b]
        # add the reverse pairs
        reverse_pairs = [(b, a) for (a, b) in index_pairs]
        # Final combinations
        index_pairs = index_pairs + reverse_pairs

        for i,pair in enumerate(index_pairs):

            # Extract information for peptide1
            id1=df_sequences.loc[pair[0], 'ID']
            frag_1= df_sequences.loc[pair[0], name_ref]
            
            id2=df_sequences.loc[pair[1], 'ID']
            frag_2= df_sequences.loc[pair[1], name_ref]

            # Extract the mutations, distances and diff in properties
            flag_prop=0
            if len(frag_1) == len(frag_2):
                flag_prop=1
                
                results['ID_A'].append(id1)
                results['ID_B'].append(id2)

                if frag_1 and frag_2:
                    mutations, distance=MatchedPairs.get_mutations(frag_1, frag_2)
                    results['Mutations'].append(mutations)    
                    results['Distance'].append(distance)
                    
            # Calculate the property differences, min and max for each property
            if flag_prop==1:
                if property_columns:
                    if not operation_columns:
                        operation_columns=['divide']*len(property_columns)
                    else:
                        if len(operation_columns)!=len(property_columns):
                            raise ValueError("The length of the operation_columns should be the same of the property_columns. Exiting ...")
                            sys.exit(1)
                        for oper in operation_columns:
                            if oper not in ['substract','divide']:
                                raise ValueError("The operation value should be substract or divide. Please correct. Exiting ...")
                                sys.exit(1)

                    for j,prop in enumerate(property_columns):
                        value_a = float(df.loc[pair[0], prop])
                        value_b = float(df.loc[pair[1], prop])

                        diff_prop, min_prop, max_prop = MatchedPairs.diff_property(value_a, value_b, operation=operation_columns[j])
                        diff_label = f'Diff_{prop}'
                        min_label = f'Min_{prop}'
                        max_label = f'Max_{prop}'
                        if diff_label not in results:
                            results[diff_label]=[diff_prop]
                        else:
                            results[diff_label].append(diff_prop)

                        if min_label not in results:
                            results[min_label]=[min_prop]
                        else:
                            results[min_label].append(min_prop)

                        if max_label not in results:
                            results[max_label]=[max_prop]
                        else:
                            results[max_label].append(max_prop)

        df_final = pd.DataFrame(results)
        
        return df_final
    
    ########################################################################################
    @staticmethod
    def get_sequences_simple(df, id_column, seq_column, seq_ref=''):
        """
        Get sequences with or without a reference sequence. The generation of gaps is allowed

        :param df: dataframe with the source data
        :param id_column: ID number
        :param seq_column: Column containing the peptide sequences
        :param seq_ref: Reference peptide to compare (optional)
        
        :return: dataframe with the sequences mapped or not to a reference sequence
        """

        if not seq_ref:
            sequences={'ID':[], 'Peptide':[]}
        else:
            sequences={'ID':[], 'Peptide':[], 'Aligned_Ref':[], 'Aligned_Pep':[]}
        
        for index, row in df.iterrows():
            id_val=row[id_column]
            if id_val not in sequences['ID']:
                try:
                    peptide = row[seq_column]
                    sequences['ID'].append(id_val)
                    sequences['Peptide'].append(peptide)

                    # Case if a reference sequence is provided
                    if seq_ref:
                        seqref_aligned, seqnew_aligned, score=MatchedPairs.alignment_pairs(seq_ref, peptide)
                        sequences['Aligned_Ref'].append(seqref_aligned)
                        sequences['Aligned_Pep'].append(seqnew_aligned)
                except:
                    pass
        
        # Generate dataframe
        df_sequences = pd.DataFrame(sequences)

        return df_sequences

    ########################################################################################
    @staticmethod
    def get_pairs_simple(df, df_sequences, property_columns, seq_ref='', operation_columns=''):
        """"
        Function to extract all the pairs based on sequence distance and difference in their properties.
        Adapted for the case where a reference is provided or not, as well as the generation of gaps.

        :param df: Dataframe with the dataset
        :param df_sequences: Dataframe with the sequences mapped to the reference peptides
        :param property_columns: List with the columns used to compare the properties
        :param seq_ref: Reference peptide to compare (optional)
        :param operation_columns: Operations used to compare the numerical values of the dataframes. 
                                  Substract and Divide are accepted
        
        :return dataframe with all the pairs calculated containing the mutations and the distances
        """

        if len(df.index) != len(df_sequences.index):
            raise ValueError("The number of sequences in both dataframes should be the same. Please correct")
            sys.exit(1)

        if seq_ref:
            results={'ID_A':[], 'ID_B':[], 'Mutations':[], 'Distance':[]}
        else:
            results={'ID_A':[], 'ID_B':[], 'Mutations':[]}
        
        # generate all possible combinations of row indices
        index_pairs = list(itertools.combinations(df.index, 2))
        # filter out pairs where both elements are the same
        index_pairs = [(a, b) for (a, b) in index_pairs if a != b]
        # add the reverse pairs
        reverse_pairs = [(b, a) for (a, b) in index_pairs]
        # Final combinations
        index_pairs = index_pairs + reverse_pairs
        
        for i,pair in enumerate(index_pairs):
            
            # Extract information for peptide1
            id1=df_sequences.loc[pair[0], 'ID']
            if seq_ref:
                frag1= df_sequences.loc[pair[0], 'Aligned_Pep']
                ref1=df_sequences.loc[pair[0], 'Aligned_Ref']
            else:
                frag1= df_sequences.loc[pair[0], 'Peptide']

            id2=df_sequences.loc[pair[1], 'ID']
            if seq_ref:
                frag2= df_sequences.loc[pair[1], 'Aligned_Pep']
                ref2= df_sequences.loc[pair[1], 'Aligned_Ref']
            else:
                frag2= df_sequences.loc[pair[1], 'Peptide']

            # Extract the mutations, distances and diff in properties
            if seq_ref:
                flag_prop=0
                if len(frag1) == len(frag2):

                        flag_prop=1

                        results['ID_A'].append(id1)
                        results['ID_B'].append(id2)
                        
                        if frag1 and frag2:
                            mutations, distance = MatchedPairs.get_mutations_simple(frag1, frag2, ref1, ref2)
                            results[f'Mutations'].append(mutations)
                            results[f'Distance'].append(distance)
                        else:
                            results[f'Mutations'].append('')
                            results[f'Distance'].append(99)
            else:
                flag_prop=1

                results['ID_A'].append(id1)
                results['ID_B'].append(id2)
                seq1_aligned, seq2_aligned, score=MatchedPairs.alignment_pairs(frag1, frag2)
                positions_pep1=[]
                for i,ele in enumerate(seq1_aligned):
                    if ele != '-':
                        positions_pep1.append(i+1)
                    else:
                        positions_pep1.append(0)
                
                positions_pep2=[]
                for i,ele in enumerate(seq2_aligned):
                    if ele != '-':
                        positions_pep2.append(i+1)
                    else:
                        positions_pep2.append(0)

                mutations=[]
                for pos, p in enumerate(seq1_aligned):
                    # Generate the tuples with the pair of amino acids to compare
                    real_pos1 = positions_pep1[pos]
                    real_pos2 = positions_pep2[pos]
                    
                    val1=p
                    val2=seq2_aligned[pos]

                    if val1 != val2:
                        mutations.append(f"({real_pos1}){val1}/({real_pos2}){val2}")
                
                # Join the mutations
                mutations_str = " | ".join(mutations)
                        
                results[f'Mutations'].append(mutations_str)
                
                # Calculate the property differences, min and max for each property
            if flag_prop==1:
                if property_columns:
                    if not operation_columns:
                        operation_columns=['divide']*len(property_columns)
                    else:
                        if len(operation_columns)!=len(property_columns):
                            raise ValueError("The length of the operation_columns should be the same of the property_columns. Exiting ...")
                            sys.exit(1)
                        for oper in operation_columns:
                            if oper not in ['substract','divide']:
                                raise ValueError("The operation value should be substract or divide. Please correct. Exiting ...")
                                sys.exit(1)

                    for j,prop in enumerate(property_columns):
                        value_a = float(df.loc[pair[0], prop])
                        value_b = float(df.loc[pair[1], prop])
                        
                        diff_prop, min_prop, max_prop = MatchedPairs.diff_property(value_a, value_b, operation=operation_columns[j])
                        diff_label = f'Diff_{prop}'
                        min_label = f'Min_{prop}'
                        max_label = f'Max_{prop}'
                        if diff_label not in results:
                            results[diff_label]=[diff_prop]
                        else:
                            results[diff_label].append(diff_prop)
                        
                        if min_label not in results:
                            results[min_label]=[min_prop]
                        else:
                            results[min_label].append(min_prop)
                        
                        if max_label not in results:
                            results[max_label]=[max_prop]
                        else:
                            results[max_label].append(max_prop)

        df_final = pd.DataFrame(results)

        return df_final
    
    ########################################################################################
    @staticmethod
    def plot_activities(filtered_df, property_columns, thresholds, color='red', out_name='activities.png'):
        """"
        Method to plot the activities for more than two properties

        :param filtered_df: dataframe used to extract the properties
        :param property_columns: list with the columns used to filter
        :param thresholds: list with theshold values for the properties
        :param color: color used to plot the activities
        :param out_name: name of the image file with the plot

        :return png file with the generated plot
        """

        if len(property_columns) < 2:
            raise ValueError(f"At least two properties should be included to generate the plot")
            sys.exit(1)
        
        # df with the properties
        plot_df=filtered_df[property_columns]

        if color=='red':
            cm = plt.cm.get_cmap('Reds')
        elif color=='blue':
            cm = plt.cm.get_cmap('Blues')
        elif color=='green':
            cm = plt.cm.get_cmap('Greens')
        else:
            raise ValueError(f"Only red, blue, and green are supported. Red will be used by default")
            cm = plt.cm.get_cmap('Reds')

        values=[]
        conditions=[]
        # Scatter plot
        for i,prop in enumerate(property_columns):
            values.append(plot_df[prop].values.tolist())
            conditions.append(plot_df[prop] > thresholds[i])
        
        # filter the DataFrame based on the conditions and export column 'A' to a list.
        if len(conditions)==2:
            filtered_values_1 = plot_df.loc[conditions[0] & conditions[1], property_columns[0]].tolist()
            filtered_values_2 = plot_df.loc[conditions[0] & conditions[1], property_columns[1]].tolist()

            plt.scatter(values[0],values[1], marker="o",s=25,alpha=1.0, cmap=cm, edgecolors='black', linewidths=0.5)
            plt.plot(filtered_values_1, filtered_values_2, marker="o", markersize=5,alpha=1.0,markeredgewidth=1.0, markeredgecolor='green',mfc='none',linestyle="None")
            
        if len(conditions)==3:
            filtered_values_1 = plot_df.loc[conditions[0] & conditions[1] & conditions[2], property_columns[0]].tolist()
            filtered_values_2 = plot_df.loc[conditions[0] & conditions[1] & conditions[2], property_columns[1]].tolist()
            
            plt.scatter(values[0],values[1],c=values[2], marker="o",s=25,alpha=1.0,cmap=cm, edgecolors='black', linewidths=0.5)
            plt.plot(filtered_values_1, filtered_values_2, marker="o", markersize=5,alpha=1.0,markeredgewidth=1.0, markeredgecolor='green',mfc='none',linestyle="None")

        
        plt.axvline(x=thresholds[0], linestyle='dashed', color='black', linewidth=0.5)
        plt.axhline(y=thresholds[1], linestyle='dashed', color='black', linewidth=0.5)
        axes = plt.gca()
        plt.title('Plot properties')
        plt.xlabel(property_columns[0])
        plt.ylabel(property_columns[1])
        cbar = plt.colorbar()
        cbar.ax.invert_yaxis() 
        plt.savefig(out_name)
        plt.show()

############################################################
## End of pairs.py
############################################################
