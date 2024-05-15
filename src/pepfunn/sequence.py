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
import re
import pickle
import math
import os
import warnings
import sys
import tempfile

# BioPython
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import rdMolDescriptors

########################################################################################
# Classes and Functions
########################################################################################

class SequenceConstants:
    """
    A class to hold defaults values
    """
    def_path = "data"
    def_lib_filename = "monomers.sdf"
    def_matrix = "matrix.txt"
    def_property = "property.txt"
    monomer_join = "-"
    chain_separator = "."
    csv_separator = ","
    aminoacids ={"A" :"ALA", "D" :"ASP", "E" :"GLU", "F" :"PHE", "H" :"HIS",
                 "I" :"ILE", "K" :"LYS", "L" :"LEU", "M" :"MET", "G" :"GLY",
                 "N" :"ASN", "P" :"PRO", "Q" :"GLN", "R" :"ARG", "S" :"SER",
                 "T" :"THR", "V" :"VAL", "W" :"TRP", "Y" :"TYR", "C" :"CYS"}

# End of Sequence class-related constants definition.
############################################################

class Sequence:
    """
    Class with functions to perform different type of analysis using a peptide sequence as an object
    """

    def __init__(self, sequence, format='fasta', report_liabilities=False):
        """
        Inititalize the class calculating some basic properties

        :param sequence: Peptide sequence
        :param format: format used to represent the sequence. Main is fasta, but  biln can be used
        :param report_liabilities: Flag to generate a file with the report of the sequence liabilities

        :return Based on the sequence, the class will generate all the possible properties based on the internal functions
        """

        if format=='fasta':
            try:
                self.sequence = sequence
                helm = ".".join(list(self.sequence))
                self.mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
                self.smiles=Chem.MolToSmiles(self.mol)
            except ValueError:
                warnings.warn(f"Please check the input sequence is in FASTA format")
                sys.exit(1)

        if format=='biln':
            self.sequence = get_peptide(sequence)
            helm = ".".join(list(self.sequence))
            self.mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" % helm)
            self.smiles=Chem.MolToSmiles(self.mol)

        self.pH = 7
        self.solubility_rules_failed = 0
        self.set_sol_rules = []
        self.length_peptide = len(self.sequence)
        self.synthesis_rules_failed = 0
        self.set_syn_rules = []

        # Calculate internal functions
        self.__compute_peptide_charges() # Add pH function
        self.__calculate_properties_from_mol()
        self.__calculate_properties_from_sequence()
        self.__solubility_rules()
        self.__synthesis_rules()

        if report_liabilities:
            # Rules
            sol_rules={1:"Discard if the number of charged and/or of hydrophobic amino acids exceeds 45%",
                    2:"Discard if the absolute total peptide charge at pH 7 is more than +1",
                    3:"Discard if the number of glycine or proline is more than one in the sequence",
                    4:"Discard if the first or the last amino acid is charged",
                    5:"Discard if any amino acid represents more than 25% of the total sequence"}

            syn_rules={6:"Discard if 2 prolines are consecutive",
                    7:"Discard if the motifs DG and DP are present in the sequence",
                    8:"Discard if the sequences ends with N or Q residues",
                    9:"Discard if there are charged residues every 5 amino acids",
                    10:"Discard if there are oxidation-sensitive amino acids (M, C or W)"}

            # Check the specific solubility rules violated
            rules_sol_violated=""
            for i,rules in enumerate(self.set_sol_rules):
                if rules==1: 
                    rules_sol_violated+=str(i+1)+","
            if rules_sol_violated: 
                sol_violated=rules_sol_violated[:-1]
            else: 
                sol_violated="None"

            # Check the specific synthesis rules violated
            rules_syn_violated=""
            for i,rules in enumerate(self.set_syn_rules):
                if rules==1: 
                    rules_syn_violated+=str(i+6)+","
            if rules_syn_violated: 
                syn_violated=rules_syn_violated[:-1]
            else: 
                syn_violated="None"

            # Report
            sequence_report=open("rules_report.txt","w")
            sequence_report.write("List of failed rules identified by number. The rule ids are explained at the end of the report.\n")
            sequence_report.write("###########################################\n")
            sequence_report.write("A total of {} solubility rules failed. The rules id(s) are: {}.\n".format(self.solubility_rules_failed,sol_violated))
            sequence_report.write("A total of {} synthesis rules failed. The rules id(s) are: {}.\n".format(self.synthesis_rules_failed,syn_violated))
            sequence_report.write("###########################################\n")
            sequence_report.write("\nThe higher the number of rules violated, the lower the probability to be solubilized or synthesized experimentally (https://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/).\n")
            sequence_report.write("\n- List of solubility rules violations:\n")
            for key in sol_rules:
                sequence_report.write("{}. {}\n".format(key,sol_rules[key]))
            sequence_report.write("\n- List of synthesis rules violations:\n")
            for key in syn_rules:
                sequence_report.write("{}. {}\n".format(key,syn_rules[key]))
            sequence_report.close()
        
    ############################################################################
    def __compute_peptide_charges(self, pH_internal=7):
        """
        Function to calculate the average net charge based on pka values

        :param pH_internal -- By default is 7

        :return net_charge: The net charge based on reported pka values
        """
        # Set the general variables and pka terms
        self.pH = pH_internal
        self.netCharge = 0.0
        pka_alpha_amino = {'G': 9.60, 'A': 9.69, 'V': 9.62, 'L': 9.60, 'I': 9.68, 'M': 9.21, 'F': 9.13, 'W': 9.39,
                           'P': 10.60, 'S': 9.15,
                           'T': 9.10, 'C': 10.78, 'Y': 9.11, 'N': 8.84, 'Q': 9.13, 'D': 9.82, 'E': 9.67, 'K': 8.95,
                           'R': 9.04, 'H': 9.17}
        pka_alpha_carboxy = {'G': 2.34, 'A': 2.34, 'V': 2.32, 'L': 2.36, 'I': 2.36, 'M': 2.28, 'F': 1.83, 'W': 2.38,
                             'P': 1.99, 'S': 2.21,
                             'T': 2.63, 'C': 1.71, 'Y': 2.2, 'N': 2.02, 'Q': 2.17, 'D': 2.09, 'E': 2.19, 'K': 2.18,
                             'R': 2.17, 'H': 1.82}
        pka_sidechain_positive = {'K': 10.79, 'R': 12.48, 'H': 6.04}
        pka_sidechain_negative = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10.07}

        # Calculate the net charge for the extreme groups (without modifications)
        amino = self.sequence[0]
        carboxy = self.sequence[-1]
        self.netCharge = self.netCharge + (math.pow(10, pka_alpha_amino[amino]) / (
                    math.pow(10, pka_alpha_amino[amino]) + math.pow(10, self.pH)))
        self.netCharge = self.netCharge - (
                    math.pow(10, self.pH) / (math.pow(10, pka_alpha_carboxy[carboxy]) + math.pow(10, self.pH)))

        # Calculate the net charge for the charged amino acid side chains
        for aa in self.sequence:
            if aa in pka_sidechain_positive:
                self.netCharge = self.netCharge + (math.pow(10, pka_sidechain_positive[aa]) / (
                            math.pow(10, pka_sidechain_positive[aa]) + math.pow(10, self.pH)))
            if aa in pka_sidechain_negative:
                self.netCharge = self.netCharge - (
                            math.pow(10, self.pH) / (math.pow(10, pka_sidechain_negative[aa]) + math.pow(10, self.pH)))

    ############################################################################
    def __calculate_properties_from_mol(self):
        """
        Function to calculate some molecular properties based on RDKit functionalities

        :return: Static physico-chemical properties like molecular weight, crippen logP, number of hydrogen bond acceptors and donors, among others
        """

        # Generate molecule from sequence
        mol = Chem.MolFromSmiles(self.smiles)
        mol.SetProp("_Name", self.sequence)

        # Calculate the descriptors
        self.num_hdonors = Lipinski.NumHDonors(mol)
        self.num_hacceptors = Lipinski.NumHAcceptors(mol)
        self.mol_weight = Descriptors.MolWt(mol)
        self.mol_logp = Crippen.MolLogP(mol)
        self.num_heteroatoms=Descriptors.NumHeteroatoms(mol)
        self.num_rotatablebonds=Descriptors.NumRotatableBonds(mol)
        self.num_heavyatoms=Descriptors.HeavyAtomCount(mol)
        self.ring_count=Descriptors.RingCount(mol)
        self.tpsa=Descriptors.TPSA(mol)


    ############################################################################
    def __calculate_properties_from_sequence(self):
        """
        Function to calculate some molecular properties based on RDKit functionalities

        :return Average Eisenberg hydrophobicity, ProtParam parameters: Isolectric point, aromaticity, instability index, amino acid percentage
        """

        # Hydrophobicity -> Eisenberg scale
        hydrophobicity = {'A': 0.620, 'R': -2.530, 'N': -0.780, 'D': -0.900, 'C': 0.290, 'Q': -0.850, 'E': -0.740,
                          'G': 0.480, 'H': -0.400, 'Y': 0.260,
                          'I': 1.380, 'L': 1.060, 'K': -1.500, 'M': 0.640, 'F': 1.190, 'P': 0.120, 'S': -0.180,
                          'T': -0.050, 'W': 0.810, 'V': 1.080}
        self.avg_hydro = sum([hydrophobicity[resi] for resi in self.sequence])

        # ProParam properties
        prot_parameters = ProteinAnalysis(self.sequence)
        self.aromaticity = prot_parameters.aromaticity()
        self.aa_percent = prot_parameters.get_amino_acids_percent()
        self.instability_index = prot_parameters.instability_index()
        self.isoelectric_point = prot_parameters.isoelectric_point()

    ############################################################################
    def __solubility_rules(self):
        """
        Function to calculate some solubility rules based on recommendations of http://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/

        :return solubility_rules_failed: return the number of rules failed based on the criteria
        """
        # Rule N1. Number of hydrophobic or charged residues
        hydro_residues = ['V', 'I', 'L', 'M', 'F', 'W', 'C']
        charged_residues = ['H', 'R', 'K', 'D', 'E']

        count_hydro_charged = 0
        for aa in self.sequence:
            if aa in hydro_residues or aa in charged_residues: count_hydro_charged += 1

        # This condition should change depending on the sequence length
        hydro_char_threshold = float(self.length_peptide) * 0.45
        if count_hydro_charged > hydro_char_threshold:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N2. Computed peptide charge
        charge_threshold = 1
        self.__compute_peptide_charges()
        if self.netCharge > 1:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N3. Glycine or Proline content in the sequence
        count_gly_pro = 0
        for aa in self.sequence:
            if aa == "G" or aa == "P": count_gly_pro += 1
        # Check threshold
        if count_gly_pro > 1:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N4. First or last amino acid charged
        count_charge = 0
        if self.sequence[0] in charged_residues:
            count_charge += 1
        if self.sequence[-1] in charged_residues:
            count_charge += 1
        # Check threshold
        if count_charge > 0:
            self.solubility_rules_failed += 1
            self.set_sol_rules.append(1)
        else:
            self.set_sol_rules.append(0)

        # Rule N5. Any amino acid represent more than 25% of the total sequence
        prot_parameters = ProteinAnalysis(self.sequence)
        aa_content = prot_parameters.get_amino_acids_percent()
        flag5 = 0
        for aa in aa_content:
            if aa_content[aa] >= 0.3:
                self.solubility_rules_failed += 1
                self.set_sol_rules.append(1)
                flag5 = 1
                break
        if flag5 == 0: self.set_sol_rules.append(0)

    ############################################################################
    def __synthesis_rules(self):
        """
        Function to check some synthesis rules based on empirical recommendations

        :return synthesis_rules_failed: return the number of rules failed based on the criteria
        """
        # Presence of forbiden motifs
        forbidden_motifs = {'2-prolines': r'[P]{3,}', 'DG-DP': r'D[GP]', 'N-Q-Nterminal': r'^[NQ]', }
        for motif in forbidden_motifs:
            if re.search(forbidden_motifs[motif], self.sequence):
                self.synthesis_rules_failed += 1
                self.set_syn_rules.append(1)
            else:
                self.set_syn_rules.append(0)

        # test if there are charged residues every 5 amino acids
        charged_residues = ['H', 'R', 'K', 'D', 'E']
        counter_charged = 0
        flag9=0
        for residue in self.sequence:
            counter_charged += 1
            if residue in charged_residues:
                counter_charged = 0
            if counter_charged >= 5:
                flag9=1
                self.synthesis_rules_failed += 1
                break
        if flag9==1:
            self.set_syn_rules.append(1)
        else:
            self.set_syn_rules.append(0)

        # Check if there are oxidation-sensitive amino acids
        aa_oxidation = ['M', 'C', 'W']
        flag10 = 0
        for aa in self.sequence:
            if aa in aa_oxidation:
                self.synthesis_rules_failed += 1
                flag10 = 1
                break
        if flag10 == 1:
            self.set_syn_rules.append(1)
        else:
            self.set_syn_rules.append(0)

    # End of Sequence class
    ############################################################

##########################################################################
# Additional functions
##########################################################################   
def blast_online(sequence):
    """
    Function to run online blast configured with parameters suitable to compare peptides

    :param sequence: Peptide sequence in FASTA format

    :return hits: List of hits with dictionary containing fields from the alignment result
    """

    with tempfile.NamedTemporaryFile(mode='w', delete=True) as temp_file:
        temp_file.write(">sequence\n{}".format(sequence))
        temp_file.flush()
      
        # Run the BLAST command
        record = SeqIO.read(temp_file.name, format="fasta")
        result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"), word_size=2, expect=20000.0,
                                        matrix_name="PAM30", gapcosts="9 1", format_object="Alignment")
        b_record = NCBIXML.read(result_handle)

        # Parse the results
        hits = []
        for alignment in b_record.alignments:
            for hsp in alignment.hsps:
                dict_hits = {}
                dict_hits["identities"] = hsp.identities
                dict_hits["positives"] = hsp.positives
                dict_hits["gaps"] = hsp.gaps
                dict_hits["align_length"] = hsp.align_length
                dict_hits["query_start"] = hsp.query_start
                dict_hits["e-value"] = hsp.expect
                dict_hits["query_sequence"] = hsp.query[0:75]
                dict_hits["match_id"] = alignment.title[:100]
                dict_hits["subject_sequence"] = hsp.sbjct[0:75]
                hits.append(dict_hits)

    return hits

##########################################################################
def get_monomer_info(path):
    """
    Convert a monomer SDF file to a Pandas dataframe object.

    :param path: os path of the monomers.sdf file
    :type path: os path

    :return: monomer dictionary as a dataframe
    """
    sdf_file = path
    df_group = PandasTools.LoadSDF(sdf_file)
    df_group = df_group.set_index('symbol')
    df_group = df_group.rename(columns={"ROMol": "m_romol"})

    return df_group

############################################################
def split_outside(string, by_element, outside, keep_marker=True):
    """
    Splits a string by delimiter only if outside of a given delimiter

    :param string: string to be split
    :param by_element: delimiter(s) by which to be split
    :param outside: only split if outside of this
    :param keep_marker: if True keep the chunk marker, remove otherwise

    :return splitChains: split string as list
    """
    # by can be more than 1 character
    by_element = list(by_element)

    # if outside is only one character (e.g. ', "), double it for start and end
    if len(outside) == 1:
        outside = outside + outside

    # Special character
    grpsep = chr(29)

    out = ''
    inside = False
    for i in string:
        if i == outside[0]:
            if inside:
                if keep_marker:
                    j = i
                else:
                    j = ''
                inside = False
            else:
                inside = True
                if keep_marker:
                    j = i
                else:
                    j = ''
        elif i == outside[1]:
            inside = False
            if keep_marker:
                j = i
            else:
                j = ''
        else:
            if not inside and i in by_element:
                j = grpsep
            else:
                j = i
        out = out + j

    # Do the final split
    split_chains = out.split(grpsep)
    return split_chains

##########################################################################
def get_peptide(biln, monomer_lib=None):
    """
    Obtain a FASTA version of the peptide. 
    Unknown NNAAs are replaced by alanine monomers.

    :param biln: BILN representation of the peptide
    :param monomer_lib: name of the monomer file

    :return: A fasta sequence of a natural peptide counterpart - str
    """
    aa_dict = {'G': 'NCC(=O)',
                'A': 'N[C@@]([H])(C)C(=O)',
                'R': 'N[C@@]([H])(CCCNC(=N)N)C(=O)',
                'N': 'N[C@@]([H])(CC(=O)N)C(=O)',
                'D': 'N[C@@]([H])(CC(=O)O)C(=O)',
                'C': 'N[C@@]([H])(CS)C(=O)',
                'E': 'N[C@@]([H])(CCC(=O)O)C(=O)',
                'Q': 'N[C@@]([H])(CCC(=O)N)C(=O)',
                'H': 'N[C@@]([H])(CC1=CN=C-N1)C(=O)',
                'I': 'N[C@@]([H])(C(CC)C)C(=O)',
                'L': 'N[C@@]([H])(CC(C)C)C(=O)',
                'K': 'N[C@@]([H])(CCCCN)C(=O)',
                'M': 'N[C@@]([H])(CCSC)C(=O)',
                'F': 'N[C@@]([H])(Cc1ccccc1)C(=O)',
                'P': 'N1[C@@]([H])(CCC1)C(=O)',
                'S': 'N[C@@]([H])(CO)C(=O)',
                'T': 'N[C@@]([H])(C(O)C)C(=O)',
                'W': 'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)',
                'Y': 'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)',
                'V': 'N[C@@]([H])(C(C)C)C(=O)'}

    # Read the monomer dataframe
    if monomer_lib is None:
        module_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(module_dir, SequenceConstants.def_path)
        file_path = os.path.join(data_dir, SequenceConstants.def_lib_filename)
        with open(file_path, 'rb') as handle:
            new_df = get_monomer_info(handle)
    else:
        module_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(module_dir, SequenceConstants.def_path)
        file_path = os.path.join(data_dir, monomer_lib)
        with open(file_path, 'rb') as handle:
            new_df = get_monomer_info(handle)

    try:
        if SequenceConstants.chain_separator in biln:
            chains=biln.split(SequenceConstants.chain_separator)
            counter_aa={}
            for i,chain in enumerate(chains):
                counter_aa[i]=0
                monomers = split_outside(chain, SequenceConstants.monomer_join, '[]')
                for mon in monomers:
                    if mon in aa_dict:
                        counter_aa[i]+=1
            chain_order = sorted(counter_aa.items(), key=lambda x: x[1], reverse=True)
            m_seq = chains[chain_order[0][0]]
        else:
            m_seq = biln
    except ValueError:
        warnings.warn(f"No main peptide was detected for peptide: {biln}")
        sys.exit(1)

    # Loop through the list of monomers of the main peptide
    total_monomers = []
    monomers = split_outside(m_seq, SequenceConstants.monomer_join, '[]')
    for res in monomers:
        mon = re.sub(r'\[', '', res)
        res = re.sub(r'\]', '', mon)
        mon = re.sub(r'\(\d+,\d+\)', '', res)
        if mon in aa_dict:
            total_monomers.append(mon)
        else:
            # Check if a natural analog is available in the dataframe
            type_mon = new_df.loc[new_df['m_abbr'] == mon, 'm_type'].item()
            if type_mon != 'cap':
                nat_analog = new_df.loc[new_df['m_abbr'] == mon, 'natAnalog'].item()
                if nat_analog != "X":
                    total_monomers.append(nat_analog)
                else:
                    # Run a similarity calculation to check the most similar AAs
                    romol = new_df.loc[new_df['m_abbr'] == mon, 'm_romol'].item()
                    mol1 = Chem.RemoveHs(romol)

                    sim_total={}
                    for aa_value in aa_dict:
                        mol2 = Chem.MolFromSmiles(aa_dict[aa_value])
                        fp1 = rdMolDescriptors.GetMorganFingerprint(mol1,2,useChirality=True)
                        fp2 = rdMolDescriptors.GetMorganFingerprint(mol2,2,useChirality=True)
                        smiles_similarity=DataStructs.DiceSimilarity(fp1,fp2)
                        sim_total[aa_value] = smiles_similarity

                    temp_order = sorted(sim_total.items(), key=lambda x: x[1])
                    aa_value = temp_order[-1][0]
                    val = float(temp_order[-1][1])
                    if val >= 0.5:
                        total_monomers.append(aa_value)
                    else:
                        # Assign an alanine if no similar AAs are found
                        total_monomers.append("A")

    fasta = ''.join(total_monomers)
    return fasta

############################################################
## End of sequence.py
############################################################