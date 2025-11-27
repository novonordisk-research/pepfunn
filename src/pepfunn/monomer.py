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
import pandas as pd
import os
import itertools
import warnings
import sys
import re
from string import ascii_uppercase as alc
import pickle
from pathlib import Path

# RDKit
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import SDWriter
from rdkit.Chem import AllChem

########################################################################################
# Classes and functions
########################################################################################

class Monomer:

    '''Class with functions to manipulate monomers'''

    ########################################################################################
    @staticmethod
    def get_abbreviation(iupac_name):
        '''
        Accept an IUPAC name and return an abbreviation suitable for monomer name
        '''
        # Remove any parenthetical information from the IUPAC name
        iupac_name = re.sub(r'\[.*?\]', '', iupac_name)
        iupac_name = re.sub(r'\(.*?\)', '', iupac_name)

        # Split the IUPAC name into words
        words = iupac_name.split()

        # Remove any common words and abbreviations
        common_words = ['and', 'of', 'in', 'the', 'with', 'for', 'on', 'as', 'an', 'at']
        abbreviations = ['meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec']
        words = [word for word in words if word.lower() not in common_words and word.lower() not in abbreviations]

        final_bag=[]
        for word in words:
            if '-' in word:
                list_w=word.split('-')
                for part in list_w:
                    if part:
                        if len(list(part))>4:
                            final_bag.append(part[0:4].capitalize())
                        else:
                            if len(part)==1:
                                final_bag.append(part)
                            else:
                                final_bag.append(part.capitalize())
            else:
                if len(list(word))>4:
                    final_bag.append(word[0:4])
                else:
                    if len(word)==1:
                        final_bag.append(word)
                    else:
                        final_bag.append(word.capitalize())

        abbreviation='-'.join(final_bag)
        
        return abbreviation
    
    ########################################################################################
    @staticmethod
    def get_natural_analog(smiles, sim_threshold=0.5):
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
        
        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=4,fpSize=1024) 
        romol = Chem.MolFromSmiles(smiles)
        mol1 = Chem.RemoveHs(romol)
        
        sim_total={}
        for aa_value in aa_dict:
            mol2 = Chem.MolFromSmiles(aa_dict[aa_value])
            fp1 = mfpgen.GetFingerprint(mol1)
            fp2 = mfpgen.GetFingerprint(mol2)
            
            smiles_similarity=DataStructs.DiceSimilarity(fp1,fp2)
            sim_total[aa_value] = smiles_similarity
        
        temp_order = sorted(sim_total.items(), key=lambda x: x[1])
        aa_value = temp_order[-1][0]
        val = float(temp_order[-1][1])
        if val >= sim_threshold:
            nat_analog=aa_value
        else:
            nat_analog='X'

        return nat_analog, sim_total
    
    ########################################################################################
    @staticmethod
    def generate_PDB_code(symbol, comp_pdb, list_codes):
        totalchar = alc + '0123456789'
        monomers = {"A": "ALA", "D": "ASP", "E": "GLU", "F": "PHE", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
                "M": "MET", "G": "GLY", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR",
                "V": "VAL", "W": "TRP", "Y": "TYR", "C": "CYS", "ac": "ACE", "Aib": "AIB", "am": "NH2", "Iva": "6ZS"}

        if symbol in monomers:
            pdb_code = monomers[symbol]
        else:
            newSymbol_pre = symbol.replace('_','')
            newSymbol = newSymbol_pre.replace('-','')

            if len(newSymbol)>=3:
                pdb_code = newSymbol[:3].upper()
                count = 0
                attempts = []
                attempts_extra = []
                while (pdb_code in list_codes) or (pdb_code in comp_pdb):

                    if count >= len(totalchar):
                        for c1 in range(0,len(totalchar)):
                            for c2 in range(0,len(totalchar)):
                                attempts.append(pdb_code[:1] + totalchar[c1] + totalchar[c2])

                        for a in attempts:
                            if a not in list_codes and a not in comp_pdb:
                                pdb_code = a
                                break

                        if (pdb_code in list_codes) or (pdb_code in comp_pdb):
                            for c1 in range(0, len(totalchar)):
                                for c2 in range(0, len(totalchar)):
                                    for c3 in range(0, len(totalchar)):
                                        attempts_extra.append(totalchar[c1] + totalchar[c2] + totalchar[c3])

                            for a in attempts_extra:
                                if a not in list_codes and a not in comp_pdb:
                                    pdb_code = a
                                    break
                    else:
                        pdb_code = pdb_code[:2] + totalchar[count]
                        count += 1
            else:
                pdb_code = newSymbol.upper()
                attempts = []
                if len(pdb_code) == 2:
                    for ch in totalchar:
                        new_code = pdb_code+ch
                        if new_code not in list_codes and new_code not in comp_pdb:
                            pdb_code = new_code
                            break
                    if len(pdb_code) == 2:
                        for c1 in range(0,len(totalchar)):
                            for c2 in range(0,len(totalchar)):
                                attempts.append(pdb_code[:1] + totalchar[c1] + totalchar[c2])
                        for a in attempts:
                            if a not in list_codes and a not in comp_pdb:
                                pdb_code = a
                                break
            print(pdb_code)

        return pdb_code
    
    ########################################################################################
    @staticmethod
    def generate_SDF(monomers, output_name='monomer_custom.sdf'):
        '''
        Generate SDF file of the monomers
        '''

        compList=[x.strip() for x in open("comp_categories_name.txt")]
        comp_PDB=[]
        for c in compList:
            info=c.split()
            code=info[0]
            comp_PDB.append(code)
        
        list_codes = ["ALA", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU", "MET", "GLY", "ASN", "PRO", "GLN", "ARG", "SER", "THR",
                    "VAL", "TRP", "TYR", "CYS", "ACE", "AIB", "NH2", "6ZS"]
        writer = SDWriter(output_name)

        for mon in monomers:
            mol = monomers[mon][0]
            naturalA = monomers[mon][2]
            attachments = monomers[mon][3]
            name = monomers[mon][4]
            mType = monomers[mon][5]
            full_smiles = monomers[mon][6]
            symbol = mon
            
            order = list(range(mol.GetNumAtoms()))
            indices = []
            for at in mol.GetAtoms():
                if at.GetSymbol()[0]=='*':
                    indices.append(at.GetIdx())
                    order.remove(at.GetIdx())
        
            for j in indices: order.append(j)
            newMol = Chem.RenumberAtoms(mol, newOrder=order)
            mol=newMol

            r1 = None; r2=None; r3=None; r4=None
            for att in attachments:
                if att[0]=='1':
                    r1=attachments[att]
                if att[0]=='2':
                    r2=attachments[att]
                if att[0]=='3':
                    r3=attachments[att]
                if att[0]=='4':
                    r4=attachments[att]
            
            rGroupIdx = [None, None, None, None]
            attachmentIdx = [None, None, None, None]
            rGroups=[r1,r2,r3,r4]

            count_anchor=1
            for at in mol.GetAtoms():
                
                if at.GetSymbol()[0]=='*':
                    #label=int(at.GetProp('molAtomMapNumber'))
                    label=count_anchor
                    rootAtom = at.GetNeighbors()[0]
                    rGroupIdx[label-1]=at.GetIdx()
                    attachmentIdx[label-1]=rootAtom.GetIdx()
                    count_anchor+=1
            
            pdb_code = Monomer.generate_PDB_code(mon, comp_PDB, list_codes)
            
            mol.SetProp('m_name', name)
            mol.SetProp('symbol', symbol)
            mol.SetProp('m_abbr', symbol)
            if full_smiles:
                mol.SetProp('m_smiles', full_smiles)
            else:
                mol.SetProp('m_smiles', 'None')
            mol.SetProp('m_type', mType)
            mol.SetProp('m_subtype', mType)
            new_rGroups=[str(x) for x in rGroups]
            mol.SetProp('m_Rgroups', ','.join(new_rGroups))
            new_rGroupIdx=[str(x) for x in rGroupIdx]
            mol.SetProp('m_RgroupIdx', ','.join(new_rGroupIdx))
            new_attachmentIdx=[str(x) for x in attachmentIdx]
            mol.SetProp('m_attachmentPointIdx', ','.join(new_attachmentIdx))
            if naturalA=='None':
                mol.SetProp('natAnalog', 'X')
            else:
                mol.SetProp('natAnalog', naturalA)

            mol.SetProp('pdbName', pdb_code)
            writer.write(mol)

        writer.close()

    ########################################################################################
    @staticmethod
    def create_boltz_ccd_record(name_mon, romol, reorder_atom=False, update_ccd=False):
        
        # Pickle conditions
        Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

        # Load ccd.pkl
        if update_ccd:
            ccd_path = Path('ccd.pkl')
            with ccd_path.open("rb") as file:
                ccd = pickle.load(file)

        # Boltz metadata
        output_ccd_code=name_mon
        smiles_mol=romol

        # Assign leaving group
        for idx, atom in enumerate(smiles_mol.GetAtoms()):
            info = atom.GetPDBResidueInfo()
            name = info.GetName().strip()

            # Set atom properties
            atom.SetProp('name', name)
            atom.SetProp('alt_name', name) 
            is_leaving = False
            if name == 'OXT':
                is_leaving = True
            atom.SetBoolProp('leaving_atom', is_leaving)

        # Reorder atoms
        mol = AllChem.RemoveHs(smiles_mol)
        
        if reorder_atom:
            curr_atom_order = {atom.GetProp('name'): idx for idx, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() != 'H'}

            list_index=[]
            list_after=[]
            for atom in curr_atom_order:
                if atom in ['N', 'C', 'CA', 'O', 'OXT']:
                    list_index.append(curr_atom_order[atom])
                else:
                    list_after.append(curr_atom_order[atom])

            list_index=list_index+list_after

            mol = Chem.RenumberAtoms(mol,newOrder=list_index)

        try:
            mol = AllChem.AddHs(mol)
            AllChem.EmbedMolecule(mol, maxAttempts=5000, randomSeed=0xF00D)
            AllChem.UFFOptimizeMolecule(mol)
        except:
            pass

        # Remove hydrogens
        trim_reordered = AllChem.RemoveHs(mol)

        # Set conformer properties
        for c in trim_reordered.GetConformers():
            c.SetProp("name", "Ideal")

        #Save file
        if update_ccd:
            ccd[output_ccd_code] = trim_reordered
            ccd_path = Path('ccd.pkl')
            with ccd_path.open("wb") as file:
                pickle.dump(ccd, file)

        ccd_path = Path(f'{output_ccd_code}.pkl')
        with ccd_path.open("wb") as file:
            pickle.dump(trim_reordered, file)
            

    ########################################################################################
    @staticmethod
    def create_boltz_input(sequences, id_values=[], modifications=[], folder_name='configurations'):
        '''
        Function to create fasta files in a new folder
        '''
        os.makedirs(folder_name, exist_ok=True)

        names=[]
        ids= open('ids.txt','w')
        counter=1
        for i,s in enumerate(sequences):
            if id_values:
                name=id_values[i]
            else:
                name=f'seq_{counter}'
            output = open(f'{folder_name}/input_{name}.yaml','w')
            output.write(f'version: 1\n')
            output.write(f'sequences:\n')
            output.write(f'  - protein:\n')
            output.write(f'      id: A\n')
            output.write(f'      sequence: {s}\n')
            if modifications:
                if modifications[i]:
                    for j,key in enumerate(modifications[i]):
                        if j==0:
                            output.write(f'      modifications:\n')
                        output.write(f'        - position: {key}\n')
                        output.write(f'          ccd: {modifications[i][key]}\n')
            output.close()
            names.append(f'{name}')
            ids.write(f'{name}\n')
            counter+=1

        return names

    ########################################################################################
    @staticmethod
    def create_boltz_complex(sequences, target_sequence, target_sequence2='', id_values=[], modifications=[], folder_name='configurations'):
        '''
        Function to create fasta files in a new folder
        '''
        os.makedirs(folder_name, exist_ok=True)
        chains = ['A','B','C','D','E']
        names=[]
        ids= open('ids.txt','w')
        counter=1
        for i,s in enumerate(sequences):
            chain_counter=0
            if id_values:
                name=id_values[i]
            else:
                name=f'seq_{counter}'
            output = open(f'{folder_name}/input_{name}.yaml','w')
            output.write(f'version: 1\n')
            output.write(f'sequences:\n')
            output.write(f'  - protein:\n')
            output.write(f'      id: {chains[chain_counter]}\n')
            output.write(f'      sequence: {target_sequence}\n')
            chain_counter+=1
            if target_sequence2:
                output.write(f'  - protein:\n')
                output.write(f'      id: {chains[chain_counter]}\n')
                output.write(f'      sequence: {target_sequence2}\n')
                chain_counter+=1
            output.write(f'  - protein:\n')
            output.write(f'      id: {chains[chain_counter]}\n')
            output.write(f'      sequence: {s}\n')
            chain_counter+=1
            if modifications:
                if modifications[i]:
                    for j,key in enumerate(modifications[i]):
                        if j==0:
                            output.write(f'      modifications:\n')
                        output.write(f'        - position: {key}\n')
                        output.write(f'          ccd: {modifications[i][key]}\n')
            output.close()
            names.append(f'{name}')
            ids.write(f'{name}\n')
            counter+=1

        return names

    ########################################################################################
    @staticmethod
    def create_boltz_input_smiles(smiles, id_values=[], folder_name='configurations'):
        '''
        Function to create fasta files in a new folder
        '''
        os.makedirs(folder_name, exist_ok=True)

        names=[]
        ids= open('ids.txt','w')
        counter=1
        for i,s in enumerate(smiles):
            if id_values:
                name=id_values[i]
            else:
                name=f'seq_{counter}'
            output = open(f'{folder_name}/input_{name}.yaml','w')
            output.write(f'version: 1\n')
            output.write(f'sequences:\n')
            output.write(f'  - ligand:\n')
            output.write(f'      id: A\n')
            output.write(f'      smiles: {s}\n')
            output.close()
            names.append(f'{name}')
            ids.write(f'{name}\n')
            counter+=1

        return names
    
    ########################################################################################
    @staticmethod
    def create_boltz_complex_smiles(smiles, target_sequence, id_values=[], affinity_flag=True, folder_name='configurations'):
        '''
        Function to create fasta files in a new folder
        '''
        os.makedirs(folder_name, exist_ok=True)

        names=[]
        ids= open('ids.txt','w')
        counter=1
        for i,s in enumerate(smiles):
            if id_values:
                name=id_values[i]
            else:
                name=f'seq_{counter}'
            output = open(f'{folder_name}/input_{name}.yaml','w')
            output.write(f'version: 1\n')
            output.write(f'sequences:\n')
            output.write(f'  - protein:\n')
            output.write(f'      id: A\n')
            output.write(f'      sequence: {target_sequence}\n')
            output.write(f'  - ligand:\n')
            output.write(f'      id: B\n')
            output.write(f'      smiles: {s}\n')
            if affinity_flag:
                output.write(f'properties:\n')
                output.write(f'  - affinity:\n')
                output.write(f'      binder: B\n')

            output.close()
            names.append(f'{name}')
            ids.write(f'{name}\n')
            counter+=1

        return names
    
