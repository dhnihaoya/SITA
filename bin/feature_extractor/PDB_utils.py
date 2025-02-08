import os
import numpy as np
from typing import Dict, List, Set
import freesasa
import shutil

def is_carbon_only_pdb(file_path):
    """
    检查PDB文件是否所有氨基酸原子只包括碳原子。
    """
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_type = line[76:78].strip()
                if atom_type != 'C':
                    return False
    return True

def orgnize_pdb_by_chains(PDB_folder: str, pdb_chain_fp: str, 
                                            # pdb_id, [chains]
                        new_pdb_path: str) -> List:
    chain_information = open(pdb_chain_fp)
    complex_informations = []
    for line in chain_information:
        line_split = line.strip().split(",")
        PDB_id = line_split[0]
        chains = line_split[1:]
        complex_informations.append((PDB_id, chains))
        dest_file_path = os.path.join(new_pdb_path, f"{PDB_id}_{'_'.join(chains)}.pdb")
        new_PDB_file = open(dest_file_path, 'w')
        
        original_PDB = os.path.join(PDB_folder, f'{PDB_id}.pdb')
        ## 判断pdb文件是否仅包括碳原子，如果仅包括碳原子 直接跳过文件
        if not is_carbon_only_pdb(original_PDB):
            for line in open(original_PDB):
                if not line.startswith('ATOM'):
                    continue
                chain = line[21]
                if chain in chains:
                    new_PDB_file.write(line)
            new_PDB_file.close()
    return complex_informations


def pdb2asa(pdb_folder: str, pdb_file_name: str, dest_folder: str):
    input_file = os.path.join(pdb_folder, pdb_file_name)
    structure = freesasa.Structure(input_file)
    sasa_result = freesasa.calc(structure)
    asa_file_name = pdb_file_name.replace('.pdb', '.asa')
    sasa_result.write_pdb(os.path.join(dest_folder, asa_file_name))

        
def get_surface_aas(asa_path: str) -> Set:
    surface_aas = set()
    for atom_line in open(asa_path):
        if not atom_line.startswith("ATOM"):
            continue
                
        aa_info = atom_line[17: 30].strip()
        surface_aas.add(aa_info)
    return surface_aas
        
def get_resds_atoms(pdb_path: str, ag_chains: List[str]
                  # 氨基酸信息：氨基酸每个原子的坐标
                  ) -> Dict[str, List[List[float]]]:
    residue_atoms = {}
    with open(pdb_path) as pdb_f:
        while 1:
            line = pdb_f.readline()
            if not line:
                break
            if not line.startswith('ATOM'):
                continue
            addinfo = line[17:30].strip()
            chain_id = line[21]
                    
            atom_coord = [float(line[30:38]), float(line[38:46]), 
                          float(line[46:54])]
            if chain_id in ag_chains:
                if addinfo not in residue_atoms:
                    residue_atoms[addinfo] = []
                residue_atoms[addinfo].append(atom_coord)
                            
    return residue_atoms

def get_aa_index(aa_idx_file_path) -> List[Dict[str, float]]:
    # 顺序定了就不要再改
    AA_types = ['ALA', 'LEU', 'ARG', 'LYS', 'ASN', 
                 'MET', 'ASP', 'PHE', 'CYS', 'PRO',
                 'GLN', 'SER', 'GLU', 'THR', 'GLY', 
                 'TRP', 'HIS', 'TYR', 'ILE', 'VAL']

    aa_idx_file = open(aa_idx_file_path)
    aa_indexes = []
    for line in aa_idx_file:
        line = line.strip().split('\t')
        aa_index = [float(x) for x in line[-1].split()]
        # 这一种aa index中，各种氨基酸对应的值
        aa_data_dict = {}
        for aa_type, aa_idx_val in zip(AA_types, aa_index):
            aa_data_dict[aa_type] = aa_idx_val
        aa_indexes.append(aa_data_dict)
    return aa_indexes
       
        
def calc_eucliean_distance(cords1: List[float], 
                           cords2: List[float]) -> float:
    cords1 = np.array(cords1)
    cords2 = np.array(cords2)
    dist = np.linalg.norm(cords2 - cords1)
    return dist


# 用np计算替代循环，这样更快
def resd_dist_between_thresh(resd_atoms_cords1: List[List[float]],
                             resd_atoms_cords2: List[List[float]],
                             dist_thresh: int=4) -> bool:
    
    resd_atoms_cords1 = np.array(resd_atoms_cords1).reshape([-1, 1, 3])
    resd_atoms_cords2 = np.array(resd_atoms_cords2).reshape([1, -1, 3])
    
    dists = np.sqrt(np.sum((resd_atoms_cords1 - resd_atoms_cords2) ** 2, axis=2))
    
    return np.any(dists <= dist_thresh)      
             
    
            