import os
from typing import Dict, List, Set
import numpy as np

from . import PDB_utils

class TriangleExtractor():
    
    _crt_path = os.path.dirname(__file__)
    _aa2subtype = {'ARG': 'B', 'LYS': 'B', 'GLU': 'J', 'ASP': 'J', 
                    'SER': 'O', 'THR': 'O', 'LEU': 'U', 'VAL': 'U', 
                    'ILE': 'U', 'GLN': 'X', 'ASN': 'X', 'TRP': 'Z', 
                    'PHE': 'Z', 'ALA': 'A', 'CYS': 'C', 'GLY': 'G', 
                    'HIS': 'H', 'MET': 'M', 'PRO': 'P', 'TYR': 'Y'}      
    
    def __init__(self, dest_path: str, verbose=False,
                 aa_index_path=None, tri_file_path=None #允许后续改动
                 ):
        self._verbose = verbose
        self._dest_path = dest_path
        # 这里有另外的aa index
        if aa_index_path is None:
            aa_index_path = os.path.join(TriangleExtractor._crt_path, 'information',
                                         'tri_AA_index.txt')
        self._aa_indexes = PDB_utils.get_aa_index(aa_index_path)
        # 初始化三角类别
        if tri_file_path is None:
            tri_file_path = os.path.join(TriangleExtractor._crt_path, 'information',
                                         '455_list.txt')
        self.triangles_in_subtype = []
        for line in open(tri_file_path):
            self.triangles_in_subtype.append(line.strip())

    
    def _get_triangles(self, resds_atoms: Dict[str, List[List[float]]],
                       surface_aas: Set[str]) -> Dict[str, List[List[str]]]:
        '''
        在一个链内一位位氨基酸的计算三角
        '''
        # aa—info：和他距离小于4埃的氨基酸列表
        dist_below4_aas = {}
        # 维护一个dict，包含所有距离小于4埃的氨基酸对,后续不用反复算了
        calced_dist_pairs = {}
        # 整个链中各个氨基酸对应的三角
        aas_triangles = {}
        for aa_info1, resd_atoms1 in resds_atoms.items():
            # 结果中不考虑非表面氨基酸
            if aa_info1 not in surface_aas or aa_info1.count('UNK') != 0:
                continue
            
            dist_below4_aas[aa_info1] = []
            for aa_info2, resd_atoms2 in resds_atoms.items():
                if aa_info1 == aa_info2 or aa_info2.count('UNK') != 0:
                    continue
                # 算过，就不再次算了
                if (aa_info1, aa_info2) in calced_dist_pairs:
                    if calced_dist_pairs[(aa_info1, aa_info2)]:
                        dist_below4_aas[aa_info1].append(aa_info2)

                else: # 没算过，先算后记录
                    dist_below_thresh = PDB_utils.resd_dist_between_thresh(resd_atoms1,
                                                                           resd_atoms2)
                    calced_dist_pairs[(aa_info2, aa_info1)] = dist_below_thresh
                    calced_dist_pairs[(aa_info1, aa_info2)] = dist_below_thresh
                    if dist_below_thresh:
                        dist_below4_aas[aa_info1].append(aa_info2)
        for aa_info in resds_atoms:
            # 检查4埃内的氨基酸有没有可以凑出三角的
            aas_triangles[aa_info] = []
            choosen_aas = dist_below4_aas[aa_info]
            for idx in range(len(choosen_aas)):
                choosen_aa1 = choosen_aas[idx]
                for choosen_aa2 in choosen_aas[idx+1: ]:
                    if calced_dist_pairs[(choosen_aa1, choosen_aa2)]:
                        aas_triangles[aa_info].append([choosen_aa1, choosen_aa2])      
        return aas_triangles
    
    
    def _output_triangle_file(self, pdb_name: str,
                              aa_triangles: Dict[str, List[List[str]]]):
        triangle_folder = os.path.join(self._dest_path, '3.2_1_Triangle')
        if not os.path.exists(triangle_folder):
            os.mkdir(triangle_folder)
        dest_path = os.path.join(triangle_folder, f'Triangle_{pdb_name}.txt')
        dest_file = open(dest_path, 'w')
        for aa_key, triangles in aa_triangles.items():
            dest_file.write(f'{aa_key}\t')
            triangle_list = []
            for triangle in triangles:
                triangle_list.append(','.join(triangle))
            dest_file.write(f"{'|'.join(triangle_list)}\n")
        dest_file.close()
          
                                    
    def _triangle2subtypes(self, aas_triangles: Dict[str, List[List[str]]]
                           ) -> Dict[str, List[str]]:
        aas_tri_subtypes = {}
        for aa_info1, triangles in aas_triangles.items():
            aa_type = aa_info1[:3]
            triangles_as_subtype = []
            for triangle in triangles:
                crt_subtypes = [
                    TriangleExtractor._aa2subtype[aa_info2[:3]] 
                    for aa_info2 in triangle
                ]
                crt_subtypes.append(TriangleExtractor._aa2subtype[aa_type])
                triangles_as_subtype.append("_".join(sorted(crt_subtypes)))
            aas_tri_subtypes[aa_info1] = triangles_as_subtype
        return aas_tri_subtypes 
    
    def _output_tri_subtype_file(self, pdb_name: str,
                                 aas_tri_subtypes: Dict[str, List[str]]):
        tri_subtype_folder = os.path.join(self._dest_path, '3.2_2_Triangle_to_subtypes')
        if not os.path.exists(tri_subtype_folder):
            os.mkdir(tri_subtype_folder)
        dest_path = os.path.join(tri_subtype_folder, f'Triangle_{pdb_name}.txt')
        dest_file = open(dest_path, 'w')
        for aa_key, tri_in_subtypes in aas_tri_subtypes.items():
            dest_file.write(f'{aa_key}\t')
            dest_file.write(f'{"|".join(tri_in_subtypes)}\n')
        dest_file.close()
                
    
    def _count_aa_triangle(self, aas_tri_subtypes: Dict[str, List[str]]
                          ) -> Dict[str, List[int]]:
        # 统计455种triangle在每一个目标氨基酸中的分布
        all_triangle_subtype_cnts = {}
        for aa_info, triangles_in_subtype in aas_tri_subtypes.items():
            tri_subtype_cnts = [
                triangles_in_subtype.count(subtype) 
                for subtype in self.triangles_in_subtype
            ]
            all_triangle_subtype_cnts[aa_info] = tri_subtype_cnts
        return all_triangle_subtype_cnts
    
    
    def _calc_mean_aa_index(self, aas_triangles: Dict[str, List[List[str]]]
                           ) -> Dict[str, List[float]]:
        # 对每一位氨基酸统计可以和他组成三角的其他氨基酸
        # 全部算清楚之后计算每一种aa_index的平均
        all_avg_aa_indexes = {}
        for aa_info1, triangles in aas_triangles.items():
            all_aas = [aa_info1]
            for triangle in triangles:
                all_aas += triangle
            all_aas = [aa_info2[:3] for aa_info2 in set(all_aas)]  
            crt_aa_indexes_vals = []
            for aa_index in self._aa_indexes:
                aa_index_avg = np.mean([aa_index[aa_type] 
                                        for aa_type in all_aas])
                crt_aa_indexes_vals.append(aa_index_avg)
            all_avg_aa_indexes[aa_info1] = crt_aa_indexes_vals
        return all_avg_aa_indexes
                
    def _output_full_descriptor(self, mean_aa_indexes: Dict[str, List[float]],
                                triangle_subtype_counts:  Dict[str, List[int]],
                                pdb_name: str) -> Dict[str, List[float]]:
        # 偷个懒，计算和输出写在一起了
        full_triangle_descriptors = {}

        for aa_info, mean_aa_index in mean_aa_indexes.items():

            triangle_subtype_cnt = triangle_subtype_counts[aa_info]
            full_descriptor = mean_aa_index + triangle_subtype_cnt
            full_triangle_descriptors[aa_info] = full_descriptor
        if self._verbose:
            descriptor_folder = os.path.join(self._dest_path, '3.2_3_AA_descriptors_Triangles')
            if not os.path.exists(descriptor_folder):
                os.mkdir(descriptor_folder)
            dest_path = os.path.join(descriptor_folder, f'{pdb_name}.txt')
            dest_file = open(dest_path, 'w')
            for aa_info, full_descriptor in full_triangle_descriptors.items():
                aa_type, chain, aa_number = aa_info.split()
                full_descriptor = [str(f'{x:.4f}') for x in full_descriptor]
                dest_file.write(f'{aa_type}_{aa_number}_{chain}\t{" ".join(full_descriptor)}\n')
            dest_file.close()
        return full_triangle_descriptors
    
    
    def get_descriptors(self, interested_chains: List[str], 
                       surface_aas: Set, pdb_folder: str, 
                       pdb_name: str) -> Dict[str, List[float]]:
        pdb_path = os.path.join(pdb_folder, f'{pdb_name}.pdb')
        # 不想单独出一个函数了，而且也便于后面迁移到更多链，就这么做了hhh
        all_aa_triangles = {}
        for chain in interested_chains:
            resds_atoms = PDB_utils.get_resds_atoms(pdb_path, [chain])
            triangles = self._get_triangles(resds_atoms, surface_aas)
            all_aa_triangles.update(triangles)
        if self._verbose:
            self._output_triangle_file(pdb_name, all_aa_triangles)
        # 每个氨基酸附近存在的三角转化为subtype
        all_tri_subtypes = self._triangle2subtypes(all_aa_triangles)
        if self._verbose:
            self._output_tri_subtype_file(pdb_name, all_tri_subtypes)
        
        # 生成三角描述符的两个部分
        triangle_subtype_counts = self._count_aa_triangle(all_tri_subtypes)
        mean_aa_indexes = self._calc_mean_aa_index(all_aa_triangles)
        
        # 合并描述符，写入文件
        full_descriptor = self._output_full_descriptor(mean_aa_indexes, 
                                                       triangle_subtype_counts,
                                                       pdb_name)
        return full_descriptor
            
                
            