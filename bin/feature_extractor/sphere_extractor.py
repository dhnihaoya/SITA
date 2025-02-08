import os
from typing import Dict, List, Set
import numpy as np

from . import PDB_utils

class SphereExtractor():
    '''
    用于根据aa index文件和pdb文件提取整个结构中每一位氨基酸的球壳描述符
    如果需要换一套aa描述符, 根据文件格式重新组织一下之后这个文件还能继续用
    '''
    _crt_path = os.path.dirname(__file__)
    
    def __init__(self, dest_path: str, verbose=False, 
                 distance_ubs=[2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
                 aa_index_path=None): #允许重新定义aa index file，但是需要用和原来的文件类似的方法组织
        self._verbose = verbose
        self._distance_ubs = distance_ubs
        self.dest_path = dest_path
        # 载入aa_index, 没定义就用之前做好了的
        if aa_index_path is None:
            aa_index_path = os.path.join(SphereExtractor._crt_path, 'information', 
                                         'sphere_AA_index.txt')
        self._aa_indexes = PDB_utils.get_aa_index(aa_index_path)
            
        
    def _calc_sphere_layer(self, resds_atoms: Dict[str, List[List[float]]],
                           surface_aas: Set[str]) -> Dict[str, List[List[str]]]: # aa信息：各个球壳层中包含的氨基酸种类
        '''
        计算整个结构中每一位alpha碳之间的距离, 然后按照各层球壳的距离标准计算
        这两位残基属于哪个范围.
        结果中应该只包含表面氨基酸的球壳情况
        返回每一位氨基酸对应的每层球壳包含了哪些氨基酸
        道理来说可以用np改写继续提升速度，但我懒得了。我感觉这样还能好读点
        '''
        already_calced_layers = {} # ab间距=ba间距。算完就记下来在哪层，不再计算第二遍
        sphere_results = {}
        for aa_info in resds_atoms:
            # 最后一个list中是超过上限的aa，后续计算完之后不做保留
            sphere_results[aa_info] = [[] for _ in range(len(self._distance_ubs)+1)]
        
        for aa_info1, resd_cords1 in resds_atoms.items():
            crt_aa_spheres = sphere_results[aa_info1]
            
            if (len(resd_cords1) == 1):
                continue
            if (aa_info1 not in surface_aas) or (aa_info1.count("UNK")!=0):
                continue
            
            for aa_info2, resd_cords2 in resds_atoms.items():
                if len(resd_cords2) == 1 or aa_info2.count("UNK") != 0:
                    continue
                # 换个顺序之前算过了。查表，不再计算第二遍
                if (aa_info1, aa_info2) in already_calced_layers:
                    if aa_info1 == aa_info2:
                        continue
                    sphere_layer = already_calced_layers[(aa_info1, aa_info2)]
                    crt_aa_spheres[sphere_layer].append(aa_info2)
                    continue
                
                C_alpha_cords1 = resd_cords1[1]
                C_alpha_cords2 = resd_cords2[1]
                C_alpha_dist = PDB_utils.calc_eucliean_distance(C_alpha_cords1, 
                                                                C_alpha_cords2) 
                # 用这个可以不再逐层检查, 直接算出在第几层的范围里
                layer = np.digitize(C_alpha_dist, self._distance_ubs)
                crt_aa_spheres[layer].append(aa_info2)
                
                already_calced_layers[(aa_info2, aa_info1)] = layer
            # 删去最后一层
            sphere_results[aa_info1] = crt_aa_spheres[:-1]     
        return sphere_results

                               
    def _output_sphere_file(self, pdb_name: str, 
                            sphere_results: Dict[str, List[List[str]]]):
        sphere_folder = os.path.join(self.dest_path, '3.1_1_Spherical_AA_model')
        if not os.path.exists(sphere_folder):
            os.mkdir(sphere_folder)
        dest_path = os.path.join(sphere_folder, f'{pdb_name}_Spherical_AA_list.txt')
        dest_file = open(dest_path, 'w')
        for aa_info, aa_spheres in sphere_results.items():
            aa_type, _, aa_number = aa_info.split()
            dest_file.write(f"{aa_number}_{aa_type}\t")
            aa_cnt_total = 0
            sphere_aas_by_layer = [] # 每层球壳含有的aa空格链接拼接为字符串后加入列表
            for sphere in aa_spheres:
                aa_cnt_total += len(sphere)
                sphere_aas_by_layer.append(" ".join(sphere))
            crt_aa_spheres = " | ".join(sphere_aas_by_layer)
            dest_file.write(f'{aa_cnt_total}\t{crt_aa_spheres}\n')
        dest_file.close()
    
        
# 绷不住了 怎么还是写4层循环
    def _spheres2descriptors(self, spheres: Dict[str, List[List[str]]]
                            ) -> Dict[str, List[float]]:
        descriptors = {}
        aa_index_cnt = len(self._aa_indexes)
        
        for aa_info1, sphere_layers in spheres.items():
            crt_aa_descriptor = []
            # 一层层计算描述符
            for sphere_layer in sphere_layers:
                # 如果整层是空的，则这层描述符全部设为0
                if sphere_layer == []:
                    crt_layer_descriptor = [0 for _ in range(aa_index_cnt)]
                    crt_aa_descriptor.extend(crt_layer_descriptor)
                    continue
                # 如果这层不是空的，那么依次计算每一个aa index
                crt_layer_descriptor = []
                for aa_index_vals in self._aa_indexes:
                    score = 0
                    for aa_info2 in sphere_layer:
                        aa_type = aa_info2.split()[0]
                        score += aa_index_vals[aa_type]
                    crt_layer_descriptor.append(score)
                crt_aa_descriptor.extend(crt_layer_descriptor)
            descriptors[aa_info1] = crt_aa_descriptor
        return descriptors    
        
                                    
    def _output_descriptor_file(self, pdb_name: str, 
                                sphere_descriptors: Dict[str, List[float]]):
        sphere_folder = os.path.join(self.dest_path, '3.1_2_Spherical_AA_descriptors')
        if not os.path.exists(sphere_folder):
            os.mkdir(sphere_folder)
        dest_path = os.path.join(sphere_folder, 
                                 f'descriptor_{pdb_name}_Spherical_AA_list.txt')
        dest_file = open(dest_path, 'w')
        for aa_info, sphere_descriptors in sphere_descriptors.items():
            aa_type, _, aa_number = aa_info.split()
            dest_file.write(f'{aa_number}_{aa_type}\t')
            formated_descriptors = [str(format(x, '.4f')) for x in sphere_descriptors]
            dest_file.write('   '.join(formated_descriptors)+'\n')
        dest_file.close()
            
    
    def get_descriptors(self, interested_chains: List[str], 
                        surface_aas: Set, 
                        # folder和文件名分开给，处理比较方便
                        # pdb_name包含了{pdb_id}_{H}_{L}
                        pdb_folder: str, pdb_name: str) -> Dict[str, List[float]]:
        # 获取整个结构中每个aa和对应的各个原子坐标
        pdb_file_path = os.path.join(pdb_folder, f'{pdb_name}.pdb')
        resds_atoms = PDB_utils.get_resds_atoms(pdb_file_path, interested_chains)
        
        spheres_info = self._calc_sphere_layer(resds_atoms, surface_aas)
        if self._verbose and self.dest_path is not None:
            self._output_sphere_file(pdb_name, spheres_info)
        sphere_descriptors = self._spheres2descriptors(spheres_info)
        if self._verbose:
            self._output_descriptor_file(pdb_name, sphere_descriptors)

        return sphere_descriptors

                    
                     
            
            
            
                
            