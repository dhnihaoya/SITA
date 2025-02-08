import os
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from tabulate import tabulate
from typing import Dict, List

# 从PDB结构中提取轻链重链PDB
from .PDB_utils import *
from .sphere_extractor import SphereExtractor
from .triangle_extractor import TriangleExtractor

class FeatureExtractor:
    _crt_path = os.path.dirname(__file__)
    
    def __init__(self, verbose=False, dest_path=None):
        self._verbose = verbose # output more detailed information
        self._dest_path = dest_path 
        self._sphere_extractor = SphereExtractor(dest_path, verbose)
        self._triangle_extractor = TriangleExtractor(dest_path, verbose)

    def _concat_descriptors(self, sphere_feas: Dict[str, List[float]],
                            triangle_feas: Dict[str, List[float]]) -> Dict[str, List[float]]:
        result = {}
        for aa_info, sphere_fea in sphere_feas.items():
            if aa_info not in triangle_feas:
                continue
            
            triangle_fea = triangle_feas[aa_info]
            result[aa_info] = sphere_fea + triangle_fea
            
        return result 
            
    
    def get_features(self, pdb_folder: str, chain_file_path: str,
                    dest_folder: str) -> Dict[str, Dict[str, List[float]]]:
        # 根据chain文件重新组织pdb文件内容，只保留需要的链
        orgnized_pdb_folder = os.path.join(dest_folder, '0.1_1_orgnized_PDB')
        if not os.path.exists(orgnized_pdb_folder):
            os.mkdir(orgnized_pdb_folder)
        antibody_complexs = orgnize_pdb_by_chains(pdb_folder, 
                                                  chain_file_path, 
                                                  orgnized_pdb_folder)
        # 如果verbose, 输出更为详细的过程内容
        if self._verbose:
            print('About to predict the immunogenicity of these antibody:')
            complex_information = [
                [antibody_complex[0], ','.join(antibody_complex[1])] 
                for antibody_complex in antibody_complexs
            ]
            header = ["PDB_id", "chains"]
            print(tabulate(complex_information, headers=header, tablefmt="grid"))
            
            print('Checking surface aa residues with NACCESS:')
        
        # 计算并产生asa文件
        asa_folder = os.path.join(dest_folder, '1.1_1_naccess_results')
        if not os.path.exists(asa_folder):
            os.mkdir(asa_folder)
        for file_name in os.listdir(orgnized_pdb_folder):
            pdb2asa(orgnized_pdb_folder, file_name, asa_folder)   

        if self._verbose:
            print(f'Finished checking surface residues, moving result files to {asa_folder}')
            print('Calculating descriptors according to pdb and asa files')
            antibody_complexs = tqdm(antibody_complexs) # 加个进度条
        
        all_descriptors = {}
        
        for ab_complex in antibody_complexs:
            # 根据asa文件获取表面氨基酸
            pdb_id, interested_chains = ab_complex
            file_name = f"{pdb_id}_{'_'.join(interested_chains)}"
            asa_path = os.path.join(asa_folder, f'{file_name}.asa')
            surface_aas = get_surface_aas(asa_path)
        
            # 计算球壳模型部分描述符
            sphere_descriptors = self._sphere_extractor.get_descriptors(interested_chains,
                                                                        surface_aas,
                                                                        orgnized_pdb_folder,
                                                                        file_name)
            # 计算三角部分描述符
            tri_descriptors = self._triangle_extractor.get_descriptors(interested_chains,
                                                                       surface_aas,
                                                                       orgnized_pdb_folder,
                                                                       file_name)
            full_descriptors = self._concat_descriptors(sphere_descriptors, tri_descriptors)
            all_descriptors[file_name] = full_descriptors
            
            if self._verbose:
                descriptor_folder = os.path.join(dest_folder, '3.3_1_Triangle+Spherial')
                if not os.path.exists(descriptor_folder):
                    os.mkdir(descriptor_folder)
                # 输出最终的描述符文件，以csv形式输出会比较好读取
                dest_path = os.path.join(descriptor_folder, f'Descriptors_{file_name}.csv')
                dest_file = open(dest_path, 'w')
                for aa_info, descriptor in full_descriptors.items():
                    aa_type, chain, aa_number = aa_info.split()
                    dest_file.write(f'{chain}_{aa_number}_{aa_type},')
                    descriptor = [str(f'{x:.4f}') for x in descriptor]
                    dest_file.write(f"{','.join(descriptor)}\n")
                dest_file.close()
            
        return all_descriptors
                    
        
        