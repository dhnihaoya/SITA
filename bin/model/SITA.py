import os
import numpy as np
import xgboost as xgb
from tqdm import tqdm
from typing import Dict, List

class SITA:
    _crt_path = os.path.dirname(__file__)
    
    def __init__(self, verbose: bool):
        model_folder = os.path.join(SITA._crt_path, 'model_files')
        self._models = []
        # 加载所有模型文件
        for file_name in os.listdir(model_folder):
            # 正好避免macos里的.DStore
            if not file_name.endswith('json'): 
                continue
            model_path = os.path.join(model_folder, file_name)
            model = xgb.XGBClassifier()
            model.load_model(model_path)
            self._models.append(model)
        self._verbose = verbose
    
    def mk_prediction(self, pdb_features: Dict[str, Dict[str, List[float]]], 
                      dest_folder: str):

        for pdb_info, aa_feas in tqdm(pdb_features.items()):

            # 拼起来所有的特征,用模型集成一次性算完，然后全输出到文件中
            aa_names = []
            all_feas = []
            for aa_info, fea in aa_feas.items():
                aa_type, chain, aa_number = aa_info.split()
                # aa_names.append(f'{chain}_{aa_number}_{aa_type}')
                # 只保留位数==781的情况（没CA会出809位）
                if len(fea) == 781: 
                    all_feas.append(fea)
                    aa_names.append(f'{chain}_{aa_number}_{aa_type}')
            aa_feas = np.array(all_feas)
            # 拼接完之后用各个模型依次预测，
            esemble_results = [
                model.predict_proba(all_feas)[:, 1]
                for model in self._models
            ]
            result = np.mean(esemble_results, axis=0)                        
            
            dest_path = os.path.join(dest_folder, f'{pdb_info}.txt')
            dest_file = open(dest_path, 'w')
            for (aa_name, score) in zip(aa_names, result):
                dest_file.write(f'{aa_name}\t{score:.6f}\n')
            dest_file.close()
        
        

        