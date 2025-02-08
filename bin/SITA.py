import os
import argparse
from colorama import init, Fore

from feature_extractor import FeatureExtractor
from model import SITA


# CLI tool for SITA
parser = argparse.ArgumentParser()
parser.add_argument(
    '-p', '--pdb', help='path to your PDB file folder', required=True)
parser.add_argument(
    '-c', '--chain', help='path to your chain file', required=True)
parser.add_argument(
    '-o', '--out', help='path to your output folder', required=True)
parser.add_argument(
    '-v', '--verbose', help='provides more detailed output if set true', action='store_true')

args = parser.parse_args()

verbose_mode = args.verbose # 是否输出详细过程
pdb_folder = args.pdb
chain_file_path = args.chain
dest_folder = args.out

if not (os.path.isdir(dest_folder) and os.path.isdir(pdb_folder) 
        and os.path.isfile(chain_file_path)):
    init()
    print(f'{Fore.RED} Wrong input format') 
    print(f'{Fore.RED} -p, -o and -c should all be provided')
    print(f'{Fore.RED} -p and -o should be an existing folder and -c should be a file')
    exit(0)    
    
print("")
print(" +------------------------------------+")
print(" |  /----\   --+--  --------   /\     |")
print(" | |           |       |      /  \    |")
print(" |  \----\     |       |     /----\   |")
print(" |        |    |       |    /      \  |")
print(" |  \----/   __|__     |   /        \ |")
print(" +------------------------------------+")
print("")

if verbose_mode: 
    print("running in verbose mode")

extractor = FeatureExtractor(verbose_mode, dest_folder)
descriptors = extractor.get_features(pdb_folder, chain_file_path, dest_folder)

print("FINISHED EXTRACTING FEATURES, START MAKING PREDICTIONS")
dest_folder = os.path.join(dest_folder, 'SITA_prediction')
if not os.path.isdir(dest_folder):
    os.mkdir(dest_folder)

model = SITA(verbose_mode)
model.mk_prediction(descriptors, dest_folder)

print('finished predicting, BYE~')
exit(0)

