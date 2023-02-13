import os
import sys

from functools import partial
from skbio import TabularMSA, Protein
from tqdm import tqdm


# read file that contains the filenames for the familes
def read_ids():
    # first argument is the id file that contains all the file name of the dataset (one of OXbench datasets)
    fn = sys.argv[1]
    if os.path.exists(fn):

        # each file specifies sequences of a family which will be msad together    
        families=[]

        ox_id_file = open(fn, 'r')
        Lines = ox_id_file.readlines()
        for line in Lines:
            families.append(line.strip())
        
        ox_id_file.close()
    
    return families

# replace this stub with accuracy measures
def stub_func(string):
    print(string)

def read_alignments(families, func):
    #read the fasta files for each family
    directory= sys.argv[2]
    print(families)
    for family in tqdm(families):
        MSA = TabularMSA.read(directory+'\\'+family, constructor=partial(Protein, lowercase=True))
        MSA.reassign_index(minter='id')
        
        positional_conservation = MSA.conservation(metric='inverse_shannon_uncertainty', degenerate_mode='nan', gap_mode='include')
        print(positional_conservation)
        func(MSA)

read_alignments(read_ids(),stub_func)