import os
import sys
import warnings

from skbio.util import SkbioWarning
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from tqdm import tqdm
import skbio.io

from PAlign import *
# read file that contains the filenames for the familes
def read_ids():
    # first argument is the id file that contains all the file name of the dataset (one of OXbench datasets)
    fn = sys.argv[1]
    # each file specifies sequences of a family which will be msad together
    families=[]

    ox_id_file = open(fn, 'r')
    Lines = ox_id_file.readlines()
    for line in Lines:
        families.append(line.strip())
        
    ox_id_file.close()
    return families

#read the fasta files for each family
def read_families(families, sequence_constructor=None):
    # second argument is the directory that contains all the original familes, each familes has a set of sequences
    directory= sys.argv[2]
    for family in tqdm(families):
        sequences= list(skbio.read(directory+'\\'+family+'.fa', format='fasta'))
        sequences= list(map(sequence_constructor, sequences))
        for i, seq in enumerate(sequences):
            seq.metadata['id'] = i
            msa_aligner = ProteinMSA(sequences)
        result = msa_aligner.progressive_msa(preprocess='one_hot')
        result.write('output\\'+family+'.fa', format='fasta')


def align_families():
    families= read_ids()
    read_families(families, Protein)
    print()

align_families()
