import os
import re
import sys
import itertools

from functools import partial
import blosum as bl
from skbio import TabularMSA, Protein
from tqdm import tqdm

from PAlign import *

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

def compute_score_pairwise_alignment(seq1, seq2, score_matrix_dict):
    score=0
    for residue_pair in list(zip(seq1, seq2)):
        residue1=residue_pair[0].__str__()
        residue2=residue_pair[1].__str__()

        # blosum matrices in the blosum modules read gaps as '*' instead of '-'
        if(residue1=='-'):
            residue1='*'

        if(residue2=='-'):
            residue2='*'

        # add cost of substuting residue 1 with 2
        score+=score_matrix_dict[residue1+residue2]
    return score

def sum_pairs_score(msa, scoringMatrix):
    score=0
    #get number of blossum matrix
    match= re.search("blossum(.*)", scoringMatrix).group(1)
    score_matrix_dict=bl.BLOSUM(int(match))
    pairwise_combinations = itertools.combinations(msa, 2)
    
    for pair in pairwise_combinations:
        score+= compute_score_pairwise_alignment(pair[0], pair[1], score_matrix_dict)

    return score

## TODO test sum_pairs_score function
## TODO implement circular sum
## TODO implement column score


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



#read_alignments(read_ids(),stub_func)

if __name__ == '__main__':
    # DNA test
    X = [
    'ACGCGATGACCGGGCCTTGTATT',
    'ATGATGACAGGGCTTGTAACT',
    'TTCATGACCGGCTTATACTTA',
    'ACGCGATGACCGGGCCTTGTA',
    'ATGATGACAGGGCTTGTAACT',
    'TTCATGACCGGCTTATACTTA',
    'CGTTGCCGTTACATTTAAGGC',
    'GCAGGCGTTTACTTACGCAGT',
    'TCGTTAGCCTGGTTTTTACCC',
    'AACGCAGGTTAGTGGTACGTT',
    'CTGCCGTGTATTTCACAGGAC',
    'GAGCGGCTCACGGGTTTGGTT',
    'TTAGCCGAGTTTTAGTTTCCG'
    ]

    ## small Protein test
    Y=[
        'ATLKEKLIAPVAQQETTIPDNKITVVGVGQVGMACAISILGKSLTDELALVDVLEDKLKGEMMDLQHGSLFLQTPKIVANKDYSVTANSKIVVVTAGVRQQEGESRLNLVQRNVNVFKFIIPQIVKYSPNCIIIVVSNPVDILTYVAWKLSGLPKHRVIG',
        'TALKDKLIGHLATSQEPRSYNKITVVGCDAVGMADAISVLMKDLADEVALVDVMEDKLKGEMMDLEHGSLFLHTAKIVSGKDYSVSAGSKLVVITAGARQQEGESRLNLVQRNVNIFKFIIPNIVKHSPDCLKELHPELGTDKNKQDWKLSGLPMHRIIGSG',
        'ATLKDKLIGHLATSQEPRSYNKITVVGVGAVGMACAISILMKDLADEVALVDVMEDKLKGEMMDLQHGSLFLHTAKIVSGKDYSVSAGSKLVVITAGARQQEGESRLNLVQRNVNIFKFIIPNIVKHSPDCIILVVSNPVDVLTYVAWKLSGLPMHRIIGSG',
    ]
    sequences = list(map(Protein,Y))

    for i,seq in enumerate(sequences):
        seq.metadata['id'] = i

    msa_aligner = ProteinMSA(sequences)
    result = msa_aligner.progressive_msa()
    print(sum_pairs_score(result, 'blossum62'))

