from skbio import TabularMSA, Protein
from itertools import combinations
from functools import partial
from tqdm import tqdm
import pandas as pd
import blosum as bl
import sys


# get id names from file
def read_ids(_full_ids: str) -> list[str]:
    '''
    Read the IDs from the file
    :param _full_ids: the path to the file containing the IDs
    :return: a list of IDs'''
    with open(_full_ids, 'r') as f:
        # Read the lines of the file and strip any whitespace
        # Return the list of lines as the result
        return [line.strip() for line in f.readlines()]
    
def get_msa_from_file(id,main_folder):
    '''
    Get the multiple sequence alignment from the file by the ID
    :param id: the ID of the alignment
    :param main_folder: the folder containing the alignments
    :return: the multiple sequence alignment
    '''
    return TabularMSA.read(f'{main_folder}/{id}.fa', constructor=partial(Protein, lowercase=True), format='fasta')


def calc_align_score(seq1, seq2, score_mat):
    '''
    Calculate the score alignment of two sequences
    :param seq1: first sequence
    :param seq2: second sequence
    :param score_mat: the substitution matrix to use usually BLOSUM or PAM
    :return: the score of the alignment
    '''
    score = 0
    for base1,base2 in list(zip(seq1, seq2)):
        base1 = str(base1).replace('-', '*')
        base2 = str(base2).replace('-', '*')
        score += -1 if (base1 == base2 == '-') else score_mat[base1+base2]
    return score

def calc_SP_score(msa,scoring_mat):
    '''
    Calculate the SP score of a multiple sequence alignment
    :param msa: the multiple sequence alignment
    :param scoring_mat: the substitution matrix to use usually BLOSUM or PAM
    :return: the SP score of the alignment
    '''
    return sum(calc_align_score(seq1, seq2, scoring_mat) for seq1,seq2 in combinations(msa,2))

def calc_CS_score(msa,scoring_mat):
    '''
    Calculate the circular sum score of a multiple sequence alignment
    :param msa: the multiple sequence alignment
    :param scoring_mat: the substitution matrix to use usually BLOSUM or PAM
    :return: the CS score of the alignment
    '''
    # sum score of i and i+1 for all i in range(len(msa)) also add the score of the -1 and 0 sequence
    # -1,0 0,1 1,2 ..... n-1,n 0,n-1
    return sum(calc_align_score(msa[i], msa[i+1], scoring_mat) for i in range(-1,len(msa)))

# TODO calc_column_score
# TODO calc_AQ_score

def calc_AQ_score(msa,scoring_mat):
    ...

def calc_column_score(column,scoring_mat):
    ...


if __name__ == '__main__':
    full_ids = r'C:\Users\Dell\Desktop\oxbench_1_3\data\full.id'
    full_folder = r'C:\Users\Dell\Desktop\oxbench_1_3\data\align\fasta'
    ids = read_ids(full_ids)

    # get the number of sequences in the alignment
    # get the SP score of the alignment
    # get the CS score of the alignment
    # add the data to the dataframe
    data = pd.DataFrame(columns=['ID','NB_sequences','SP','CS'])

    for id in tqdm(ids):
        msa_fam = get_msa_from_file(id,full_folder)
        row_data = {
            'ID': id,
            'NB_sequences': len(msa_fam),
            'SP': calc_SP_score(msa_fam, bl.BLOSUM(62)),
            'CS': calc_CS_score(msa_fam, bl.BLOSUM(62))
        }
        data = data.append(row_data, ignore_index=True)

    data.to_csv('benchmark_SP_CS.csv', index=False)
