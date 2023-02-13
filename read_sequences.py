import os
import sys
import fastaparser
import warnings

from skbio.util import SkbioWarning
from progressive_align import *
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from tqdm import tqdm


# init function for MSA
def example_guide_tree_generation(sequences, metric=hamming_distance):
    guide_dm = DistanceMatrix.from_iterable(sequences, metric=metric,
                                            key='id')  # maybe NN here to get the distance matrix
    # print(guide_dm) # uncomment to see the distance matrix
    guide_lm = sp.cluster.hierarchy.average(
        guide_dm.condensed_form())  # kmeans or some other hierarchical clustering algorithm
    guide_tree = TreeNode.from_linkage_matrix(guide_lm,
                                              guide_dm.ids)  # i'm figuring out how to get the tree from the linkage matrix
    return guide_tree

aligner = partial(global_pairwise_align_protein, gap_open_penalty=4, gap_extend_penalty=1)


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

#read the fasta files for each family
def read_families(families):
    # second argument is the directory that contains all the original familes, each familes has a set of sequences
    directory= sys.argv[2]
    for family in tqdm(families):
        with open(directory+"\\"+family+".fa") as fasta_file:
            family_sequences = fastaparser.Reader(fasta_file)
            family_sequences_aligned=[]
            for seq in family_sequences:
            # seq is a FastaSequence object
                family_sequences_aligned.append( Protein(seq.sequence_as_string(), metadata={'id':seq.id, 'Description': seq.description}) )
                #print('ID:', seq.id)
                #print('Description:', seq.description)
                #print('Sequence:', seq.sequence_as_string())
                #print()
            guide_tree= example_guide_tree_generation(family_sequences_aligned, metric=hamming_distance)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=SkbioWarning)
                MSA = progressive_msa(family_sequences_aligned, pairwise_aligner=aligner, guide_tree=guide_tree)
            MSA.write(family+'.fa', format='fasta')
            #print('family : ', family)
            #print(MSA)
            #print(MSA[0])
            #print(MSA[1])
            #print('//////////////////////////////////////////////////////////////')
    list_files(directory)
    print()



def align_families():
    families= read_ids()
    read_families(families)
    print()

align_families()


def read_fasta_from_id_file():
    print()



