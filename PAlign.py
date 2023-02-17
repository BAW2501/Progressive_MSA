import numpy as np
import warnings
from functools import partial
from sklearn.cluster import AgglomerativeClustering, ward_tree
from hdbscan import HDBSCAN
from pathlib import Path
from skbio import TreeNode , DNA, RNA, Protein ,io
from skbio.alignment import global_pairwise_align_nucleotide, global_pairwise_align_protein
from tqdm import tqdm



warnings.filterwarnings('ignore',message="You're using skbio's python implementation of Needleman-Wunsch ")
def progressive_msa_func(sequences, pairwise_aligner, guide_tree):
    '''
    Perform progressive multiple sequence alignment using the guide tree.
    '''

    seq_lookup = {s.metadata['id']: s for s in sequences}
    c1, c2 = guide_tree.children
    if c1.is_tip():
        c1_aln = seq_lookup[c1.name]
    else:
        c1_aln = progressive_msa_func(sequences, pairwise_aligner, c1)

    if c2.is_tip():
        c2_aln = seq_lookup[c2.name]
    else:
        c2_aln = progressive_msa_func(sequences, pairwise_aligner, c2)

    alignment, _, _ = pairwise_aligner(c1_aln, c2_aln)

    return alignment


def get_linkage_matrix(children, distances, n_samples):
    '''
    Create linkage matrix (scipy format) from children and distances.
    '''

    # Create linkage matrix and return it
    # create the counts of samples under each node
    counts = np.zeros(children.shape[0])
    for i, merge in enumerate(children):
        counts[i] = sum(1 if child_idx < n_samples else counts[child_idx - n_samples]for child_idx in merge)
    linkage_matrix = np.column_stack([children, distances, counts]).astype(float)
    return linkage_matrix


class DnaMSA:
    def __init__(self, sequences, clustering_algo='Aglo') -> None:
        self.alphabet = list('ACGT')
        self.sequences = sequences
        self.ids = [x.metadata['id'] for x in sequences]
        self.n_sequences = len(sequences)
        self.cluster_algo_to_use = clustering_algo
        self.pair_aligner = partial(global_pairwise_align_nucleotide,match_score = 6,mismatch_score = -2, gap_open_penalty=4, gap_extend_penalty=1)

    def generate_guidetree_aglo(self):
        '''
        Generate guide tree using agglomerative clustering.
        '''
        onehot_sequences_flat = self.one_hot_sequences()
        model = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(onehot_sequences_flat)
        children = model.children_
        distances = model.distances_
        n_samples = len(model.labels_)
        lm = get_linkage_matrix(children, distances, n_samples)
        gt = TreeNode.from_linkage_matrix(lm, self.ids)
        return gt

    def generate_guidetree_ward(self):
        '''
        Generate guide tree using Ward's method.
        '''
        onehot_sequences_flat = self.one_hot_sequences()
        children, _, _, _, distances = ward_tree(onehot_sequences_flat)
        lm = get_linkage_matrix(children, distances, self.n_sequences)
        gt = TreeNode.from_linkage_matrix(lm, self.ids)
        return gt
    
    def generate_guidetree_hdbscan(self):
        '''
        Generate guide tree using HDBSCAN.
        '''
        onehot_sequences_flat = self.one_hot_sequences()
        model = HDBSCAN(min_cluster_size=2, metric='euclidean', cluster_selection_method='eom')
        clusterer = model.fit(onehot_sequences_flat)
        lm = clusterer.single_linkage_tree_.to_numpy()
        gt = TreeNode.from_linkage_matrix(lm, self.ids)
        return gt

    def one_hot_sequences(self):
        max_len = max(len(seq) for seq in self.sequences)
        padded_sequences = np.array(
            [list(str(x) + '_' * (max_len - len(x))) for x in self.sequences])
        onehot_sequences = np.array(
            [
                [1 if x == y else 0 for x in self.alphabet]
                for y in padded_sequences.flatten()
            ]
        )
        return onehot_sequences.reshape((self.n_sequences, max_len * len(self.alphabet)))

    def progressive_msa(self):

        if self.cluster_algo_to_use == 'Aglo':
            gt = self.generate_guidetree_aglo()
        elif self.cluster_algo_to_use == 'Ward':
            gt = self.generate_guidetree_ward()
        elif self.cluster_algo_to_use == 'HDBSCAN':
            gt = self.generate_guidetree_hdbscan()
        else:
            raise ValueError('Invalid clustering algorithm')
        return progressive_msa_func(self.sequences, self.pair_aligner, gt)

#  inherit from DnaMSA to make RnaMSA


class RnaMSA(DnaMSA):
    def __init__(self, sequences, clustering_algo='Aglo') -> None:
        super().__init__(sequences, clustering_algo)
        self.alphabet = list('ACGU')


class ProteinMSA(DnaMSA):
    def __init__(self, sequences, clustering_algo='Aglo') -> None:
        super().__init__(sequences, clustering_algo)
        self.alphabet = list('ACDEFGHIKLMNPQRSTVWY')
        self.pair_aligner = global_pairwise_align_protein

class SequenceFactory:
    '''
    helper class to create MSA objects
    '''
    classes = {'DNA': DnaMSA, 'RNA': RnaMSA, 'Protein': ProteinMSA}
    sequence_types = {'DNA': DNA, 'RNA':RNA, 'Protein':Protein}
    # i should remove some cause they aren't sequence filetypes

    supported_file_types = [
        'gff3' ,
        'fasta' ,
        'qseq' ,
        'fastq' ,
        'embl' ,
        'genbank'
        ]
    def __init__(self,msa_type,clustering_algo):
        self.msa_type = msa_type
        self.clustering_algo = clustering_algo
    
    def init_msa_object_from_strings(self,sequences, ids):
        if self.msa_type not in SequenceFactory.classes:
            raise ValueError('Invalid MSA type')
        class_constructor = SequenceFactory.classes[self.msa_type]
        sequence_constructor = SequenceFactory.sequence_types[self.msa_type]
        seqs = [sequence_constructor(x, metadata={'id': y}) for x, y in zip(sequences, ids)]
        return class_constructor(seqs, self.clustering_algo)
    
    def init_msa_object_from_file(self,file_path,file_type):
        if file_type not in SequenceFactory.supported_file_types:
            raise ValueError('Invalid file type')
        if self.msa_type not in SequenceFactory.classes:
            raise ValueError('Invalid MSA type')
        
        class_constructor = SequenceFactory.classes[self.msa_type]
        sequence_constructor = SequenceFactory.sequence_types[self.msa_type]
        seqs = list(map(sequence_constructor,list(io.read(file_path, format = file_type))))
        return class_constructor(seqs, self.clustering_algo)
            




        


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
    factory = SequenceFactory('DNA','Aglo')
    msa_aligner = factory.init_msa_object_from_strings(X,range(len(X)))
    result = msa_aligner.progressive_msa()
    print(result)
    # RNA test
    rna_X = [seq.replace('T', 'U') for seq in X]
    factory = SequenceFactory('RNA','Aglo')
    msa_aligner = factory.init_msa_object_from_strings(rna_X,range(len(rna_X)))
    result = msa_aligner.progressive_msa()
    print(result)
    # Protein test
    X = [
     'DVVMTQTPLSLPVSLGDQASISCRSSQSLVHSNG',
     'SLPVSLGDQSISCRSSQSLVHSNGNTYLHWYLQKPG',
     'TILDMDVVEGSAARFDCKVEGYPDPE',
     'DVVMQTPLSLPVLGNTYLYWYLQKPG']
    factory = SequenceFactory('Protein','Aglo')
    msa_aligner = factory.init_msa_object_from_strings(X,range(len(X)))
    result = msa_aligner.progressive_msa()
    print(result)
