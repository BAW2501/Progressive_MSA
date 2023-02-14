import numpy as np
from sklearn.cluster import AgglomerativeClustering, ward_tree
from skbio import TreeNode
from skbio.alignment import global_pairwise_align_nucleotide, global_pairwise_align_protein


def progressive_msa_func(sequences, pairwise_aligner, guide_tree):

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
    # sourcery skip: inline-immediately-returned-variable
    # Create linkage matrix and then plot the dendrogram
    # create the counts of samples under each node
    counts = np.zeros(children.shape[0])
    for i, merge in enumerate(children):
        counts[i] = sum(1 if child_idx < n_samples else counts[child_idx - n_samples]
                        for child_idx in merge)
    linkage_matrix = np.column_stack(
        [children, distances, counts]).astype(float)
    return linkage_matrix


class DnaMSA:
    def __init__(self, sequences, clustering_algo='Aglo') -> None:
        self.alphabet = list('ACGT')
        self.sequences = sequences
        self.ids = [x.metadata['id'] for x in sequences]
        self.n_sequences = len(sequences)
        self.cluster_algo_to_use = clustering_algo
        self.pair_aligner = global_pairwise_align_nucleotide

    def generate_guidetree_aglo(self):
        onehot_sequences_flat = self.one_hot_sequences()
        model = AgglomerativeClustering(
            distance_threshold=0, n_clusters=None).fit(onehot_sequences_flat)
        children = model.children_
        distances = model.distances_
        n_samples = len(model.labels_)
        return get_linkage_matrix(children, distances, n_samples)

    def generate_guidtree_ward(self):

        onehot_sequences_flat = self.one_hot_sequences()
        children, _, _, _, distances = ward_tree(onehot_sequences_flat)
        return get_linkage_matrix(children, distances, self.n_sequences)

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
            linkage_matrix = self.generate_guidetree_aglo()
            gt = TreeNode.from_linkage_matrix(linkage_matrix, self.ids)
        elif self.cluster_algo_to_use == 'Ward':
            linkage_matrix = self.generate_guidtree_ward()
            gt = TreeNode.from_linkage_matrix(linkage_matrix, self.ids)
        else:
            raise ValueError('Invalid clustering algorithm')
        return progressive_msa_func(self.sequences, self.pair_aligner, gt)

#  inherit from DnaMSA to make RnaMSA


class RnaMSA(DnaMSA):
    def __init__(self, sequences, clustering_algo='Aglo') -> None:
        super().__init__(sequences, clustering_algo)
        self.alphabet = list('ACGU')


# inherit from dnaMSA to make ProteinMSA by changing alphabet and modifying progressive_msa method, to use global_pairwise_align_protein

class ProteinMSA(DnaMSA):
    def __init__(self, sequences, clustering_algo='Aglo') -> None:
        super().__init__(sequences, clustering_algo)
        self.alphabet = list('ACDEFGHIKLMNPQRSTVWY')
        self.pair_aligner = global_pairwise_align_protein
