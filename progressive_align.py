from functools import partial
import scipy as sp
import numpy as np
from skbio import DNA, TabularMSA, TreeNode, DistanceMatrix, Sequence
from skbio.alignment import global_pairwise_align_nucleotide


def kmer_distance(sequence1, sequence2, k=7, overlap=False):
    sequence1_kmers = set(map(str, sequence1.iter_kmers(k, overlap)))
    sequence2_kmers = set(map(str, sequence2.iter_kmers(k, overlap)))
    all_kmers = sequence1_kmers | sequence2_kmers
    shared_kmers = sequence1_kmers & sequence2_kmers
    number_unique = len(all_kmers) - len(shared_kmers)
    fraction_unique = number_unique / len(all_kmers)
    return fraction_unique


def hamming_distance(sequence1, sequence2):
    n, m = len(sequence1), len(sequence2)
    return sum(sequence1[i] != sequence2[i] for i in range(min(n, m))) + abs(n-m)


def minimum_edit_distance(sequence1, sequence2):
    # dynamic programming
    n, m = len(sequence1), len(sequence2)
    dp = np.zeros((n+1, m+1))
    for i in range(n+1):
        dp[i][0] = i
    for j in range(m+1):
        dp[0][j] = j
    for i in range(1, n+1):
        for j in range(1, m+1):
            dp[i][j] = min(dp[i-1][j]+1, dp[i][j-1]+1, dp[i-1]
                           [j-1]+(sequence1[i-1] != sequence2[j-1]))
    return dp[n][m]


def alignement_score_distance(sequence1, sequence2):
    if sequence1 == sequence2:
        return 0
    _, score, _ = global_pairwise_align_nucleotide(sequence1, sequence2)
    return score


def example_guide_tree_generation(sequences, metric=hamming_distance):
    guide_dm = DistanceMatrix.from_iterable(sequences, metric=metric, key='id')
    print(guide_dm)
    guide_lm = sp.cluster.hierarchy.average(guide_dm.condensed_form())
    guide_tree = TreeNode.from_linkage_matrix(guide_lm, guide_dm.ids)
    return guide_tree


def progressive_msa(sequences, pairwise_aligner, guide_tree=None, metric=kmer_distance):

    if guide_tree is None:
        guide_dm = DistanceMatrix.from_iterable(
            sequences, metric=metric, key='id')
        guide_lm = sp.cluster.hierarchy.average(guide_dm.condensed_form())
        guide_tree = TreeNode.from_linkage_matrix(guide_lm, guide_dm.ids)

    seq_lookup = {s.metadata['id']: s for i, s in enumerate(sequences)}
    c1, c2 = guide_tree.children
    c1_aln = seq_lookup[c1.name] if c1.is_tip() else progressive_msa(
        sequences, pairwise_aligner, c1)
    c2_aln = seq_lookup[c2.name] if c2.is_tip() else progressive_msa(
        sequences, pairwise_aligner, c2)

    alignment, _, _ = pairwise_aligner(c1_aln, c2_aln)

    return alignment
