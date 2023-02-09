from progressive_align import *

# just an example guide tree for testing

def example_guide_tree_generation(sequences, metric=hamming_distance):
        guide_dm = DistanceMatrix.from_iterable(sequences, metric=metric, key='id') # maybe NN here to get the distance matrix 
        # print(guide_dm) # uncomment to see the distance matrix
        guide_lm = sp.cluster.hierarchy.average(guide_dm.condensed_form()) # kmeans or some other hierarchical clustering algorithm
        guide_tree = TreeNode.from_linkage_matrix(guide_lm, guide_dm.ids) # i'm figuring out how to get the tree from the linkage matrix
        return guide_tree

if "__main__" == __name__:
    # sequences from our TP
    query_of_sequences = [
                            DNA("ACGCGATGACCGGGCCTTGTAAAAT", metadata={'id': 'seq1'}),
                            DNA("ATGATGACAGGGCTTGTAACT", metadata={'id': 'seq2'}),
                            DNA("TTCATGACCGGCTTATACTTAT", metadata={'id': 'seq3'}),
                            DNA("ACCCTACCTGTCGTATTGTAAT", metadata={'id': 'seq4'}),
                        ]
    guide_tree1 = example_guide_tree_generation(query_of_sequences,metric=kmer_distance)
    guide_tree2 = example_guide_tree_generation(query_of_sequences,metric=hamming_distance)
    guide_tree3 = example_guide_tree_generation(query_of_sequences,metric=minimum_edit_distance)
    guide_tree4 = example_guide_tree_generation(query_of_sequences,metric=alignement_score_distance)
    # global_pairwise_align_nucleotide this function is faster than our TP implementation and works with profiles
    # but it's not nearly as fast as we need it to be( dies at length 1000+ sequences)
    # TODO: supress the warnings
    # TODO: Faster Aligner ( preferably with C and bindings to python )
    aligner = partial(global_pairwise_align_nucleotide,match_score = 6,mismatch_score = -2, gap_open_penalty=4, gap_extend_penalty=1)
    # you can use any of the guide trees above
    MSA = progressive_msa(query_of_sequences, pairwise_aligner=aligner,guide_tree=guide_tree4)
    print(MSA)
    MSA.write("MSA.fasta",format="fasta")
