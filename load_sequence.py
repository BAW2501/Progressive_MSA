import sys
import os
import functools
import collections

from matplotlib import pyplot as plt

from skbio import TabularMSA
import fastaparser
from sklearn import manifold
from sklearn.manifold import MDS
from skbio import Protein
from skbio.sequence.distance import kmer_distance
from skbio import DistanceMatrix

def read_sequence_file_fasta():
    fn = sys.argv[1]
    if os.path.exists(fn):
        with open(fn) as fasta_file:
            sequences_objects = fastaparser.Reader(fasta_file)

            # sequences
            sequences = collections.OrderedDict()

            for seq in sequences_objects:
                sequences[seq.id] = Protein(sequence=seq.sequence_as_string(), metadata={'id': seq.id})
        
    return sequences

def embed_sequences():
    sequences=read_sequence_file_fasta()
    #build distace marix
    kmer_distances = DistanceMatrix.from_iterable(sequences.values(), metric=functools.partial(kmer_distance, k=8),
                                                  key='id')
     # init mds model
    mds= manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, dissimilarity="precomputed", n_jobs=1,
        normalized_stress="auto",
    )

    pos = mds.fit(kmer_distances.data).embedding_
    s = 100
    plt.scatter(pos[:, 0], pos[:, 1], color="turquoise", s=s, lw=0, label="MDS")
    plt.show()
    print(pos)


embed_sequences()