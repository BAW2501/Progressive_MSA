from functools import partial
import scipy as sp
from skbio import DistanceMatrix
from skbio import DNA


def kmer_distance(sequence1, sequence2, k=7, overlap=False):
    sequence1_kmers = set(map(str, sequence1.iter_kmers(k, overlap)))
    sequence2_kmers = set(map(str, sequence2.iter_kmers(k, overlap)))
    all_kmers = sequence1_kmers | sequence2_kmers
    shared_kmers = sequence1_kmers & sequence2_kmers
    number_unique = len(all_kmers) - len(shared_kmers)
    fraction_unique = number_unique / len(all_kmers)
    return fraction_unique


def guide_tree_from_sequences(sequences, metric=kmer_distance, display_tree=False):
    guide_dm = DistanceMatrix.from_iterable(sequences, metric=metric, key='id')
    guide_lm = sp.cluster.hierarchy.average(guide_dm.condensed_form())
    guide_tree = sp.cluster.hierarchy.to_tree(guide_lm)
    if display_tree:
        sp.cluster.hierarchy.dendrogram(
            guide_lm, labels=guide_dm.ids, orientation='right', link_color_func=lambda x: 'black')
    return guide_tree


filenames = [
    "APOE Variant 1.txt",
    "APOE Variant 2.txt",
    "APOE Variant 3.txt",
    "APOE Variant 4.txt",
    "Cloning Vector APOE2.txt",
    "Cloning Vector APOE3.txt",
    "ABCA7 - Brown Rat.txt",
    "ABCA7 - House Mouse V1.txt",
    "ABCA7 House Mouse V2.txt",
    "ABCA7 Human.txt"]

texts = []
for filename in filenames:
    with open(f"Sequences/{filename}", encoding="utf8") as myfile:
        texts.append(myfile.read().replace('\n', ''))

v3, txtL, txtBRat, txtHV1, txtHV2, txtHum, ratCut, m1Cut, m2Cut, hCut = texts

s1 = DNA("ACCGGTGACCAGTTGACCAGT")
s2 = DNA("ATCGGTACCGGTAGAAGT")
s3 = DNA("GGTACCAAATAGAA")
BrownRat = DNA(txtBRat)
HouseM1 = DNA(txtHV1)
HouseM2 = DNA(txtHV2)
HumanA = DNA(txtHum)
ratCut = DNA(ratCut)
m1Cut = DNA(m1Cut)
m2Cut = DNA(m2Cut)
hCut = DNA(hCut)
print(kmer_distance(HouseM2, HouseM1, k=7, overlap=False))
print(HumanA.distance(HouseM1, kmer_distance))
print(HouseM1.distance(HouseM2, kmer_distance))

fivemer_distance = partial(kmer_distance, k=5)

s1 = DNA("ACCGGTGACCAGTTGACCAGT")
s2 = DNA("ATCGGTACCGGTAGAAGT")
s3 = DNA("GGTACCAAATAGAA")

print(s1.distance(s2, fivemer_distance))
print(s1.distance(s3, fivemer_distance))

query_of_sequences = [DNA(BrownRat, {"id": "s1"}),
                      DNA(HouseM1, {"id": "s2"}),
                      DNA(HouseM2, {"id": "s3"}),
                      DNA(HumanA, {"id": "s4"})
                      ]

t = guide_tree_from_sequences(query_of_sequences, display_tree=True)
