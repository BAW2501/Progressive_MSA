import skbio
import os
import sys

from skbio import Protein
from skbio.alignment import global_pairwise_align_protein

sequences= list(skbio.read(sys.argv[1], format='fasta'))
sequences= list(map(Protein, sequences))

alignment = global_pairwise_align_protein(sequences[0], sequences[1], penalize_terminal_gaps=True, gap_open_penalty=int(sys.argv[3]), gap_extend_penalty=int(sys.argv[4]))

alignment[0].write(sys.argv[2], format='fasta')


