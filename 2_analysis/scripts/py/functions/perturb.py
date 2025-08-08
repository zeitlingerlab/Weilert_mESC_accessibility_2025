"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to applying perturbations
"""

import os
import sys
import json
import pandas as pd
import numpy as np
import re

#Custom functions
sys.path.insert(0, f'/n/projects/mw2098/shared_code/bpreveal/functions')
from functional import one_hot_decode_sequence, one_hot_encode_sequence


def generate_random_seq(seqlen, weights = [.25, .25, .25, .25]):
    """
    Purpose: Generate a random DNA sequence of a specified length.
    """
    import random
    seq = random.choices(['A','C','G','T'], weights = weights, k=seqlen)
    return(''.join(seq))

def motif_coords(motif, position):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L19
    Purpose: Given motif (string) and a center position, find the motif boundaries.
    """
    start = position - len(motif) // 2
    end = start + len(motif)
    # print(start, end)
    return start, end

def insert_motif(seq, motif, position):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L25
    Purpose: Given a sequence, inject a motif centered on a position.
    """
    assert position < len(seq)
    start, end = motif_coords(motif, position)
    new_seq = seq[:start] + motif + seq[end:]
    assert len(new_seq) == len(seq)
    return new_seq

def generate_injected_seq(primary_motif, primary_position,
                          secondary_motif='', secondary_distances=[], seqlen=1000,
                          weights = [.25, .25, .25, .25]):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L37
    Purpose: Given a sequence, inject a motif centered on a position.
        + If you want to inject 2 motifs to compare distances, use:
            + side_motif to define the second motif sequence
            + side_distances to define the second motif center
    """
    random_seq = generate_random_seq(seqlen = seqlen, weights = weights)
    injected_seq = insert_motif(seq = random_seq, motif = primary_motif, position = primary_position)
    if len(secondary_distances)>1:
        print('Warning! You have entered multiple side distances.',
              'This will inject multiple secondary motifs into the SAME sequence.',
              'If you do not want this, loop this function through single side distances.')
    for d in secondary_distances:
        injected_seq = insert_motif(injected_seq, secondary_motif, d)
    return injected_seq

def generate_alt_sequences(ref_seq, motifs, motif_unique_col, sequence_window_start_col, sequence_window_end_col,
                           comb_max = None, comb_min = None):
    """
    Purpose: Generate all combinations of mutated sequences based on given sequence and desired mutant coordinates.
    Given:
        ref_seq: np.array with shape [l x 4] of the one-hot encoded reference sequence
        motifs: pd.df of the motifs across the reference sequence
        motif_unique_col: column name designating the unique motif label
        sequence_window_start_col: column name designating where in the input sequence the motifs start
        sequence_window_end_col: column name designating where in the input sequence the motifs end
        comb_max: max number of simultaneous mutations allowed
        comb_min: min number of simultaneous mutations allowed
    Output:
        perturb_seqs_dict_list = trials x -> ref/mutation -> seqlen x 4
    """
    import sys
    from itertools import combinations, chain


    #Insert relevant mutant sequences in each combination
    def get_alt_seq_from_combo(combo, ref_seq, motifs, motif_unique_col, sequence_window_start_col, sequence_window_end_col, mut_seqs_dict):
        alt_seq = ref_seq.copy()
        for mut in combo:
            motif = motifs[motifs[motif_unique_col]==mut]
            alt_seq[motif[sequence_window_start_col].iloc[0]:motif[sequence_window_end_col].iloc[0]] = mut_seqs_dict[mut]
        return alt_seq

    #Mark the max depth of combinations
    c_min_depth = 1
    c_max_depth = motifs.shape[0]
    if(comb_max is not None): c_max_depth = comb_max
    if(comb_min is not None): c_min_depth = comb_min

    #Ensure these coordinate positions are integers
    motifs[sequence_window_start_col] = motifs[sequence_window_start_col].astype(int)
    motifs[sequence_window_end_col] = motifs[sequence_window_end_col].astype(int)

    #Determine all combinations of unique motifs
    mut_combos = [list(combinations(motifs[motif_unique_col], d)) for d in range(c_min_depth, c_max_depth+1)]
    mut_combos = list(chain(*mut_combos))

    #Record mutant sequences used in this combination: motif -> motif_len x 4
    mut_seqs_dict = {motifs[motif_unique_col].iloc[idx]:
                          one_hot_encode_sequence(generate_random_seq(motifs[sequence_window_end_col].iloc[idx] - motifs[sequence_window_start_col].iloc[idx]))
                     for idx in range(len(motifs))}

    # Generate: combo -> 1000 x 4
    alt_seqs_dict = {'_'.join(list(combo)): get_alt_seq_from_combo(combo = combo,
                                                                   ref_seq = ref_seq,
                                                                   motifs = motifs,
                                                                   motif_unique_col = motif_unique_col,
                                                                   sequence_window_start_col = sequence_window_start_col,
                                                                   sequence_window_end_col = sequence_window_end_col,
                                                                   mut_seqs_dict = mut_seqs_dict)
                     for combo in mut_combos}

    #Add reference sequence
    perturb_seqs_dict = {**{'Reference': ref_seq}, **alt_seqs_dict}

    #Clean up results
    return(perturb_seqs_dict)
