"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions related to genomics
"""

import os
import json
import pandas as pd
import numpy as np
import scipy


def one_hot_decode_sequence(array):
    """
    Purpose: Given an array [position x 4], decode sequence to a string.
    """
    onehot_decoder = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
    }

    idxs = np.where(array)[1]
    return (''.join([onehot_decoder[i] for i in idxs]))

def one_hot_encode_sequence(sequence):
    """
    Kudos to Charles: /n/projects/cm2363/bpnet-nucleosomes/work/localimportance/allLocalImportances.py
    Purpose: Given a SINGLE sequence string, one-hot encode the data.
        + default control_profiles and control_logcounts is to be zeroed out
        + naively detects whether the sequence is one-hot-encoded.
    """
    onehot_mapping = {
    'A': [1,0,0,0],
    'C': [0,1,0,0],
    'G': [0,0,1,0],
    'T': [0,0,0,1],
    'a': [1,0,0,0],
    'c': [0,1,0,0],
    'g': [0,0,1,0],
    't': [0,0,0,1],
    'N': [0,0,0,0]
    }
    return np.array([onehot_mapping[x] for x in sequence])

def one_hot_encode_sequences(sequences):
    """
    Purpose: Given an array of sequences, one-hot-encode into a [region x position x 4] array.
    """
    return(np.stack([one_hot_encode_sequence(s) for s in sequences]))

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

def swap_motif(seq, positionA, positionB, width):
    """
    Purpose: Given a set of positions and widths, swap the sequences.
    """
    assert (positionA < len(seq)) & (positionB < len(seq)), 'Positions need to be within the sequence margins.'
    assert ((positionA + width) <= positionB), 'PositionA needs to be smaller than positionB and not overlap with the swaps.'
    startA, endA = motif_coords('N'*width, positionA)
    startB, endB = motif_coords('N'*width, positionB)
    new_seq = seq[:startA] + seq[startB:endB] + seq[endA:startB] + seq[startA:endA] + seq[endB:]
    assert len(new_seq) == len(seq)
    return new_seq

def string_to_char_array(seq):
    """
    Kudos: https://alextseng.net/blog/posts/20201122-kmer-shuffles/
    Converts an ASCII string to a NumPy array of byte-long ASCII codes.
    e.g. "ACGT" becomes [65, 67, 71, 84].
    """
    return np.frombuffer(bytes(seq, "utf8"), dtype=np.int8)
        # precis_buckets = np.zeros_like(precis)
        # np.put_along_axis(precis_buckets, thresh_buckets, precis, -1)

def char_array_to_string(arr):
    """
    Kudos: https://alextseng.net/blog/posts/20201122-kmer-shuffles/
    Converts a NumPy array of byte-long ASCII codes into an ASCII string.
    e.g. [65, 67, 71, 84] becomes "ACGT".
    """
    assert arr.dtype == np.int8
    return arr.tostring().decode("ascii")


def one_hot_to_tokens(one_hot):
    """
    Kudos: https://alextseng.net/blog/posts/20201122-kmer-shuffles/
    Converts an L x D one-hot encoding into an L-vector of integers in the range
    [0, D], where the token D is used when the one-hot encoding is all 0. This
    assumes that the one-hot encoding is well-formed, with at most one 1 in each
    column (and 0s elsewhere).
    """
    tokens = np.tile(one_hot.shape[1], one_hot.shape[0])  # Vector of all D
    seq_inds, dim_inds = np.where(one_hot)
    tokens[seq_inds] = dim_inds
    return tokens


def tokens_to_one_hot(tokens, one_hot_dim):
    """
    Kudos: https://alextseng.net/blog/posts/20201122-kmer-shuffles/
    Converts an L-vector of integers in the range [0, D] to an L x D one-hot
    encoding. The value `D` must be provided as `one_hot_dim`. A token of D
    means the one-hot encoding is all 0s.
    """
    identity = np.identity(one_hot_dim + 1)[:, :-1]  # Last row is all 0s
    return identity[tokens]


def shuffle_seqs(seq, num_shufs, k=2, rng=None):
    """
    Kudos: https://alextseng.net/blog/posts/20201122-kmer-shuffles/
    Creates shuffles of the given sequence, in which dinucleotide frequencies
    are preserved.
    Arguments:
        `seq`: either a string of length L, or an L x D NumPy array of one-hot
            encodings
        `num_shufs`: the number of shuffles to create, N
        `k`: the length k-mer whose frequencies are to be preserved; defaults
            to k = 2 (i.e. preserve dinucleotide frequencies)
        `rng`: a NumPy RandomState object, to use for performing shuffles
    If `seq` is a string, returns a list of N strings of length L, each one
    being a shuffled version of `seq`. If `seq` is a 2D NumPy array, then the
    result is an N x L x D NumPy array of shuffled versions of `seq`, also
    one-hot encoded.
    """
    # Convert the sequence (string or one-hot encoded array) into a 1D array of
    # numbers (for simplicity)
    if type(seq) is str:
        arr = string_to_char_array(seq)
    elif type(seq) is np.ndarray and len(seq.shape) == 2:
        seq_len, one_hot_dim = seq.shape
        arr = one_hot_to_tokens(seq)
    else:
        raise ValueError("Expected string or one-hot encoded array")

    if not rng:
        rng = np.random.RandomState()

    if k == 1:
        # Do simple shuffles of `arr`
        if type(seq) is str:
            all_results = []
        else:
            all_results = np.empty(
                (num_shufs, seq_len, one_hot_dim), dtype=seq.dtype
            )
        for i in range(num_shufs):
            rng.shuffle(arr)
            if type(seq) is str:
                all_results.append(char_array_to_string(arr))
            else:
                all_results[i] = tokens_to_one_hot(arr, one_hot_dim)
        return all_results

    # Tile `arr` from a 1D array to a 2D array of all (k-1)-mers (i.e.
    # "shortmers"), using -1 as a "sentinel" for the last few values
    arr_shortmers = np.empty((len(arr), k - 1), dtype=arr.dtype)
    arr_shortmers[:] = -1
    for i in range(k - 1):
        arr_shortmers[:len(arr) - i, i] = arr[i:]

    # Get the set of all shortmers, and a mapping of which positions start with
    # which shortmers; `tokens` is the mapping, and is an integer representation
    # of the original shortmers (range [0, # unique shortmers - 1])
    shortmers, tokens = np.unique(arr_shortmers, return_inverse=True, axis=0)

    # For each token, get a list of indices of all the tokens that come after it
    shuf_next_inds = []
    for token in range(len(shortmers)):
        # Locations in `arr` where the shortmer exists; some shortmers will have
        # the sentinel, but that's okay
        mask = tokens == token
        inds = np.where(mask)[0]
        shuf_next_inds.append(inds + 1)  # Add 1 to indices for next token

    if type(seq) is str:
        all_results = []
    else:
        all_results = np.empty(
            (num_shufs, seq_len, one_hot_dim), dtype=seq.dtype
        )

    for i in range(num_shufs):
        # Shuffle the next indices
        for t in range(len(shortmers)):
            inds = np.arange(len(shuf_next_inds[t]))
            inds[:-1] = rng.permutation(len(inds) - 1)  # Keep last index same
            shuf_next_inds[t] = shuf_next_inds[t][inds]

        counters = [0] * len(shortmers)

        # Build the resulting array
        ind = 0
        result = np.empty_like(tokens)
        result[0] = tokens[ind]
        for j in range(1, len(tokens)):
            t = tokens[ind]
            ind = shuf_next_inds[t][counters[t]]
            counters[t] += 1
            result[j] = tokens[ind]

        shuffled_arr = shortmers[result][:, 0]  # First character of each shortmer
        # (this leaves behind the sentinels)

        if type(seq) is str:
            all_results.append(char_array_to_string(shuffled_arr))
        else:
            all_results[i] = tokens_to_one_hot(shuffled_arr, one_hot_dim)
    return all_results

def logitsToProfile(logitsAcrossSingleRegion, logCountsAcrossSingleRegion):
    """
    Purpose: Given a single task and region sequence prediction (position x channels),
        convert output logits/logcounts to human-readable representation of profile prediction.
    """
    assert len(logitsAcrossSingleRegion.shape)==2 #Logits will have shape (output-width x numTasks)
    assert len(logCountsAcrossSingleRegion.shape)==1 #Logits will be a scalar value

    profileProb = scipy.special.softmax(logitsAcrossSingleRegion)
    profile = profileProb * np.exp(logCountsAcrossSingleRegion)
    return profile
