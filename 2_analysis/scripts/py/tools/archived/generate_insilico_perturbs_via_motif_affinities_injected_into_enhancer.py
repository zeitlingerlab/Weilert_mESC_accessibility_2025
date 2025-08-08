
"""
Melanie Weilert
May 2023
Purpose: Given an enhancer coordinate, generate injections into a motif group.
"""

##################################################################
# Computational setup
##################################################################
#Packages
import os
import sys
import keras
import json
import pandas as pd
from optparse import OptionParser
import numpy as np
import itertools
import tensorflow as tf
import keras.backend as K
from keras.models import load_model
from tqdm import tqdm

# Settings
sys.path.insert(0, f'/n/projects/mw2098/shared_code/bpreveal/functions/')
from functional import shuffle_seqs, one_hot_encode_sequence, one_hot_encode_sequences, \
    one_hot_decode_sequence, insert_motif, logitsToProfile
from motifs import extract_seqs_from_df, resize_coordinates
from perturb import generate_random_seq

sys.path.insert(0, f'/home/mw2098/bin/bpreveal/src')
import losses

#Set up options
parser = OptionParser()
parser.add_option("--model_dir",
                  help="Directory of the BPReveal model")
parser.add_option("--tasks_separated_by_commas",
                  help="Comma separated tasks aligned with the given BPReveal model")
parser.add_option("--output_tsv_path",
                  help="Output .tsv.gz filepath to save the predictions to")
parser.add_option("--motifA_metadata_tsv",
                  help=".tsv file of [motif, seq] columns to loop through to inject motifA")
parser.add_option("--motifB_metadata_tsv",
                  help=".tsv file of [motif, seq] columns to loop through to inject motifB")
parser.add_option("--motifA_enhancer_output_window_position", type = "int",
                  help="WRT the output window, where is the motifA to be injected, wrt motif center?")
parser.add_option("--motifB_enhancer_output_window_position", type = "int",
                  help="WRT the output window, where is the motifB to be injected, wrt motif center?")
parser.add_option("--enhancer_chrom", type = "str",
                  help="Enhancer chromosome name")
parser.add_option("--enhancer_output_start_position", type = "int",
                  help="Enhancer start coordinate, wrt output length")
parser.add_option("--enhancer_output_end_position", type = "int",
                  help="Enhancer end coordinate, wrt output length")
parser.add_option("-f", "--fasta_file",
                  help="Path to the .fasta file to extract sequences from.")
parser.add_option("--null_trials", default = 16, type = "int",
                  help="Number of null trials to average mutated motifs across [default: %default]")
parser.add_option("--input_length", default = 2114, type = "int",
                  help="Input sequence length of model. [default: %default]")
parser.add_option("--output_length", default = 1000, type = "int",
                  help="Output sequence length of model. [default: %default]")
(options, args) = parser.parse_args()

def bpreveal_generate_insilico_perturbs_via_motif_affinities_injected_into_enhancer(
model_dir,
tasks_separated_by_commas,
output_tsv_path,
motifA_metadata_tsv,
motifB_metadata_tsv,
motifA_enhancer_output_window_position,
motifB_enhancer_output_window_position,
enhancer_chrom,
enhancer_output_start_position,
enhancer_output_end_position,
fasta_file,
null_trials = 16,
input_length = 2114,
output_length = 1000
):
    assert (enhancer_output_end_position - enhancer_output_start_position)==output_length

    #Import information and arrange enhancer sequence
    enhancer_df = pd.DataFrame([enhancer_chrom, enhancer_output_start_position, enhancer_output_end_position, '*']).transpose()
    enhancer_df.columns = ['chrom', 'start', 'end', 'strand']
    enhancer_input_df = resize_coordinates(enhancer_df, input_length, 'center')
    seqs = extract_seqs_from_df(coords_df = enhancer_input_df, fasta_path = fasta_file,
                                chrom_column = 'chrom', start_column = 'start', end_column = 'end')

    #Get motif positions relative to the input seqlen
    flank = np.floor((input_length - output_length) // 2)
    motifA_enhancer_input_window_position = int(motifA_enhancer_output_window_position + flank)
    motifB_enhancer_input_window_position = int(motifB_enhancer_output_window_position + flank)

    #Read in required data
    tasks = tasks_separated_by_commas.split(',')
    motifA_seqs = pd.read_csv(motifA_metadata_tsv, sep = '\t').reset_index()
    motifB_seqs = pd.read_csv(motifB_metadata_tsv, sep = '\t').reset_index()
    motifA_width = len(motifA_seqs['seq'].values[0])
    motifB_width = len(motifB_seqs['seq'].values[0])

    output_columns = ['state', 'motifA', 'motifA_seq', 'motifB', 'motifB_seq']
    ##################################################
    #Prepare enhancer predictions with "null" features (shuffle null examples 16 times).
    ##################################################
    null_seqs = []
    null_df = pd.DataFrame([['null']*null_trials, ['none']*null_trials, ['none']*null_trials, ['none']*null_trials, ['none']*null_trials]).transpose()
    for t in range(null_trials):
        motifA_null = generate_random_seq(motifA_width)
        motifB_null = generate_random_seq(motifB_width)
        null_seq = insert_motif(seq =  seqs[0], motif = motifA_null, position = motifA_enhancer_input_window_position)
        null_seq = insert_motif(seq = null_seq, motif = motifB_null, position = motifB_enhancer_input_window_position)
        null_seqs = null_seqs + [null_seq]
    null_seqs_1he = one_hot_encode_sequences(null_seqs)
    null_df.columns = output_columns

    #Predict results
    K.clear_session()
    model = load_model(model_dir, custom_objects = {'multinomialNll' : losses.multinomialNll})
    null_raw_arr = model.predict(null_seqs_1he)

    #Collect counts predictions
    for k,task in enumerate(tasks):
        assert len(null_raw_arr[k + len(tasks)].shape)==2
        null_df[task] = null_raw_arr[k + len(tasks)]
    null_df = null_df.groupby(output_columns)[tasks].mean().reset_index()

    ##################################################
    #Prepare enhancer predictions with "A" features (shuffle null examples 16 times).
    ##################################################
    for t in range(null_trials):
        A_seqs = []
        A_seqs_n = motifA_seqs.shape[0]
        A_left_df = pd.DataFrame()
        A_right_df = pd.DataFrame()

        for i,row in motifA_seqs.iterrows():
            left_df = pd.DataFrame([['A']*1, [row.motif]*1, [row.seq]*1, ['none']*1, ['none']*1]).transpose()
            right_df = pd.DataFrame([['B']*1, ['none']*1, ['none']*1, [row.motif]*1, [row.seq]*1]).transpose()
            A_seq = insert_motif(seq = null_seqs[t], motif = row.seq, position = motifA_enhancer_input_window_position)
            A_seqs = A_seqs + [A_seq]
            A_left_df = pd.concat([A_left_df, left_df])
            A_right_df = pd.concat([A_right_df, right_df])

        A_seqs_1he = one_hot_encode_sequences(A_seqs)
        A_left_df.columns = output_columns
        A_right_df.columns = output_columns

        #Predict results
        K.clear_session()
        model = load_model(model_dir, custom_objects = {'multinomialNll' : losses.multinomialNll})
        preds_raw_arr = model.predict(A_seqs_1he)

        #Collect counts predictions
        for k,task in enumerate(tasks):
            assert len(preds_raw_arr[k + len(tasks)].shape)==2
            A_left_df[task] = preds_raw_arr[k + len(tasks)]
            A_right_df[task] = preds_raw_arr[k + len(tasks)]
        A_preds_df = pd.concat([A_left_df, A_right_df])
        A_preds_df = A_preds_df.groupby(output_columns)[tasks].mean().reset_index()

    ##################################################
    #Prepare enhancer predictions with "A" features (shuffle null examples 16 times).
    ##################################################
    for t in range(null_trials):
        B_seqs = []
        B_seqs_n = motifB_seqs.shape[0]
        B_left_df = pd.DataFrame()
        B_right_df = pd.DataFrame()

        for i,row in motifB_seqs.iterrows():
            left_df = pd.DataFrame([['A']*1, [row.motif]*1, [row.seq]*1, ['none']*1, ['none']*1]).transpose()
            right_df = pd.DataFrame([['B']*1, ['none']*1, ['none']*1, [row.motif]*1, [row.seq]*1]).transpose()
            B_seq = insert_motif(seq = null_seqs[t], motif = row.seq, position = motifB_enhancer_input_window_position)
            B_seqs = B_seqs + [B_seq]
            B_left_df = pd.concat([B_left_df, left_df])
            B_right_df = pd.concat([B_right_df, right_df])

        B_seqs_1he = one_hot_encode_sequences(B_seqs)
        B_left_df.columns = output_columns
        B_right_df.columns = output_columns

        #Predict results
        K.clear_session()
        model = load_model(model_dir, custom_objects = {'multinomialNll' : losses.multinomialNll})
        preds_raw_arr = model.predict(B_seqs_1he)

        #Collect counts predictions
        for k,task in enumerate(tasks):
            assert len(preds_raw_arr[k + len(tasks)].shape)==2
            B_left_df[task] = preds_raw_arr[k + len(tasks)]
            B_right_df[task] = preds_raw_arr[k + len(tasks)]
        B_preds_df = pd.concat([B_left_df, B_right_df])
        B_preds_df = B_preds_df.groupby(output_columns)[tasks].mean().reset_index()

    ##################################################
    #Prepare enhancer predictions with "A and B" features, no shuffling needed.
    ##################################################
    AB_seqs = []
    AB_left_df = pd.DataFrame()
    AB_right_df = pd.DataFrame()
    for i,rowA in motifA_seqs.iterrows():
        for j,rowB in motifB_seqs.iterrows():
            left_df =  pd.DataFrame([['AB']*1, [rowA.motif]*1, [rowA.seq]*1, [rowB.motif]*1, [rowB.seq]*1]).transpose()
            right_df = pd.DataFrame([['AB']*1, [rowB.motif]*1, [rowB.seq]*1, [rowA.motif]*1, [rowA.seq]*1]).transpose()
            AB_seq = insert_motif(seq = seqs[0], motif = rowA.seq, position = motifA_enhancer_input_window_position)
            AB_seq = insert_motif(seq = AB_seq, motif = rowB.seq, position = motifB_enhancer_input_window_position)

            #Add on to collect all sequences
            AB_seqs = AB_seqs + [AB_seq]
            AB_left_df = pd.concat([AB_left_df, left_df])
            AB_right_df = pd.concat([AB_right_df, right_df])

    AB_seqs_1he = one_hot_encode_sequences(AB_seqs)
    AB_left_df.columns = output_columns
    AB_right_df.columns = output_columns
    AB_left_df = AB_left_df.reset_index(drop = True)
    AB_right_df = AB_right_df.reset_index(drop = True)

    #Predict results in chunks so that the GPU can handle it
    # step_size = 100000
    step_size = 100000
    AB_preds_df = pd.DataFrame()
    for i,s in enumerate(range(0, AB_seqs_1he.shape[0], step_size)):
        K.clear_session()
        model = load_model(model_dir, custom_objects = {'multinomialNll' : losses.multinomialNll})
        upper_bound = int(np.min(np.array([(s + step_size), AB_seqs_1he.shape[0]])))
        seqs = AB_seqs_1he[s:upper_bound]
        preds_raw_arr = model.predict(seqs)
        left_df = AB_left_df[s:upper_bound]
        right_df = AB_right_df[s:upper_bound]

        #Collect counts predictions
        for k,task in enumerate(tasks):
            assert len(preds_raw_arr[k + len(tasks)].shape)==2
            left_df[task] = preds_raw_arr[k + len(tasks)]
            right_df[task] = preds_raw_arr[k + len(tasks)]
        preds_df = pd.concat([left_df, right_df])
        AB_preds_df = pd.concat([AB_preds_df, preds_df])

    ##################################################
    #Combine all predictions
    ##################################################
    preds_all_df = pd.concat([null_df, A_preds_df, B_preds_df, AB_preds_df])
    preds_all_df.to_csv(output_tsv_path, sep = '\t', index = False)

#Run featured function
bpreveal_generate_insilico_perturbs_via_motif_affinities_injected_into_enhancer(
model_dir = options.model_dir,
tasks_separated_by_commas = options.tasks_separated_by_commas,
output_tsv_path = options.output_tsv_path,
motifA_metadata_tsv = options.motifA_metadata_tsv,
motifB_metadata_tsv = options.motifB_metadata_tsv,
motifA_enhancer_output_window_position = options.motifA_enhancer_output_window_position,
motifB_enhancer_output_window_position = options.motifB_enhancer_output_window_position,
enhancer_chrom = options.enhancer_chrom,
enhancer_output_start_position = options.enhancer_output_start_position,
enhancer_output_end_position = options.enhancer_output_end_position,
fasta_file = options.fasta_file,
null_trials = options.null_trials,
input_length = options.input_length,
output_length = options.output_length)
