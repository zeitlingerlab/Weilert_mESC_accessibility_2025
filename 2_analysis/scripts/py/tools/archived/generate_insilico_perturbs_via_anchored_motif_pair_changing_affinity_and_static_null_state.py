
"""
Melanie Weilert
May 2023
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
sys.path.insert(0, f'/home/mw2098/bin/bpreveal/src')
import losses

#Set up options
parser = OptionParser()
parser.add_option("--model_dir",
                  help="Directory of the BPReveal model")
parser.add_option("--null_sequence_np_filepath",
                  help="Filepath to the .npz file that contains 1he null sequences to run trials across")
parser.add_option("--tasks_separated_by_commas",
                  help="Comma separated tasks aligned with the given BPReveal model")
parser.add_option("--output_tsv_path",
                  help="Output .tsv.gz filepath to save the predictions to")
parser.add_option("--anchored_seqA",
                  help="Sequence of the anchored motifA to be injected")
parser.add_option("--anchored_nameA",
                  help="Name of the anchored motifA to be injected")
parser.add_option("--anchored_output_positionA", default = 470, type = "int",
                  help="WRT the output window, where is the anchored motifA to be injected?")
parser.add_option("--anchored_seqB",
                  help="Sequence of the anchored motifB to be injected")
parser.add_option("--anchored_nameB",
                  help="Name of the anchored motifB to be injected")
parser.add_option("--anchored_output_positionB", default = 500, type = "int",
                  help="WRT the output window, where is the anchored motifB to be injected?")
parser.add_option("--affinity_txt_filepath",
                  help="Filepath of tab-separated file that contains the column 'seq' to represent all motif affinities to be injected into null.")
parser.add_option("--affinity_name",
                  help="Name of the affinity-changed motif to be injected")
parser.add_option("--affinity_output_position", default = 530, type = "int",
                  help="WRT the output window, where is the affinity-changing motif to be injected?")
parser.add_option("--input_length", default = 2114, type = "int",
                  help="Input sequence length of model. [default: %default]")
parser.add_option("--output_length", default = 1000, type = "int",
                  help="Output sequence length of model. [default: %default]")
parser.add_option("--null_sequence_np_prefix", default = 'seqs_1he', type = "str",
                  help="Npz subarray prefix to access the 1he seqquences [default: %default]")
(options, args) = parser.parse_args()

def bpreveal_generate_insilico_perturbs_via_anchored_motif_pair_and_changing_null_state(
model_dir,
null_sequence_np_filepath,
tasks_separated_by_commas,
output_tsv_path,
anchored_seqA,
anchored_nameA,
anchored_seqB,
anchored_nameB,
affinity_txt_filepath,
affinity_name,
anchored_output_positionA = 470,
anchored_output_positionB = 500,
affinity_output_position = 530,
input_length = 2114,
output_length = 1000,
null_sequence_np_prefix = 'seqs_1he'
):
    #Read in required data
    null_seqs = np.load(null_sequence_np_filepath)[null_sequence_np_prefix]
    tasks = tasks_separated_by_commas.split(',')
    anchored_input_positionA = int((input_length-output_length)/2 + anchored_output_positionA)
    anchored_input_positionB = int((input_length-output_length)/2 + anchored_output_positionB)
    affinity_input_position = int((input_length-output_length)/2 + affinity_output_position)

    #Collect predictions by trial
    preds_by_trial_df = pd.DataFrame()
    for t in tqdm(range(null_seqs.shape[0])):

        #Initialize features
        injected_seqs = []
        left_is_motifA_df = pd.DataFrame()
        right_is_motifA_df = pd.DataFrame()

        #Read in affinity sequences
        affinity_seqs_df = pd.read_csv(affinity_txt_filepath, sep = '\t')

        for s in affinity_seqs_df.seq.values:
            #Prepare null sequence and anchored sequence information
            null_seq = one_hot_decode_sequence(null_seqs[t])

            #seqA has "two" sequences inserted
            seqA = insert_motif(seq = null_seq, motif = anchored_seqA, position = anchored_input_positionA)
            seqA = insert_motif(seq = seqA, motif = anchored_seqB, position = anchored_input_positionB)

            #Seq B has the changing affinity feature
            seqB = insert_motif(seq = null_seq, motif = s, position = affinity_input_position)

            #Seq AB has all three motifs
            seqAB = insert_motif(seq = seqA, motif = s, position = affinity_input_position)
            injected_seqs = injected_seqs + [null_seq, seqA, seqB, seqAB]

            #Create indexes of how the prediction arrays were assigned
            left_df = pd.DataFrame([[s]*4, [anchored_nameA + '_' + anchored_nameB]*4, [affinity_name]*4,  ['null', 'A_grouped', 'B','A_groupedB']]).transpose()
            right_df = pd.DataFrame([[s]*4, [affinity_name]*4, [anchored_nameA + '_'  + anchored_nameB]*4, ['null', 'B', 'A_grouped','A_groupedB']]).transpose()
            left_is_motifA_df = pd.concat([left_is_motifA_df, left_df])
            right_is_motifA_df = pd.concat([right_is_motifA_df, right_df])

        injected_seqs_1he = one_hot_encode_sequences(injected_seqs)
        left_is_motifA_df.columns = ['affinity_seq','motifA', 'motifB', 'state']
        right_is_motifA_df.columns = ['affinity_seq','motifA', 'motifB', 'state']

        #Predict results
        K.clear_session()
        model = load_model(model_dir, custom_objects = {'multinomialNll' : losses.multinomialNll})
        preds_raw_arr = model.predict(injected_seqs_1he)

        #Collect counts predictions
        for k,task in enumerate(tasks):
            assert len(preds_raw_arr[k + len(tasks)].shape)==2
            left_is_motifA_df[task] = preds_raw_arr[k + len(tasks)]
            right_is_motifA_df[task] = preds_raw_arr[k + len(tasks)]

        preds_df = pd.concat([left_is_motifA_df, right_is_motifA_df])
        preds_df['trial'] = t
        preds_by_trial_df = pd.concat([preds_by_trial_df, preds_df])

    preds_all_df = preds_by_trial_df.groupby(['affinity_seq', 'motifA', 'motifB', 'state'])[tasks].mean().reset_index()
    preds_all_df['affinity_name'] = affinity_name
    preds_all_df['anchored_nameA'] = anchored_nameA
    preds_all_df['anchored_nameB'] = anchored_nameB
    preds_all_df['anchored_seqA'] = anchored_seqA
    preds_all_df['anchored_seqB'] = anchored_seqB
    preds_all_df.to_csv(output_tsv_path, sep = '\t', index = False)

#Run featured function
bpreveal_generate_insilico_perturbs_via_anchored_motif_pair_and_changing_null_state(
model_dir = options.model_dir,
null_sequence_np_filepath = options.null_sequence_np_filepath,
tasks_separated_by_commas = options.tasks_separated_by_commas,
output_tsv_path = options.output_tsv_path,
anchored_seqA = options.anchored_seqA,
anchored_nameA = options.anchored_nameA,
anchored_seqB = options.anchored_seqB,
anchored_nameB = options.anchored_nameB,
affinity_txt_filepath = options.affinity_txt_filepath,
affinity_name = options.affinity_name,
anchored_output_positionA = options.anchored_output_positionA,
anchored_output_positionB = options.anchored_output_positionB,
affinity_output_position = options.affinity_output_position,
input_length = options.input_length,
output_length = options.output_length,
null_sequence_np_prefix = options.null_sequence_np_prefix)
