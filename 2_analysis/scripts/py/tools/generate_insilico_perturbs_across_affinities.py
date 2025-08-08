#!/home/mw2098/anaconda3/envs/bpreveal_404/bin/python

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
import pandas as pd
import keras
import json
from optparse import OptionParser
import numpy as np
import itertools
import tensorflow as tf
import keras.backend as K
from keras.models import load_model
from tqdm import tqdm

# Settings
sys.path.insert(0, f'/n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/scripts/py/functions')
from functional import shuffle_seqs, one_hot_encode_sequence, one_hot_encode_sequences, \
    one_hot_decode_sequence, insert_motif, logitsToProfile
from motifs import extract_seqs_from_df, resize_coordinates
sys.path.insert(0, f'/n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404/src')
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
parser.add_option("--seqA_txt_filepath",
                  help="Filepath of tab-separated file that contains the column 'seq' to represent all motif affinities to be injected.")
parser.add_option("--motifA_name",
                  help="Name of the anchored motif to be injected")
parser.add_option("--motifA_output_position", default = 470, type = "int",
                  help="WRT the output window, where is the anchored motif to be injected?")
parser.add_option("--seqB_txt_filepath",
                  help="Filepath of tab-separated file that contains the column 'seq' to represent all motif affinities to be injected.")
parser.add_option("--motifB_name",
                  help="Name of the anchored motif to be injected")
parser.add_option("--motifB_output_position", default = 530, type = "int",
                  help="WRT the output window, where is the anchored motif to be injected?")
parser.add_option("--input_length", default = 2114, type = "int",
                  help="Input sequence length of model. [default: %default]")
parser.add_option("--output_length", default = 1000, type = "int",
                  help="Output sequence length of model. [default: %default]")
parser.add_option("--null_sequence_np_prefix", default = 'seqs_1he', type = "str",
                  help="Npz subarray prefix to access the 1he seqquences [default: %default]")
(options, args) = parser.parse_args()

def bpreveal_generate_insilico_perturbs_across_affinities(
model_dir,
null_sequence_np_filepath,
tasks_separated_by_commas,
output_tsv_path,
seqA_txt_filepath,
motifA_name,
seqB_txt_filepath,
motifB_name,
motifA_output_position = 470,
motifB_output_position = 530,
input_length = 2114,
output_length = 1000,
null_sequence_np_prefix = 'seqs_1he'
):
    #Read in required data
    null_seqs = np.load(null_sequence_np_filepath)[null_sequence_np_prefix]
    motifA_seqs_df = pd.read_csv(seqA_txt_filepath, sep = '\t')
    if motifB_name != 'empty':
        motifB_seqs_df = pd.read_csv(seqB_txt_filepath, sep = '\t')
    else:
        print('Marginalization without motif pairs running...')
    tasks = tasks_separated_by_commas.split(',')
    motifA_input_position = int((input_length-output_length)/2 + motifA_output_position)
    motifB_input_position = int((input_length-output_length)/2 + motifB_output_position)
    colnames_of_df = ['seqA', 'motifA', 'seqB', 'motifB', 'state']

    #Collect predictions by trial
    preds_by_trial_df = pd.DataFrame()
    for t in tqdm(range(null_seqs.shape[0])):

        #Prepare null sequence and anchored sequence information
        null_seq = one_hot_decode_sequence(null_seqs[t])
        injected_seqs = [null_seq]
        left_is_motifA_df = pd.DataFrame([['none']*1, [motifA_name]*1, ['none']*1, [motifB_name]*1, ['null']]).transpose()
        right_is_motifA_df = pd.DataFrame([['none']*1, [motifB_name]*1, ['none']*1, [motifA_name]*1, ['null']]).transpose()

        for sA in motifA_seqs_df.seq.values:
            seqA = insert_motif(seq = null_seq, motif = sA, position = motifA_input_position)
            injected_seqs = injected_seqs + [seqA]

            left_df = pd.DataFrame([[sA]*1, [motifA_name]*1, ['none']*1, [motifB_name]*1, ['A']]).transpose()
            left_is_motifA_df = pd.concat([left_is_motifA_df, left_df])
            right_df = pd.DataFrame([['none']*1, [motifB_name]*1, [sA]*1, [motifA_name]*1, ['B']]).transpose()
            right_is_motifA_df = pd.concat([right_is_motifA_df, right_df])

        if motifB_name != 'empty':
            for sB in motifB_seqs_df.seq.values:
                seqB = insert_motif(seq = null_seq, motif = sB, position = motifB_input_position)
                injected_seqs = injected_seqs + [seqB]

                left_df = pd.DataFrame([['none']*1, [motifA_name]*1, [sB]*1, [motifB_name]*1, ['B']]).transpose()
                right_df = pd.DataFrame([[sB]*1, [motifB_name]*1, ['none']*1, [motifA_name]*1, ['A']]).transpose()
                left_is_motifA_df = pd.concat([left_is_motifA_df, left_df])
                right_is_motifA_df = pd.concat([right_is_motifA_df, right_df])
            
            for sA in motifA_seqs_df.seq.values:
                for sB in motifB_seqs_df.seq.values:
                    seqAB = insert_motif(seq = null_seq, motif = sA, position = motifA_input_position)
                    seqAB = insert_motif(seq = seqAB, motif = sB, position = motifB_input_position)
                    injected_seqs = injected_seqs + [seqAB]

                    left_df = pd.DataFrame([[sA]*1, [motifA_name]*1, [sB]*1, [motifB_name]*1, ['AB']]).transpose()
                    right_df = pd.DataFrame([[sB]*1, [motifB_name]*1, [sA]*1, [motifA_name]*1, ['AB']]).transpose()
                    left_is_motifA_df = pd.concat([left_is_motifA_df, left_df])
                    right_is_motifA_df = pd.concat([right_is_motifA_df, right_df])

        injected_seqs_1he = one_hot_encode_sequences(injected_seqs)
        left_is_motifA_df.columns = colnames_of_df
        right_is_motifA_df.columns = colnames_of_df
        

        #Predict results in chunks so that the GPU can handle it
        # step_size = 100000
        print('Predicting sequences...')
        step_size = 100000
        all_preds_df = pd.DataFrame()
        for s in range(0, injected_seqs_1he.shape[0], step_size):
            K.clear_session()
            model = load_model(model_dir, custom_objects = {'multinomialNll' : losses.multinomialNll, 'reweightableMse': losses.dummyMse})
            upper_bound = int(np.min(np.array([(s + step_size), injected_seqs_1he.shape[0]])))
            seqs = injected_seqs_1he[s:upper_bound]
            preds_raw_arr = model.predict(seqs)
            left_df = left_is_motifA_df[s:upper_bound]
            right_df = right_is_motifA_df[s:upper_bound]

            #Collect counts predictions
            for k,task in enumerate(tasks):
                assert len(preds_raw_arr[k + len(tasks)].shape)==2
                left_df[task] = preds_raw_arr[k + len(tasks)]
                right_df[task] = preds_raw_arr[k + len(tasks)]
            preds_df = pd.concat([left_df, right_df])
            if motifB_name != 'empty':
                all_preds_df = pd.concat([all_preds_df, preds_df])
            else:
                all_preds_df = pd.concat([all_preds_df, left_df])
        all_preds_df['trial'] = t
        preds_by_trial_df = pd.concat([preds_by_trial_df, all_preds_df])

    print('Writing predictions...')
    preds_all_df = preds_by_trial_df.groupby(colnames_of_df)[tasks].mean().reset_index()
    preds_all_df.to_csv(output_tsv_path, sep = '\t', index = False)

#Run featured function
bpreveal_generate_insilico_perturbs_across_affinities(
model_dir = options.model_dir,
null_sequence_np_filepath = options.null_sequence_np_filepath,
tasks_separated_by_commas = options.tasks_separated_by_commas,
output_tsv_path = options.output_tsv_path,
seqA_txt_filepath = options.seqA_txt_filepath,
motifA_name = options.motifA_name,
seqB_txt_filepath = options.seqB_txt_filepath,
motifB_name = options.motifB_name,
motifA_output_position = options.motifA_output_position,
motifB_output_position = options.motifB_output_position,
input_length = options.input_length,
output_length = options.output_length,
null_sequence_np_prefix = options.null_sequence_np_prefix)
