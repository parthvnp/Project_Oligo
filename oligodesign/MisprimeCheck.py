"""
Author: Parth Patel
This file contains functions and classes that are useful for testing the presence
homo and heterodimers in oligos.
Created: 18th June 18
"""

from itertools import groupby
from operator import itemgetter
from collections import OrderedDict
import more_itertools as mit
import numpy as np
import math
from oligodesign import CustomExceptions

class Thermodynamics():

    def __init__(self, sequence, temperature=55, oligo_conc=0.25e-6, monovalent=50e-3, divalent=0.0, correction_type="DEFAULT"):
        self.sequence = sequence.upper()
        self.temperature = temperature + 273.15
        self.sequence_length = len(self.sequence)
        self.monovalent = monovalent
        self.divalent = divalent
        self.oligo_conc = oligo_conc
        if correction_type is not None:
            self.correction_type = correction_type
        else:
            self.correction_type = "DEFAULT"

        self.deltaH_internal = {"AA": -7900,
                              "AC": -8400,
                              "AG": -7800,
                              "AT": -7200,
                              "CA": -8500,
                              "CC": -8000,
                              "CG": -10600,
                              "CT": -7800,
                              "GA": -8200,
                              "GC": -9200,
                              "GG": -8000,
                              "GT": -8400,
                              "TA": -7200,
                              "TC": -8200,
                              "TG": -8500,
                              "TT": -7900,
                              }

        self.deltaS_internal = {"AA": -22.2,
                          "AC": -22.4,
                          "AG": -21.0,
                          "AT": -20.4,
                          "CA": -22.7,
                          "CC": -19.9,
                          "CG": -27.2,
                          "CT": -21.0,
                          "GA": -22.2,
                          "GC": -24.4,
                          "GG": -19.9,
                          "GT": -22.4,
                          "TA": -21.3,
                          "TC": -22.2,
                          "TG": -22.7,
                          "TT": -22.2,
                          }

        self.deltaH_external = {"A": 2300,
                              "T": 2300,
                              "C": 100,
                              "G": 100,
                              }

        self.deltaS_external = {"A": 4.1,
                          "T": 4.1,
                          "C": -2.8,
                          "G": -2.8
                          }

    def __str__(self):
        return self.sequence

    def dna_complement(self, sequence):
        complement_rules = {"A": "T", "T": "A", "G":"C", "C":"G"}
        complement = "".join([complement_rules[base] for base in sequence])
        return complement

    def reverse_complement(self, sequence):
        reverse_complement = self.dna_complement(sequence)[::-1]
        return reverse_complement

    def GC_content(self,sequence):
        content  = 0
        for base in sequence:
            if base == "G":
                content = content + 1
            elif base == "C":
                content = content + 1
        content_percentage = content/self.sequence_length
        return content_percentage

    def Tm_deltaH(self, sequence):
        TmH_internal = []
        i = 0
        j = 1
        for s in sequence:
            if j < len(sequence):
                Tm_internal_value = self.deltaH_internal[sequence[i:j+1]]
                TmH_internal.append(Tm_internal_value)
                i = i + 1
                j = j + 1

        TmH_external = [self.deltaH_external[sequence[0]], self.deltaH_external[sequence[-1]]]
        deltaH = sum(TmH_internal) + sum(TmH_external)
        return deltaH

    def Tm_deltaS(self, sequence):
        TmS_internal = []
        i = 0
        j = 1
        for s in sequence:
            if j < len(sequence):
                Tm_internal_value = self.deltaS_internal[sequence[i:j+1]]
                TmS_internal.append(Tm_internal_value)
                i = i + 1
                j = j + 1
        TmS_external = [self.deltaS_external[sequence[0]], self.deltaS_external[sequence[-1]]]
        deltaS = sum(TmS_internal) + sum (TmS_external)
        return deltaS

    def Tm_DNA(self,sequence):

        TmH_internal = []
        temp_i = 0
        temp_j = 1
        for s in sequence:
            if temp_j < len(sequence):
                Tm_internal_value = self.deltaH_internal[sequence[temp_i:temp_j+1]]
                TmH_internal.append(Tm_internal_value)
                temp_i = temp_i + 1
                temp_j = temp_j + 1
        TmH_external = [self.deltaH_external[sequence[0]], self.deltaH_external[sequence[-1]]]
        Tm_deltaH = sum(TmH_internal) + sum(TmH_external)

        TmS_internal = []
        temp_i = 0
        temp_j = 1
        for s in sequence:
            if temp_j < len(sequence):
                Tm_internal_value = self.deltaS_internal[sequence[temp_i:temp_j+1]]
                TmS_internal.append(Tm_internal_value)
                temp_i = temp_i + 1
                temp_j = temp_j + 1

        TmS_external = [self.deltaS_external[sequence[0]], self.deltaS_external[sequence[-1]]]
        Tm_deltaS = sum(TmS_internal) + sum (TmS_external)


        if self.monovalent == 1:
            Tm_celsius = Tm_deltaH / (Tm_deltaS + 1.987 * (math.log(self.oligo_conc)))
            Tm_celsius = Tm_celsius - 273.15
            return Tm_celsius
        elif self.monovalent != 1:
            if self.correction_type == "DEFAULT":
                    Tm_deltaS = Tm_deltaS + (0.368 * ((len(sequence) -1) * math.log(self.monovalent)))
                    Tm_celsius = Tm_deltaH / (Tm_deltaS + 1.987 * ((math.log(self.oligo_conc))))
                    corrected_Tm_celsius = Tm_celsius - 273.15
                    return corrected_Tm_celsius

def is_complement(char1, char2):
    if char1 == "A" and char2 == "T":
        return True
    elif char1 == "T" and char2 == "A":
        return True
    elif char1 == "C" and char2 == "G":
        return True
    elif char1 == "G" and char2 == "C":
        return True
    False

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def consecutive_array(data):
    consecutive_array_length = 0
    for d in data:
        if len(d) >=3:
            consecutive_array_length = consecutive_array_length + len(d)
    return consecutive_array_length

def dna_complement(sequence,reverse=False):
    complement_rules = {"A": "T", "T": "A", "G":"C", "C":"G"}
    complement = "".join([complement_rules[base] for base in sequence])
    if reverse == False:
        return complement
    else:
        return complement[::-1]



def dna_heterodimers_score(sequence1, sequence2, scores_returned=1):
    """Scores for the sequences are derived by calculating the length of consecutive base pairs """

    global_sequence_matrix = {}
    top_score_seqs = {}
    sequence_of_interest = {}

    # Determine the sliding conditions

    if len(sequence1) == len(sequence2):
        max_len_seq = sequence1
        min_len_seq = sequence2

    elif len(sequence1) > len(sequence2):
        max_len_seq  = sequence1
        min_len_seq = sequence2

    elif len(sequence2) > len(sequence1):
        max_len_seq = sequence2
        min_len_seq = sequence1

    # Initiate the sliding algorithm
    min_len_start = -1
    while min_len_start >= -len(min_len_seq):
        min_sub_seq = min_len_seq[min_len_start:]
        max_sub_seq = max_len_seq[0:len(min_sub_seq)]
        remainder_min_sub_seq = min_len_seq[0: min_len_start]
        remainder_max_sub_seq = max_len_seq[len(min_sub_seq):]
        score_matrix = []
        mismatch_score = 0
        gaps_min = len(remainder_min_sub_seq)
        gaps_min = " " * gaps_min
        interaction = [" "] * len(min_sub_seq)
        for min_char, max_char in zip(range(len(min_sub_seq)), range(len(max_sub_seq))):
            if is_complement(min_sub_seq[min_char], max_sub_seq[max_char]):
                score_matrix.append(min_char)
                interaction[min_char] = "|"
            else:
                mismatch_score = mismatch_score + 1
                interaction[min_char] = "x"
        interaction = "".join(interaction)
        global_sequence = f"{gaps_min}{max_sub_seq}{remainder_max_sub_seq}\n{gaps_min}{interaction}\n{remainder_min_sub_seq}{min_sub_seq}"
        if len(score_matrix) >=3:
            analysis_score_matrix = np.array(score_matrix, dtype=np.int16)
            consecutive_matrix = consecutive(analysis_score_matrix)
            consecutive_array_length = consecutive_array(data=consecutive(analysis_score_matrix))
            if consecutive_array_length >=3:
                analysis_seq = max_sub_seq  + "|" + min_sub_seq
                sequence_of_interest[analysis_seq] = consecutive_array_length, consecutive_matrix, mismatch_score
                global_sequence_matrix[global_sequence] = consecutive_array_length, score_matrix

        min_len_start = min_len_start - 1

    max_len_start = 1
    while max_len_start < len(max_len_seq):
        max_sub_seq_2 = max_len_seq[max_len_start: max_len_start + len(min_len_seq)]
        min_sub_seq_2 = min_len_seq[0:len(max_sub_seq_2)]
        remainder_max_sub_seq_2 = max_len_seq[0: max_len_start]
        remainder_min_sub_seq_2 = min_len_seq[len(max_sub_seq_2):]
        score_matrix = []
        mismatch_score = 0
        gaps_min = len(remainder_max_sub_seq_2)
        gaps_min = " " * gaps_min
        interaction = [""] * len(min_sub_seq_2)
        for min2_char, max2_char in zip(range(len(min_sub_seq_2)), range(len(max_sub_seq_2))):

            if is_complement(min_sub_seq[min2_char], max_sub_seq_2[max2_char]):
                score_matrix.append(min2_char)
                interaction[min2_char] = "|"
            else:
                mismatch_score = mismatch_score + 1
                interaction[min2_char] = "x"

        interaction = "".join(interaction)
        global_sequence = f"{remainder_max_sub_seq_2}{max_sub_seq_2}\n{gaps_min}{interaction}\n{gaps_min}{min_sub_seq_2}{remainder_min_sub_seq_2}"
        if len(score_matrix) >=3:
            analysis_score_matrix = np.array(score_matrix, dtype=np.int16)
            consecutive_matrix = consecutive(analysis_score_matrix)
            consecutive_array_length = consecutive_array(data=consecutive(analysis_score_matrix))
            if consecutive_array_length >=3:
                analysis_seq = max_sub_seq_2  + "|" + min_sub_seq_2
                sequence_of_interest[analysis_seq] = consecutive_array_length, consecutive_matrix, mismatch_score
                global_sequence_matrix[global_sequence] = consecutive_array_length, score_matrix

        max_len_start = max_len_start + 1

    top_score_seqs = sorted(sequence_of_interest.items(), reverse=True, key= lambda kv: kv[1][0])[0:scores_returned]
    global_sequence_matrix = dict(sorted(global_sequence_matrix.items(), reverse=True, key = lambda kv: kv[1][0])[0:scores_returned])

    return (top_score_seqs, global_sequence_matrix)

def dna_heterodimers_energy(sequence1, sequence2):
    top_scores, sequnences_scores = dna_heterodimers_score(sequence1=sequence1, sequence2=sequence2)
    top_scores_sequences, top_scores_matrix, mismatch_penalty = [], [], []
    sequence_energy  = {}

    for score_seq in top_scores:
        top_scores_sequences.append(score_seq[0])
        top_scores_matrix.append(score_seq[1][1])
        mismatch_penalty.append(score_seq[1][2])

    for sequence, score_matrix, mismatch_score in zip(top_scores_sequences, top_scores_matrix, mismatch_penalty):
        sense, anti_sense = sequence.split("|")
        deltaH, deltaS = 0 , 0
        complete_sequence = " "
        for seq_array in score_matrix:
            if len(seq_array) >= 3:
                sequence_segment = Thermodynamics(sense[seq_array[0]: seq_array[-1] + 1])
                sequence_deltaH = sequence_segment.Tm_deltaH(sequence=str(sequence_segment))
                sequence_deltaS = sequence_segment.Tm_deltaS(sequence=str(sequence_segment))
                deltaH = deltaH + sequence_deltaH
                deltaS = deltaS + sequence_deltaS
                complete_sequence = complete_sequence+str(sequence_segment)

        deltaG = deltaH - (sequence_segment.temperature * (sequence_segment.Tm_deltaS(sequence=str(sequence_segment))) + (mismatch_score * (438 + 509)))
        deltaG = deltaG/1000
        if sequence_segment.monovalent == 1:
            melting_temperature = (deltaH + ( mismatch_score * (438 + 509) ) ) / (deltaS + (1.987*math.log(sequence_segment.oligo_conc)))
            melting_temperature = melting_temperature - 273.15
        elif sequence_segment.monovalent != 1:
            deltaS = deltaS + (0.368 * ((len(str(sequence_segment)) -1) * math.log(sequence_segment.monovalent)))
            Tm_celsius = deltaH / (deltaS + 1.987 * ((math.log(sequence_segment.oligo_conc))))
            melting_temperature = Tm_celsius - 273.15

        # A generalized penalty is applied for internal mismatches while calculating melting temperature and deltaG values for homo and heterodimers.
        # This penalty value is obtained from the following paper:
        # AutoDimer: a screening tool for primer-dimer and hairpin structures Biotechniques Vol 37, 2
        sequence_energy[complete_sequence] = deltaG, melting_temperature
    return sequence_energy,sequnences_scores

def primer_dimer_representation(sequence1, sequence2, Tm_value):
    sequence_thermodynamics, sequence_scores = dna_heterodimers_energy(sequence1=sequence1, sequence2=sequence2)
    sequence_score_matrix_score = []
    representation = ""
    for value in sequence_scores.values():
        sequence_score_matrix_score.append(value[0])
    for sequence, score, energy in zip(sequence_scores.keys(), sequence_score_matrix_score, sequence_thermodynamics.values()):
        if energy[1] >= Tm_value:
            representation = representation + f"Score:{score}  deltaG:{energy[0]}  Tm:{energy[1]}\n{sequence}\n\n"
    return representation


def homodimer_check(sequence, Tm_value):
    sequence1 = sequence
    sequence2 = dna_complement(sequence1)
    sequence2 = sequence2[::-1]
    outcome = primer_dimer_representation(sequence1=sequence1, sequence2=sequence2, Tm_value=Tm_value)
    return outcome

def main_misprime_check(sequence1, sequence2, Tm_value, heterodimer=True):
    try:
        if heterodimer == True:
            return primer_dimer_representation(sequence1=sequence1, sequence2=sequence2, Tm_value=Tm_value)
        else:
            return homodimer_check(sequence=sequence1)

    except Exception as e:

        return (f"An error occured while processing the input.")

