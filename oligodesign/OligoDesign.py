import numpy as np
import re
import os
import traceback as tb
import time
import math
from oligodesign import CustomExceptions
from oligodesign import MisprimeCheck
from colorama import Fore, Style, Back
import time
from random import choice
ABSOLUTE_TEMPERATURE = 273.15

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    COLORFUL_FORWARD = "\033[0;30;41m"
    COLORFUL_REVERSE = " \033[0;37;45m"

class EquiTmOligo(object):

    def __init__(self, sequence, min_tm=55.0, min_oligo_length=15, max_oligo_length=60, oligo_conc=0.25e-6, monovalent=50e-3, divalent=0.0, correction_type="DEFAULT"):
        self.sequence = sequence.upper()
        self.sequence_length = len(self.sequence)
        self.min_tm = min_tm
        self.min_oligo_length = min_oligo_length
        self.max_oligo_length  = max_oligo_length
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

    def sequence_validation(self):
        is_valid = True
        try:

            for base in self.sequence:
                if base not in "ATCGU":
                    raise CustomExceptions.InvalidSequence

            if self.sequence_length < 60 or self.sequence_length > 600:
                raise CustomExceptions.InvalidLength

            if self.min_oligo_length > self.max_oligo_length:
                raise CustomExceptions.InvalidLength

            if self.min_tm < 55 or self.min_tm > 70:
                raise CustomExceptions.InvalidTemperature

        except CustomExceptions.InvalidSequence:
            print("Invalid DNA sequence.")
            is_valid = False

        except CustomExceptions.InvalidLength:
            print("The length of provided is invalid. The minimum acceptable length is 60 and the maximum acceptable length is 600.")
            is_valid = False

        except CustomExceptions.InvalidTemperature:
            print("The acceptable temperature range for the program is 55-70" + u"\u00B0" + ".")
            is_valid = False

        return is_valid


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


    def Tm_DNA(self,sequence):

        # TODO: Add support for creating oligos within user specificed length.

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
                    Tm_celsius = Tm_deltaH / (Tm_deltaS + 1.987 * (math.log(self.oligo_conc)))
                    Tm_celsius = Tm_celsius - 273.15
                    Tm_deltaS = Tm_deltaS + (0.368 * ((len(sequence) -1) * math.log(self.monovalent)))
                    Tm_celsius = Tm_deltaH / (Tm_deltaS + 1.987 * ((math.log(self.oligo_conc))))
                    corrected_Tm_celsius = Tm_celsius - 273.15
                    return corrected_Tm_celsius
            # TODO: Include Support for Owczry type ionic correction
            elif self.correction_type == "OWCZARY":
                Tm_celsius = Tm_deltaH / (Tm_deltaS + 1.987 * (math.log(self.oligo_conc)))
                Tm_celsius = Tm_celsius - 273.15
                f_GC = self.GC_content(sequence=sequence)
                Tm = Tm_deltaH / (Tm_deltaS + 1.987 * (math.log(self.oligo_conc)))
                Tm_celsius = Tm_celsius - 273.15
                corrected_Tm = (1/Tm_celsius) + (
                    (4.29*f_GC - 3.95) *
                    (math.log(self.monovalent)) +
                    (0.940*math.pow(math.log(self.monovalent),2))
                    * 1e-5
                )
                return corrected_Tm

    def primer_index_generation(self, alter_avg_size=None):
        if alter_avg_size:
            avg_size = alter_avg_size
        else:
            avg_size = int((self.max_oligo_length + self.min_oligo_length)/2)
        if self.min_tm < 60:
            overlap = 10
        elif self.min_tm < 70:
            overlap = 12
        else:
            overlap = 14

        minimum_oligos = math.ceil((self.sequence_length - avg_size)/(avg_size - overlap))
        minimum_oligos = minimum_oligos + 1

        try:
            if minimum_oligos <= 2:
                raise CustomExceptions.InvalidLength
        except CustomExceptions.InvalidLength:
            print("The sequence must be long enough to generate at least two oligos.")

        temp_i = 0
        temp_j = avg_size
        is_success = True
        primers_index = np.zeros(shape=(minimum_oligos, 3), dtype=np.int16)
        for j in range(1, minimum_oligos + 1):
            if j % 2 != 0:
                primers_index[j-1, 0] = temp_i
                primers_index[j-1, 1] = temp_j
                primers_index[j-1, 2] = 1
                temp_i = temp_j - overlap
                temp_j = temp_i + avg_size
            if j % 2 == 0:
                primers_index[j-1, 0] = temp_i
                primers_index[j-1, 1] = temp_j
                primers_index[j-1, 2] = -1
                temp_i = temp_j - overlap
                temp_j = temp_i + avg_size
        if primers_index[-1, 1] == self.sequence_length:
            primers_index[-1,1] = self.sequence_length - 1
            Warning = "NA"
            return (Warning, primers_index)
        elif primers_index[-1,1] != self.sequence_length:
            primers_index[-1,1] = self.sequence_length - 1
            if primers_index[-1,1] - primers_index[-1, 0] < self.min_oligo_length:
                Warning = 0
                return (Warning, primers_index)
            elif primers_index[-1, 1] - primers_index[-1,0] > self.max_oligo_length:
                Warning = 1
                return (Warning, primers_index)
            elif (primers_index[-1,1] - primers_index[-1,0] >= self.min_oligo_length) and (primers_index[-1,1] - primers_index[-1,0] <= self.max_oligo_length):
                Warning = "PERFECT"
                return (Warning, primers_index)

    def primers_index_readjustment(self):
        output, array_index = self.primer_index_generation()

        if output == 0:
            try:
                avg_size = int((self.max_oligo_length + self.min_oligo_length)/2)
                avg_alt_size = avg_size + 1
                while output == 0 and avg_alt_size >= avg_size + 5:
                    output, primer_array = self.primer_index_generation(alter_avg_size=avg_alt_size)
                    avg_alt_size = avg_alt_size + 1
                return primer_array

            except:
                pass

        if output == 1:
            try:
                avg_size = int((self.max_oligo_length + self.min_oligo_length)/2)
                avg_alt_size = avg_size - 1
                while output == 0 and avg_alt_size >= avg_size - 5:
                    output, primer_array = self.primer_index_generation(alter_avg_size=avg_alt_size)
                    avg_alt_size = avg_alt_size - 1
                return primer_array

            except:
                pass

        if output == "NA" or "PERFECT":
            return array_index


    def dna_fragmentation(self):
        primers_index = self.primers_index_readjustment()
        oligo_overlaps = []
        if primers_index.shape[0] >= 2:
            for i in range(primers_index.shape[0]-1):
                overlap = [primers_index[i+1,0], primers_index[i,1]]
                oligo_overlaps.append(overlap)

        oligos_tm = {}
        oligos_overlap_index = []
        default_tm = self.min_tm
        for oligo_index in oligo_overlaps:
            temp_i = 1
            melt_tm = self.Tm_DNA(sequence=self.sequence[int(oligo_index[0]): int(oligo_index[1]+1)])
            if melt_tm > default_tm:
                while default_tm + 3 <= melt_tm:
                    if temp_i%2 != 0:
                        oligo_index[0] = oligo_index[0] + 1
                        melt_tm = self.Tm_DNA(sequence=self.sequence[int(oligo_index[0]): int(oligo_index[1]+1)])
                        temp_i = temp_i + 1
                    elif temp_i%2 == 0:
                        if oligo_index[1] < self.sequence_length-1:
                            oligo_index[1] = oligo_index[1] - 1
                            melt_tm = self.Tm_DNA(sequence=self.sequence[int(oligo_index[0]):int(oligo_index[1]+1)])
                        temp_i = temp_i + 1
                oligos_overlap_index.append([oligo_index[0], oligo_index[1]])
                oligos_tm[self.sequence[oligo_index[0], oligo_index[1]+1]] = melt_tm

            elif melt_tm < default_tm:
                while melt_tm <= default_tm - 3:
                    if temp_i % 2 != 0:
                        oligo_index[0] = oligo_index[0] - 1
                        melt_tm = self.Tm_DNA(sequence=self.sequence[int(oligo_index[0]):int(oligo_index[1]+1)])
                        temp_i = temp_i + 1
                    elif temp_i%2 == 0:
                        if oligo_index[1] < self.sequence_length-1:
                            oligo_index[1] = oligo_index[1] + 1
                            melt_tm = self.Tm_DNA(sequence=self.sequence[int(oligo_index[0]): int(oligo_index[1]+1)])
                        temp_i = temp_i + 1
                oligos_overlap_index.append([oligo_index[0], oligo_index[1]])
                oligos_tm[self.sequence[int(oligo_index[0]):int(oligo_index[1] + 1)]] = melt_tm

            elif self.min_tm == melt_tm:
                oligos_overlap_index.append([oligo_index[0], oligo_index[1]])
                sequence_segment = self.sequence[int(oligo_index[0]):int(oligo_index[1]+1)]
                oligos_tm[sequence_segment] = melt_tm
        for i,j in zip(oligos_overlap_index, range(1, len(oligos_overlap_index) + 1)):
            primers_index[j,0] = i[0]
            primers_index[j-1,1] = i[1]
        return (oligos_overlap_index, oligos_tm, primers_index)


    def oligos_representation(self):
        overlap_index, oligo_dictionary, primers_index = self.dna_fragmentation()
        sense_oligos, anti_sense_oligos, five_three_oligos = [], [], []
        oligo_profile = []
        sequence_complement  = self.dna_complement(sequence=self.sequence)
        for oligo in primers_index:
            if oligo[2] == 1:
                sense_oligos.append(self.sequence[oligo[0]: oligo[1] + 1])
                five_three_oligos.append((self.sequence[oligo[0]: oligo[1] + 1], 1))
            if oligo[2] == -1:
                five_three = sequence_complement[oligo[0]: oligo[1]+1][::-1]
                anti_sense_oligos.append(five_three)
                five_three_oligos.append((five_three, -1))

        forward_counter = 1
        reverse_counter = 1

        for oligo in five_three_oligos:
            if oligo[1] == -1:
                oligo_length =  len(oligo[0])
                oligo_melt_Tm = self.Tm_DNA(sequence=oligo[0])
                oligo_GC_content = self.GC_content(sequence=oligo[0])
                oligo_name = str(forward_counter)+"_FORWARD"
                oligo_profile.append((oligo_name, oligo[0], oligo_length, oligo_melt_Tm, oligo_GC_content))
                forward_counter = forward_counter + 1
            elif oligo[1] == 1:
                oligo_length =  len(oligo[0])
                oligo_melt_Tm = self.Tm_DNA(sequence=oligo[0])
                oligo_GC_content = self.GC_content(sequence=oligo[0])
                oligo_name = str(reverse_counter)+"_REVERSE"
                oligo_profile.append((oligo_name, oligo[0],oligo_length, oligo_melt_Tm, oligo_GC_content ))
                reverse_counter = reverse_counter + 1
        print(oligo_profile)
        return oligo_profile


    def dna_representation(self):
        ( overlap_index, overlap_dictionary, primer_set) = self.dna_fragmentation()
        forward, reverse, = [], [],
        sequence_index = []
        forward_primers = []
        reverse_primers = []
        forward_seq = [" "] * self.sequence_length 
        interaction = [" "] * self.sequence_length 
        fwd_arrow = [" "] * self.sequence_length
        rev_arrow = [" "] * self.sequence_length
        sequence_complement = self.dna_complement(self.sequence)
        reverse_seq = [" "] * len(sequence_complement)

        for primer in primer_set:
            seq_start, seq_stop, direction = primer[0], primer[1], primer[2]

            if primer[2] == 1:
                for fwd in range(seq_start, seq_stop+1):
                    forward_seq[fwd] =  self.sequence[fwd]
                    forward_seq[seq_start] = self.sequence[seq_start]
                forward_primers.append(self.sequence[seq_start:seq_stop+1])

            if primer[2] == 1:
                if seq_stop < self.sequence_length:
                    forward_seq[seq_stop] =  self.sequence[seq_stop]

            if primer[2] == -1:
                for rev in range(seq_start, seq_stop+1):
                    reverse_seq[rev] = sequence_complement[rev]
                reverse_primers.append(sequence_complement[seq_start:seq_stop+1])

        counter = 0
        interaction2 = [" "] * self.sequence_length
        for inter, temp in zip(overlap_index, overlap_dictionary.values()):
            inter_start, inter_stop = inter[0], inter[1]
            inter_average = int((inter_start+inter_stop)/2)
            interaction2[inter_average-counter] = " " + f"{temp:.1f}"
            for i in range(inter_start, inter_stop + 1):
                interaction[i] = "|"
            counter = counter + 4

            # interaction[inter_stop] = interaction[inter_stop] + " " + f"{int(temp):.2f}"


        forward_seq = "".join(forward_seq)
        interaction = "".join(interaction)
        reverse_seq = "".join(reverse_seq)
        interaction2 = "".join(interaction2[0: -counter])
        forward_seq = "5'- " + forward_seq + " -3'"
        interaction = "    " + interaction + "    "
        interaction2 = "   " + interaction2 + "     " 
        reverse_seq = "3'- " + reverse_seq + " -5'"

        fwd_primers = ""
        rev_primers = ""
        for index, fwd in enumerate(forward_primers):
            if len(fwd)> self.max_oligo_length:
                print(f"Sequence Length: {len(fwd)} {fwd}")
            assert len(fwd) <= self.max_oligo_length
            fwd_primers = fwd_primers + "\n" + str(index+1) + "F" + "  " + fwd + "   " + str(len(fwd))

        for index, rev in enumerate(reverse_primers):
            if len(rev)> self.max_oligo_length:
                print(f"Sequence Length: {len(rev)} {rev}")
            assert len(rev) <= self.max_oligo_length
            rev_primers = rev_primers + "\n" + str(index+1) + "R" + "  " + rev + "  " + str(len(rev))


        representation = f"{forward_seq}\n{interaction}\n{reverse_seq}\n{interaction2}"
        representation = "\n" + representation + "\n\n" + fwd_primers +"\n" +  rev_primers
        return (forward_seq, interaction, reverse_seq, interaction2)

    def main(self):
        try:
            self.sequence_validation()
            return (self.dna_representation())

        except AssertionError:
            print("Could not find the oligos with the desired parameter")

    def primer_dimer_check(self):
        primer_dimer_list = []
        oligo_list = self.oligos_representation()
        for o1 in range(len(oligo_list)):
            for o2 in range(len(oligo_list)):
                if o1 != o2:
                    reverse_sequence = oligo_list[o2][1][::-1]
    
                    primer_dimer = MisprimeCheck.main_misprime_check(oligo_list[o1][1], reverse_sequence, Tm_value=self.min_tm)
                    if primer_dimer:
                        primer_dimer_list.append(primer_dimer)
                        print(primer_dimer)

        return primer_dimer_list

sequence = EquiTmOligo(sequence="AGCTAGCGATCGCGGCCGACGATCGATCGGCGCATGCAGCTACGATGCGCGCGCTAGCGCTAGCGCGCTAGCGCGCTAGCGCGCTACGCGTACGGCGCGCTAGCGC")
sequence.primer_dimer_check()
