"""

Author: Parth Patel

Implementation of a modifed and enhanced Nussinov Algorithm to predict
secondary structure of oligos and determine their melting temperature and free energy. We also use
RNAstructure's Python bindings for secondary structure analysis.


"""
import numpy as np
import forgi.graph.bulge_graph as fgb
from oligodesign.SecondaryStructureCheck.Zuker import RNAstructure
from itertools import groupby
from operator import itemgetter

class NussinovPredictor:
    def __init__(self, sequence):
        self.seq = sequence
        self.seq_len = len(self.seq)
        self.min_loop_length = 3
        self.score_matrix = [[0 for i in range(self.seq_len)] for j in range(self.seq_len)]
        self.opt_struct = []
        self.score = 0
        self.opt_struct = []
        self.score = 0

    def score_function(self, x, y):
        """
        This function takes in two nucleotides and returns
        a rough free energy value for their base pairing.
        """
        if x>= y - self.min_loop_length:
            self.score = 0
        else:
            if (self.seq[x] == "G" and self.seq[y] == "C") or (self.seq[x] == "C" and self.seq[y] == "G"):
                self.score = 1
            elif (self.seq[x] == "A" and self.seq[y] == "T") or (self.seq[x] == "T" and self.seq[y] == "A"):
                self.score = 1
            else:
                self.score = 0

        return self.score

    def bifurcation(self, x, y):
        tmp_value = 0
        for k in range(x+1, y-1):
            tmp_value = max(tmp_value, self.score_matrix[x][k] + self.score_matrix[k+1][y])
        return tmp_value

    def gamma(self, x, y):
        return (
            max(self.score_matrix[x+1][y],
                    self.score_matrix[x][y-1],
                    self.score_matrix[x+1][y-1] + self.score_function(x,y),
                    self.bifurcation(x,y),)
            )
    def bactrack(self, x, y):
        if x < y:
            if self.score_matrix[x][y] == self.score_matrix[x+1][y]:
                self.bactrack(x+1, y)
            elif self.score_matrix[x][y] == self.score_matrix[x][y-1]:
                self.bactrack(x, y-1)
            elif self.score_matrix[x][y] == self.score_matrix[x+1][y-1] + self.score_function(x,y):
                self.opt_struct.append((x, y))
                self.bactrack(x+1, y-1)
            else:
                for i in range(x+1, y-1):
                    if self.score_matrix[x][y] == self.score_matrix[x][i] + self.score_matrix[i+1][y]:
                        self.bactrack(x, i)
                        self.bactrack(i+1, y)
                        break

    def vienna_format(self, max_score):
        output = ""
        open_brackets = [i[0] for i in self.opt_struct]
        closing_brackets = [i[1] for i in self.opt_struct]

        for j in range(self.seq_len):
            if j in open_brackets:
                output = output + "("
            elif j in closing_brackets:
                output = output + ")"
            else:
                output = output + "."

        return output

    def printMatrix(self):
        for i in range(len(self.score_matrix)):
            print (self.score_matrix[i])

    def main_nussinov(self):
        for k in range(1, self.seq_len):
            for i in range(self.seq_len - k):
                j = i + k
                self.score_matrix[i][j] = self.gamma(i,j)
        self.bactrack(i,j)
        max_score = self.score_matrix[0][self.seq_len - 1]
        return (self.vienna_format(max_score))


class ZukerPredictor(object):

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def mfe_fold(self):
        raw_sequence = RNAstructure.RNA.fromString(self.sequence, backbone="DNA")
        raw_sequence.FoldSingleStrand(maximumstructures=1)
        mfe = raw_sequence.GetFreeEnergy(1)
        stem_pairs = [(i,raw_sequence.GetPair(i)) for i in raw_sequence.iterIndices()]
        output = ""
        open_paren, close_paren, stacking_interactions = [], [], []
        for paren in stem_pairs:
            if paren[1] != 0 and paren[1] > paren[0]:
                open_paren.append(paren[0] - 1)
                close_paren.append(paren[1] - 1)
                stacking_interactions.append((paren[0] - 1, paren[1] - 1))
        for bracket in range(len(raw_sequence)):
            if bracket in open_paren:
                output = output + "("
            elif bracket in close_paren:
                output = output + ")"
            else:
                output = output + "."
        pairs = [(paren[0] - 1, paren[1] - 1) for paren in stem_pairs if paren[1] !=0]
        return (stacking_interactions, output, mfe, pairs)
    

    def mfe_fold_output(self):
        mfe = self.mfe_fold()
        mfe_output = mfe[1]
        mfe_energy = mfe[2]
        return [self.sequence, mfe_output, mfe_energy]

class SecondaryStructEnergy(object):

    def __init__(self, sequence, algorithm="Zuker"):
        self.sequence = sequence.upper()
        self.algorithm = algorithm.title()
        if self.algorithm == "Zuker":
            Zuker_sequence = ZukerPredictor(sequence=self.sequence)
            self.stacking_interactions, self.struct, self.free_energy, self.pairs = Zuker_sequence.mfe_fold()
        elif self.algorithm == "Nussinov":
            Nussinov_sequence = NussinovPredictor(sequence=self.sequence)
            self.struct = Nussinov_sequence.main_nussinov()

        self.output_string = fgb.BulgeGraph(dotbracket_str=self.struct)
        self.element_struc = fgb.BulgeGraph(dotbracket_str=self.struct).to_element_string()
        self.element_outputs = fgb.BulgeGraph(dotbracket_str=self.struct).to_element_string(with_numbers=True)
        self.output_domains = self.output_string.get_domains()
        self.output_strings = [i+j for i,j in zip(self.element_outputs[0], self.element_outputs[1])]
        self.hairpin_loops_indices = []
        self.multiloops_indices = []
        self.deltaG = 0
        self.deltaH = 0
        self.deltaS = 0

    def find_index(self, item, iterable):
        return [i for i,e in enumerate(iterable) if e == item]

    def find_pairs(self, index, pair_lst):
        for i in pair_lst:
            if i[0] == index:
                return (i[0], i[1])

    def multiloops(self):
        helix_basepair = []

        for key, value in self.output_domains.items():
            if key == "multiloops":
                for i in value:
                    temp_lst = []
                    for j in i:
                        if j in self.output_strings:
                            temp_lst.extend(self.find_index(j, self.output_strings))
                    self.multiloops_indices.append(temp_lst)

        true_indices = []

        for m_loop in self.multiloops_indices:
            temp_loop_lst = []
            for key, group in groupby(enumerate(m_loop), lambda x: x[0] - x[1]):
                collected = (map(itemgetter(1),group))
                collected = list(map(int, collected))
                if collected[0] != collected[-1]:
                    temp_loop_lst.append((collected[0], collected[-1]))
                elif collected[0] == collected[-1]:
                    temp_loop_lst.append(collected[0])

            true_indices.append(temp_loop_lst)

        acceptable_base_pairs = []

        for true_index in true_indices:
            temp_basepairs = []
            for i in true_index:
                if isinstance(i, int):
                    left_bp = i + 1
                    right_bp = i - 1
                    temp_basepairs.append(self.find_pairs(right_bp, self.pairs))
                    temp_basepairs.append(self.find_pairs(left_bp, self.pairs))
                elif isinstance(i, tuple) and len(i) > 1:
                    if i[0] < len(self.sequence):
                        left_bp = i[0] + 1
                        right_bp = i[0] - 1
                        temp_basepairs.append(self.find_pairs(left_bp, self.pairs))
                        temp_basepairs.append(self.find_pairs(right_bp, self.pairs))
                    elif i[1] < len(sequence) - 1:
                        left_bp = i[0] + 1
                        right_bp = i[1] - 1
                        temp_basepairs.append(self.find_pairs(left_bp, self.pairs))
                        temp_basepairs.append(self.find_pairs(right_bp, self.pairs))
            acceptable_base_pairs.append((temp_basepairs))

        seen_base_pairs = []
        unique_base_pairs = []
        for temp in acceptable_base_pairs:
            temp_unique = []
            for j in temp:
                if j is not None:
                    if j[0] not in seen_base_pairs:
                        temp_unique.append(j)
                        seen_base_pairs.append(j[1])
            unique_base_pairs.append(temp_unique)

        for unpaired_bp, helix in zip(self.multiloops_indices, unique_base_pairs):
            helix_basepair.append([len(unpaired_bp), len(helix)])

        return helix_basepair

    def hairpin_loops(self):
        loop_lengths = []

        for index, seq_id in enumerate(self.element_struc):
            if seq_id == "h":
                self.hairpin_loops_indices.append(index)

        for key, group in groupby(enumerate(self.hairpin_loops_indices), lambda x:x[0] - x[1]):
            collected = (map(itemgetter(1),group))
            collected = list(map(int, collected))
            loop_lengths.append(collected[-1] - collected[0] + 1)

        return loop_lengths

    def internal_bulge_loops (self):
        internal_bulges = []
        for index, loop in enumerate(self.output_strings):
            if loop[0] == "i" and loop not in internal_bulges:
                internal_bulges.append(loop)

        internal_loop_sizes, bulge_loop_sizes = [], []

        for loop in internal_bulges:
            size = self.output_string.get_bulge_dimensions(loop)
            print(size)


    def stacking_interactions(self):
        pass

    def free_energy_tm(self):
        pass


