import forgi.graph.bulge_graph as fgb
from itertools import groupby
from operator import itemgetter
from Zuker import RNAstructure
import numpy as np

sequence = "GCATGTTGGAGTTATTGCCAACAGCGCAGCAGCAGCTAGCGCTACGCGGCATGCGCGCGATCGGCGCATCCTTAAGCGCACAGTCGCGCATAGCGGCGCATGCATACGATC"

p = RNAstructure.RNA.fromString(sequence, backbone="DNA")
p.FoldSingleStrand()
p.WriteDotBracket("foo.dbn")
x = p.GetFreeEnergy(1)
print(x)
stem_pairs = [(i,p.GetPair(i)) for i in p.iterIndices()]
output = ""
open_paren, close_paren, stacking_interactions = [], [], []
for paren in stem_pairs:
    if paren[1] != 0 and paren[1] > paren[0]:
            open_paren.append(paren[0]-1)
            close_paren.append(paren[1]-1)
            stacking_interactions.append((paren[0]- 1, paren[1] - 1))

for bracket in range(len(p)):

    if bracket in open_paren:
        output = output + "("
    elif bracket in close_paren:
        output = output + ")"
    else:
        output = output + "."

pairs = [(paren[0] - 1, paren[1] -1) for paren in stem_pairs if paren[1] != 0]
print(f"{output}\n{sequence}")
# Structure Decomposition

output_string = fgb.BulgeGraph(dotbracket_str=output)
output_element = fgb.BulgeGraph(dotbracket_str=output).to_element_string()
output_elements = fgb.BulgeGraph(dotbracket_str=output).to_element_string(with_numbers=True)
ext_loops = output_string.find_external_loops()
hairpins = []
hairpin_indices = []
multiloops = []
multiloops_indices = []
internal_loops = []
internal_loops_indices = []
five_primes = []
five_primes_indices = []
three_primes = []
three_primes_indices = []

for o,j in enumerate(output_element):
    if j == "h":
        hairpins.append(o)
    if j == "i":
        internal_loops.append(o)
    if j == "m":
        multiloops.append(o)
    if j == "f":
        five_primes.append(o)
    if j == "t":
        three_primes.append(o)

for k,g in groupby(enumerate(hairpins), lambda x: x[0] - x[1]):
    group = (map(itemgetter(1),g))
    group = list(map(int, group))
    hairpin_indices.append((group[0], group[-1]))

for k,g in groupby(enumerate(internal_loops), lambda x: x[0] - x[1]):
    group = (map(itemgetter(1),g))
    group = list(map(int, group))
    internal_loops_indices.append((group[0], group[-1]))

output_domains  = output_string.get_domains()

print(f"{output_elements}")

output_strings = [i+j for i,j in zip(output_elements[0], output_elements[1])]
for j,i in enumerate(output_strings):
    if i[0] == "m":
        print(f"{j} {i[0]}:{i[1]}")
print(f"output_strings: {output_strings}")

# MultiLoop Decomposition

def find_index(item, iterable):
    return [i for i,e in enumerate(iterable) if e == item]
for key,value in output_domains.items():
    if key == "multiloops":
        for i in value:
            y = []
            for j in i:
                if j in output_strings:
                    y.extend(find_index(j, output_strings))
            multiloops_indices.append(y)

print(output_domains)
print(f"Multiloops: {multiloops_indices}")

# Group the indices so that we can figure out the pairs that are part of M-loop
true_indices = []

for m in multiloops_indices:
    v = []
    for k,g in groupby(enumerate(m), lambda x: x[0] - x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int, group))
        if group[0] != group[-1]:
            v.append((group[0], group[-1]))
        elif group[0] == group[-1]:
            v.append(group[0])
    true_indices.append(v)

def find_pairs(index,pair_lst):
    for i in pair_lst:
        if i[0] == index:
            return (i[0], i[1])

print(f"True_Indices: {true_indices}")
temp = []
for t in true_indices:
    base_pairs = []
    for i in t:
        if isinstance(i, int):
            left = i + 1
            right = i - 1
            base_pairs.append(find_pairs(right, pairs))
            base_pairs.append(find_pairs(left, pairs))
        if isinstance(i, tuple) and len(i) > 1:
            if i[0] < len(sequence):
                left = i[0] + 1
                right = i[0] - 1
                base_pairs.append(find_pairs(left, pairs))
                base_pairs.append(find_pairs(right, pairs))
            elif i[1] < len(sequence) - 1:
                left = i[0] + 1
                right = i[0] - 1
                base_pairs.append(find_pairs(left, pairs))
                base_pairs.append(find_pairs(right, pairs))
    temp.append((base_pairs))

# print(f"\n\nBase Pairs: {temp}")
# Make sure the base pairs are not repeating
seen_values = []
unique_base_pairs = []
for t in temp:
    temp_unique = []
    for j in t:
        if j:
            if j[0] not in seen_values:
                temp_unique.append(j)
                seen_values.append(j[1])
    unique_base_pairs.append((temp_unique))

# print(f"Unique base pairs: {unique_base_pairs}")

# Calculate the energy of the multiloop in kcal/mol
deltaG = 0
deltaH = 0
deltaS = 0

for base_pairs, helix  in zip(multiloops_indices, unique_base_pairs):
    temp_deltaG, temp_deltaH = 0, 0
    temp_deltaG = 3.0 + (0.2 * len(base_pairs)) + (0.2 * len(helix))
    temp_deltaH = 9.0
    print(f"{len(base_pairs)}: {len(helix)}")
    deltaG = deltaG + temp_deltaG
    deltaH = deltaH + temp_deltaH
    print(F"dG: {deltaG} dH: {deltaH}")

# Hairpin loop Deconstruction

print("Hairpin  Decomposition")
print(hairpins)

ranges = []
length_loop = []

for k,g in groupby(enumerate(hairpins),lambda x:x[0]-x[1]):
    range_temp = []
    group = (map(itemgetter(1),g))
    group = list(map(int,group))
    range_temp.append((group[0],group[-1]))
    length_loop.append(group[-1] - group[0] + 1)
    ranges.append(range_temp)
print(length_loop)

# Internal Loop Decomposition
ip_id = []
for index, loop in enumerate(output_strings):
    if loop[0] == "i":
        if loop not in ip_id:
            ip_id.append(loop)
        # print(f"{index}: {loop}")
        # print(output_string.get_bulge_dimensions(loop))
        # print(ip_id)
# Collect the energy value of base pairs and determine the number of unpaired base pairs
# The columns correspond to number of residues in SS and the and rows
# correspond to number of helices. The rows start from three helices. These values were take from the following paper :
# The Thermodynamics of DNA structural motifs
#Vol. 33:415-440 (Volume publication date 9 June 2004)

multiloop_energy = [[2.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
                            [-1.0, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6],
                            [2.0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
                            [2.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2],
                            [2.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]]
