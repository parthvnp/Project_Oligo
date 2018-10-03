import numpy as np


def complement(sequence):
    ACCEPTABLE_BASES = {"A": "T", "G": "C", "C":"G", "T": "A"}
    complementSequence = ""
    for base in sequence:
        complementSequence += ACCEPTABLE_BASES[base]
    return complementSequence

def primer_fragmentation(sequence, sequence_length, max_length, overlap,):
    if sequence_length % (max_length - overlap) == 0:
        minFrags = sequence_length//(max_length-overlap)
    else:
        minFrags = (sequence_length//(max_length - overlap)) + 1 
        print(minFrags)
    primersIndex = np.zeros(shape=(minFrags, 3), dtype=np.int16)
    sequenceStart = 0
    sequenceStop  = max_length
    minAcceptableDistance = 5
    assert max_length - 2 * overlap > minAcceptableDistance
    for i in range(minFrags):
        if i % 2 == 0:
            primersIndex[i] = [sequenceStart, sequenceStop, 1]
            sequenceStart = sequenceStop - overlap
            sequenceStop = sequenceStart + max_length
        else:
            primersIndex[i] = [sequenceStart, sequenceStop, -1]
            sequenceStart = sequenceStop - overlap
            sequenceStop = sequenceStart + max_length
    forwardSequence = ""
    reverseSequence = ""
    complementSequence = complement(sequence)
    for i in range(primersIndex.shape[0]):
        if i % 2 == 0:
            print(sequence[primersIndex[i][0]:primersIndex[i][1]])
        if i % 2 != 0:
            print(complementSequence[primersIndex[i][0]:primersIndex[i][1]])
    
    primersIndex[primersIndex > sequence_length] = sequence_length - 1 # to make sure sequence length is within the acceptable range
        
    return primersIndex

randSequence = "AGACGAACAAATTACTAGAGTGCCGCTTTCAGCCCCCCTGTCGTCGCCGACGTCTGTAATATGGCGTTGTTGTGATTCGACTCTATTGAGGCATCAACTGATGCGTAAGGAGATCTGGAATGAATTGGCCTATGTCACTGAAACTGTCCAAACACCCAATGTCGTTAGTGAAGGTTCTGACGCATACCTCCTTCGTTGAG"
print(primer_fragmentation(randSequence, len(randSequence), 40, 10))
print(len(randSequence))

            
