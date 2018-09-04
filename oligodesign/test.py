import re
def consecutive_basepair_check(sequence):
    possible_pairs = [match[0] for match in re.findall(r'((\w)\2{3,})', sequence)]
    return possible_pairs

print(consecutive_basepair_check("GGGGTGCTGCATGCAAAAA"))
