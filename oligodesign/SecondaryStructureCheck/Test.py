from Zuker import RNAstructure
import forgi.graph.bulge_graph as fgb
from itertools import groupby
from operator import itemgetter

def add_bulge(bulges, bulge, context, message):
    # bulge = (context, bulge)
    bulges[context] = bulges.get(context, []) + [bulge]
    return bulges
def condense_stem_pairs(stem_pairs):

    stem_pairs.sort()

    prev_pair = (-10, -10)

    stems = []
    start_pair = None

    for pair in stem_pairs:
        # There's a potential bug here since we don't check the direction
        # but hopefully it won't bite us in the ass later
        if abs(pair[0] - prev_pair[0]) != 1 or abs(pair[1] - prev_pair[1]) != 1:
            if start_pair is not None:
                stems += [(start_pair, prev_pair)]
            start_pair = pair

        prev_pair = pair

    if start_pair is not None:
        stems += [(start_pair, prev_pair)]

    return stems
def find_bulges_and_stems(brackets):
    prev = 'x'
    context = 0

    bulges = dict()
    finished_bulges = []
    context_depths = dict()

    opens = []
    stem_pairs = []

    dots_start = 0
    context_depths[0] = 0

    i = 0
    for i in range(len(brackets)):
        if brackets[i] == '(':
            opens.append(i)

            if prev == '(':
                context_depths[context] = context_depths.get(context, 0) + 1
                continue
            else:
                context += 1
                context_depths[context] = 1

            if prev == '.':
                dots_end = i - 1
                bulges = add_bulge(bulges, (dots_start, dots_end), context, "4")

        if brackets[i] == ')':
            if len(opens) == 0:
                raise Exception("Unmatched close bracket")

            stem_pairs.append((opens.pop(), i))

            context_depths[context] -= 1

            if context_depths[context] == 0:
                if context in bulges:
                    finished_bulges += bulges[context]
                bulges[context] = []
                context -= 1

            if prev == '.':
                dots_end = i - 1
                bulges = add_bulge(bulges, (dots_start, dots_end), context, "2")

        if brackets[i] == '.':
            if prev == '.':
                continue

            dots_start = i

        prev = brackets[i]

    if prev == '.':
        dots_end = i
        bulges = add_bulge(bulges, (dots_start, dots_end), context, "7")
    elif prev == '(':
        print ("Unmatched bracket at the end", file=sys.stderr)
        sys.exit(1)
    """
    elif prev == ')':
        bulges = add_bulge(bulges, (i+1, i+1), context, "8")
    """

    if context in bulges.keys():
        finished_bulges += bulges[context]

    if len(opens) > 0:
        raise Exception("Unmatched open bracket")

    stem_pairs.sort()
    stems = condense_stem_pairs(stem_pairs)

    return finished_bulges, stems

struct = "..........(((....)))....((.(((.....))).))...((((.((...)).))))."
struct_2 = find_bulges_and_stems(struct)
# strcut_3 = condense_stem_pairs()
print(struct_2)
str_elements = fgb.BulgeGraph(dotbracket_str=struct)
x1 = str_elements.to_element_string(with_numbers=True)
x3 = str_elements.find_external_loops()
x4 = str_elements.describe_multiloop(x3)
lst = []

for o,j in enumerate(x1):
    if j == "h":
        lst.append(o)
print(f"{x1}\nhairpins: {lst}\nmultiloops: {x3}\nMultiloops-description:{x4}")

ranges = []

for k,g in groupby(enumerate(lst), lambda x: x[0] - x[1]):
    group = (map(itemgetter(1),g))
    group = list(map(int, group))
    ranges.append((group[0], group[-1]))

print(ranges)

output_strings = ['s0', 's0', 'm0', 's1', 's1', 's1', 's1', 's1', 's1', 'h0', 'h0', 'h0', 'h0', 'h0', 'h0', 'h0', 'h0', 's1', 's1', 's1', 's1', 's1', 's1', 's2', 's2', 'h1', 'h1', 'h1', 's2', 's2', 'm2', 's0', 's0', 'm3', 's3', 's3', 'h2', 'h2', 'h2', 'h2', 'h2', 'h2', 'h2', 's3', 's3', 't0', 't0', 't0']

multiloops = [['m0', 'm1', 'm2'], ['m3', 't0']]
super_indices = []
for i in multiloops:
    indices = []
    for j in i:
        if j in output_strings:
            indices.append(output_strings.index(j))
    super_indices.append(indices)
print(super_indices)
