def random_sequnece_generation(number_sequence=40, length_increment=20):
    random_length = [i * length_increment for i in range(1, number_sequence + 1)]
    acceptable_bases = ("A", "T", "G", "C")
    random_sequence_collection = []
    for i in random_length:
        random_sequence = ""
        for j in range(i):
            random_sequence = random_sequence + choice(acceptable_bases)
        random_sequence_collection.append(random_sequence)
    return random_sequence_collection
