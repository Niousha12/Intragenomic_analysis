import collections

import math
import numpy as np
from tqdm import tqdm


class CGR:
    def __init__(self, seq, k_mer):
        self.data = seq
        self.k_mer = k_mer

    def count_kmers(self):
        kmer_dict = collections.defaultdict(int)
        for i in range(len(self.data) - (self.k_mer - 1)):
            kmer_dict[self.data[i:i + self.k_mer]] += 1

        for kmer_key in list(kmer_dict.keys()):
            if "N" in kmer_key:
                del kmer_dict[kmer_key]
        return kmer_dict

    def probabilities(self):  # k-mers count / length of the whole dna
        probabilities = collections.defaultdict(float)
        if self.k_mer == 0:
            probabilities[self.data] = 1.0
        else:
            kmer_dict = self.count_kmers()
            kmer_dict_sum = sum(kmer_dict.values())
            for key, value in kmer_dict.items():
                probabilities[key] = float(value) / kmer_dict_sum  # / (len(self.data) - self.k_mer + 1)
        return probabilities

    def get_fcgr(self):
        probabilities = self.count_kmers()
        array_size = int(math.sqrt(4 ** self.k_mer))
        chaos = []
        for _ in range(array_size):
            chaos.append([0] * array_size)
        max_x = array_size
        max_y = array_size
        pos_x = 1
        pos_y = 1
        for key, value in probabilities.items():
            for char in reversed(key):
                if char == "G" or char == "g":  # T
                    pos_x += max_x / 2
                elif char == "A" or char == "a":  # C
                    pos_y += max_y / 2
                elif char == "T" or char == "t":  # G
                    pos_x += max_x / 2
                    pos_y += max_y / 2
                max_x /= 2
                max_y /= 2

            chaos[int(pos_y - 1)][int(pos_x - 1)] = value
            max_x = array_size
            max_y = array_size
            pos_x = 1
            pos_y = 1

        c = np.array(chaos)
        return c

    def get_cgr(self):
        seq = list(self.data)
        d = {'C': (0, 0), 'A': (0, 1), 'T': (1, 1), 'G': (1, 0)}
        size = 2 ** self.k_mer
        mat = np.zeros((2 ** self.k_mer, 2 ** self.k_mer))
        x, y = 0.5, 0.5

        for nucleotide in tqdm(seq):
            if nucleotide in d:
                x += 0.5 * (d[nucleotide][0] - x)
                y += 0.5 * (d[nucleotide][1] - y)
                mat[int(y * size) - 1][int(x * size) - 1] = 1
        return mat

    def get_scatter_cgr(self):
        seq = list(self.data)
        d = {'A': (0, 0), 'C': (0, 1), 'G': (1, 1), 'T': (1, 0)}
        x, y = 0.5, 0.5
        X, Y = [], []
        i = 0

        while seq:
            nucleotide = seq.pop(0)
            if nucleotide in d:
                x += 0.5 * (d[nucleotide][0] - x)
                y += 0.5 * (d[nucleotide][1] - y)
                X.append(x)
                Y.append(y)
                i += 1
        return X, Y

    @staticmethod
    def normalize(input_matrix):
        matrix = input_matrix - np.min(input_matrix)
        matrix = matrix / np.max(matrix)
        return matrix

    @staticmethod
    def array2img(array, bits=8, resolution=4):
        m, M = array.min(), array.max()
        # rescale to [0,1]
        img_rescaled = (array - m) / (M - m)

        max_color = resolution ** bits - 1

        # invert colors black->white
        img_array = np.ceil(max_color - img_rescaled * max_color)
        dtype = eval(f"np.int{bits}")
        img_array = np.array(img_array, dtype=dtype)

        return img_array
