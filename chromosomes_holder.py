import os
import pickle
import re

import numpy as np
from natsort import natsorted
import random


# TODO: where to add random seed

class AnnotationRecord:
    def __init__(self, chr_name, start, end, name, group, color=None):
        self.name = name
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.group = group
        self.color = color

    def __repr__(self):
        return f"{self.chr_name}_{self.name}, {self.start}, {self.end}, {self.group}"


class ChromosomesHolder:
    def __init__(self, species, root_path='./Data/species', use_cache=True):
        self.species = species
        self.root_path = root_path
        self.use_cache = use_cache

        self._chromosomes_path = {}
        self._fill_chromosomes_path()

        self.reverse_complement = {}
        self._fill_reverse_complement_info()

        self.cytobands = {}
        self._fill_cytobands_info()

        self._chromosome_sequence_cache = {}

    def get_all_chromosomes_name(self):
        return natsorted(list(self._chromosomes_path.keys()))

    def get_chromosome_sequence(self, chromosome_name):
        if chromosome_name in self._chromosome_sequence_cache:
            return self._chromosome_sequence_cache[chromosome_name]
        sequence = ""
        chr_path = self._chromosomes_path[chromosome_name]
        with open(chr_path) as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    pass
                else:
                    sequence += line
        if self.use_cache:
            self._chromosome_sequence_cache[chromosome_name] = sequence
        return sequence

    def get_segment(self, chromosome_name, start_of_sequence, sequence_length):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        assert start_of_sequence < len(chromosome_sequence), "Start of sequence is out of range."
        assert start_of_sequence + sequence_length < len(chromosome_sequence), "End of sequence is out of range."
        return chromosome_sequence[start_of_sequence:start_of_sequence + sequence_length]

    def get_random_segment(self, length, chromosome_name=None, remove_outlier=False, return_dict=False):
        if chromosome_name:
            chosen_chromosome_name = chromosome_name
        else:
            chosen_chromosome_name = random.choice(self.get_all_chromosomes_name())
        chosen_chromosome_sequence = self.get_chromosome_sequence(chosen_chromosome_name)

        random_start = random.randint(0, len(chosen_chromosome_sequence) - length)

        if remove_outlier:
            chosen_segment_cytoband = self.get_annotation_of_segment(chosen_chromosome_name, random_start, length)
            if chosen_segment_cytoband:
                while chosen_segment_cytoband.color in ['pi', 'pu']:
                    random_start = random.randint(0, len(chosen_chromosome_sequence) - length)
                    chosen_segment_cytoband = self.get_annotation_of_segment(chosen_chromosome_name, random_start,
                                                                             length)
            else:
                raise Exception("Remove outlier is only valid when cytoband annotation exists.")

        random_sequence = chosen_chromosome_sequence[random_start:random_start + length]
        if self.reverse_complement[chosen_chromosome_name]:
            random_sequence = self.get_reverse_complement_sequence(random_sequence)
        if return_dict:
            return {'chromosome_name': chosen_chromosome_name, 'sequence': random_sequence, 'start': random_start}
        else:
            return random_sequence

    def get_random_segments_list(self, length, num_segments, chromosome_name=None, overlap=False):
        if chromosome_name:
            chosen_chromosome_name = chromosome_name
        else:
            chosen_chromosome_name = random.choice(self.get_all_chromosomes_name())
        chosen_chromosome_sequence = self.get_chromosome_sequence(chosen_chromosome_name)

        retry_count = 0
        while num_segments * length > len(chosen_chromosome_sequence) and retry_count < 100:
            if chromosome_name:
                raise Exception("Number of segments times length is greater than chromosome")
            else:
                chosen_chromosome_name = random.choice(self.get_all_chromosomes_name())
                chosen_chromosome_sequence = self.get_chromosome_sequence(chosen_chromosome_name)
            retry_count += 1
        if retry_count == 100:
            raise ValueError("Sequence is too short to fit the desired number of segments of the specified length.")

        segments_list = []
        start_of_segments_list = []

        if overlap:
            for i in range(num_segments):
                random_start = random.randint(0, len(chosen_chromosome_sequence) - length)
                random_sequence = chosen_chromosome_sequence[random_start:random_start + length]
                if self.reverse_complement[chosen_chromosome_name]:
                    random_sequence = self.get_reverse_complement_sequence(random_sequence)
                segments_list.append(random_sequence)
                start_of_segments_list.append(random_start)
            return {'segments_sequences': segments_list, 'starts': start_of_segments_list,
                    'chromosome_name': chosen_chromosome_name}
        else:
            available_positions = list(range(len(chosen_chromosome_sequence) - length + 1))
            for i in range(num_segments):
                if not available_positions:
                    raise Exception("Not enough non-overlapping positions available")
                random_start = random.choice(available_positions)
                random_sequence = chosen_chromosome_sequence[random_start:random_start + length]
                if self.reverse_complement[chosen_chromosome_name]:
                    random_sequence = self.get_reverse_complement_sequence(random_sequence)
                segments_list.append(random_sequence)
                start_of_segments_list.append(random_start)
                # Remove positions that would overlap with the chosen segment
                available_positions = [pos for pos in available_positions if
                                       ((pos < random_start - length) and (pos - (random_start + length) > 0))
                                       or pos >= random_start + length]

            return {'segments_sequences': segments_list, 'starts': start_of_segments_list,
                    'chromosome_name': chosen_chromosome_name}

    def get_cytoband_segment(self, chromosome_name, cytoband_name):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        cytoband = self.cytobands[chromosome_name][cytoband_name]
        assert cytoband, "Cytoband not found."
        cytoband_segment = chromosome_sequence[cytoband.start:cytoband.end]
        if self.reverse_complement[chromosome_name]:
            cytoband_segment = self.get_reverse_complement_sequence(cytoband_segment)
        return cytoband_segment

    def get_annotation_of_segment(self, chromosome_name, start_of_sequence, sequence_length):
        for key, value in self.cytobands[chromosome_name].items():
            end_of_sequence = start_of_sequence + sequence_length
            if start_of_sequence <= value.end:
                if start_of_sequence >= value.start and end_of_sequence <= value.end:
                    return value
                if start_of_sequence >= value.start and end_of_sequence > value.end:
                    midpoint = (end_of_sequence - start_of_sequence) / 2
                    if start_of_sequence + midpoint <= value.end:
                        return value
                if start_of_sequence <= value.start and end_of_sequence < value.end:
                    midpoint = (end_of_sequence - start_of_sequence) / 2
                    if start_of_sequence + midpoint >= value.start:
                        return value
                if start_of_sequence <= value.start and end_of_sequence >= value.end:
                    return value

    def clear_cache(self):
        self._chromosome_sequence_cache = {}

    ''' Private Methods '''

    def _fill_chromosomes_path(self):
        chromosomes_path = os.path.join(self.root_path, self.species, 'chromosomes')
        for filename in os.listdir(chromosomes_path):
            if filename.lower().endswith('.fna'):
                chr_path = os.path.join(chromosomes_path, filename)
                chr_name = self._extract_chromosome_number(filename)
                self._chromosomes_path[chr_name] = chr_path

    def _fill_reverse_complement_info(self):
        info_path = os.path.join(self.root_path, self.species, 'extra', 'reverse_complement_info.pkl')
        rc_info = {}
        if os.path.exists(info_path):
            with open(info_path, 'rb') as handle:
                rc_info = pickle.load(handle)

        for chromosome_name in self.get_all_chromosomes_name():
            if chromosome_name in rc_info.keys():
                self.reverse_complement[chromosome_name] = rc_info[chromosome_name]
            else:
                self.reverse_complement[chromosome_name] = False

    def _fill_cytobands_info(self):
        for chromosome_name in self.get_all_chromosomes_name():
            self.cytobands[chromosome_name] = {}
        if self.species == "human":
            file_path = f"./Data/species/{self.species}/bedfiles/chm13v2.0_telomere.bed"
            switch = 0
            with open(file_path) as file:
                for line in file:
                    line = line.strip()
                    parts = line.split("\t")
                    chr_name = parts[0][3:]
                    self.cytobands[chr_name][f"tel_{switch + 1}"] = AnnotationRecord(chr_name=parts[0],
                                                                                     start=int(parts[1]),
                                                                                     end=int(parts[2]),
                                                                                     name=f"tel_{switch + 1}",
                                                                                     group=f"telomere",
                                                                                     color=None)

                    switch = 1 - switch

            file_path = f"./Data/species/{self.species}/bedfiles/chm13v2.0_cytobands_allchrs_color.bed"
            with open(file_path) as file:
                for line in file:
                    line = line.strip()
                    parts = line.split("\t")
                    chr_name = parts[0][3:]
                    self.cytobands[chr_name][f"{parts[3]}"] = AnnotationRecord(chr_name=parts[0],
                                                                               start=int(parts[1]),
                                                                               end=int(parts[2]),
                                                                               name=f"{parts[3]}",
                                                                               group="cytoband",
                                                                               color=parts[4])

    @staticmethod
    def _extract_chromosome_number(string):
        match = re.search(r'(chromosome\s*|scaffold_)(\d+|[A-Za-z]+)', string, re.IGNORECASE)
        if match:
            return match.group(2)
        return None

    @staticmethod
    def get_reverse_complement_sequence(sequence):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        bases = [complement[base] for base in sequence]
        bases = reversed(bases)
        return ''.join(bases)

    @staticmethod
    def split_chromosome_fasta_files(whole_genome_file_path, output_directory):
        with open(whole_genome_file_path) as fasta_file:
            f = None
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if not (f is None):
                        f.close()
                    filename = line[1:].replace('.', "").replace("/", "") + '.fna'
                    print(filename)
                    fasta_file_path = os.path.join(output_directory, filename)
                    f = open(fasta_file_path, "w")
                    f.write(line + "\n")
                else:
                    line = line.upper()
                    f.write(line)
            f.close()


if __name__ == '__main__':
    genome = ChromosomesHolder("human")
    genome.get_random_segment(length=100, remove_outlier=True)
    print("")
