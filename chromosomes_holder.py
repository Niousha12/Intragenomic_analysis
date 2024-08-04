import os
import pickle
import re

import numpy as np
from matplotlib import pyplot as plt
from natsort import natsorted
import random

from tqdm import tqdm

from PIL import Image
from chaos_game_representation import CGR

# TODO: where to add random seed
random.seed(46)
np.random.seed(46)


class AnnotationRecord:
    def __init__(self, chr_name, start, end, name, group, color=None):
        self.name = name
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.group = group
        self.color = color
        self.color_name = None
        self.display_color = None
        color_name = {"w": "White", "lg": "Light Gray", "g": "Gray", "dg": "Dark Gray", "b": "Black", "pi": "Pink",
                      "pu": "Purple", "R": "Representative"}

        if self.color is not None:
            self.color_name = color_name[self.color]

    def __repr__(self):
        return f"{self.chr_name}_{self.name}, {self.start}, {self.end}, {self.group}"

    @staticmethod
    def get_color_display_dict():
        return {"White": "#C0C0C0", "Light Gray": "#808080", "Gray": "#696969", "Dark Gray": "#505050",
                "Black": "#000000", "Pink": "#fc9ea3", "Purple": "#c89efc", "Representative": "#ff0000"}


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
        # TODO: does not work on bacteria and Archaea
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

    def get_all_chromosome_length(self):
        return {chromosome_name: len(self.get_chromosome_sequence(chromosome_name)) for chromosome_name in
                self.get_all_chromosomes_name()}

    def get_largest_chromosome_length(self):
        return max(self.get_all_chromosome_length().values())

    def get_appropriate_segment_length(self, scale='genome'):
        """
        :param scale: chromosome or genome
        """
        if scale == 'genome':
            total_genome_length = sum(self.get_all_chromosome_length().values())
            app_length = total_genome_length // 6234
            if app_length > 1000:
                app_length //= 1000
                app_length *= 1000
        else:
            app_length = self.get_largest_chromosome_length() // 500
            if app_length > 1000:
                app_length //= 1000
                app_length *= 1000
        return app_length

    def get_segment(self, chromosome_name, start_of_segment, sequence_length):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        assert start_of_segment > len(chromosome_sequence), "Start of segment is out of range."
        assert start_of_segment + sequence_length > len(chromosome_sequence), "End of segment is out of range."
        return chromosome_sequence[start_of_segment:start_of_segment + sequence_length]

    def get_random_segment(self, length, chromosome_name=None, remove_outlier=False, return_dict=False):
        # Find the chromosome name and its sequence
        if chromosome_name:
            chosen_chromosome_name = chromosome_name
            determined_chromosome = True
        else:
            chosen_chromosome_name = random.choice(self.get_all_chromosomes_name())
            determined_chromosome = False
        chosen_chromosome_sequence = self.get_chromosome_sequence(chosen_chromosome_name)

        # This condition only happens when the sequence is too short to fit the desired length, does not happen in human
        if determined_chromosome and len(chosen_chromosome_sequence) < (length // 2):
            raise ValueError("Sequence is too short to fit the desired length.")

        while len(chosen_chromosome_sequence) < (length // 2):
            chosen_chromosome_name = random.choice(self.get_all_chromosomes_name())
            chosen_chromosome_sequence = self.get_chromosome_sequence(chosen_chromosome_name)
        if len(chosen_chromosome_sequence) < length:
            chosen_chromosome_sequence += 'N' * (length - len(chosen_chromosome_sequence))

        # Find the random start position where it doesn't return a sequence with all N
        random_start = -1

        processed_chosen_seq = ""
        while len(processed_chosen_seq) < (length // 2):
            random_start = random.randint(0, len(chosen_chromosome_sequence) - length)
            random_sequence = chosen_chromosome_sequence[random_start:random_start + length]
            processed_chosen_seq = random_sequence.replace("N", "")

            if remove_outlier:
                chosen_segment_cytoband = self.get_annotation_of_segment(chosen_chromosome_name, random_start, length,
                                                                         group="cytoband")
                if chosen_segment_cytoband:
                    if chosen_segment_cytoband.color in ['pi', 'pu']:
                        processed_chosen_seq = ""  # invalidate the sequence
                else:
                    raise Exception("Remove outlier is only valid when cytoband annotation exists.")

        random_sequence = chosen_chromosome_sequence[random_start:random_start + length]
        if self.reverse_complement[chosen_chromosome_name]:
            random_sequence = self.get_reverse_complement(random_sequence)

        if return_dict:
            return {'chromosome_name': chosen_chromosome_name, 'sequence': random_sequence, 'start': random_start}
        else:
            return random_sequence

    def choose_random_fragment_index(self, chromosome_name, length, outlier=True):
        candidates = []
        for i in range(len(self.get_chromosome_sequence(chromosome_name)) // length):
            color = self.get_annotation_of_segment(chromosome_name, i * length, length, group='cytoband').color
            if outlier:
                if color in ['pi', 'pu']:
                    candidates.append(i)
            else:
                if color not in ['pi', 'pu']:
                    candidates.append(i)
        return np.random.choice(candidates)

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

        if overlap:
            start_of_segments_list = []
            for i in range(num_segments):
                random_start = random.randint(0, len(chosen_chromosome_sequence) - length)
                random_segment = chosen_chromosome_sequence[random_start:random_start + length]
                if self.reverse_complement[chosen_chromosome_name]:
                    random_segment = self.get_reverse_complement(random_segment)
                segments_list.append(random_segment)
                start_of_segments_list.append(random_start)
            return {'segments_sequences': segments_list, 'starts': start_of_segments_list,
                    'chromosome_name': chosen_chromosome_name}
        else:
            start_of_segments_list = set()
            while len(segments_list) < num_segments:
                start_position = random.randint(0, len(chosen_chromosome_sequence) - length)

                found_overlap = False
                for pos in start_of_segments_list:
                    if pos < start_position < pos + length or start_position < pos < start_position + length:
                        found_overlap = True
                        break

                if not found_overlap:
                    segments_list.append(chosen_chromosome_sequence[start_position:start_position + length])
                    start_of_segments_list.add(start_position)

            return {'segments_sequences': segments_list, 'starts': start_of_segments_list,
                    'chromosome_name': chosen_chromosome_name}

    def get_cytoband_segment(self, chromosome_name, cytoband_name):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        cytoband = self.cytobands[chromosome_name][cytoband_name]
        assert cytoband, "Cytoband not found."
        cytoband_segment = chromosome_sequence[cytoband.start:cytoband.end]
        if self.reverse_complement[chromosome_name]:
            cytoband_segment = self.get_reverse_complement(cytoband_segment)
        return cytoband_segment

    def get_annotation_of_segment(self, chromosome_name, start_of_segment, segment_length, group=None):
        for key, value in self.cytobands[chromosome_name].items():
            if group and value.group != group:
                continue
            end_of_segment = start_of_segment + segment_length
            if start_of_segment <= value.end:
                if start_of_segment >= value.start and end_of_segment <= value.end:
                    return value
                if start_of_segment >= value.start and end_of_segment > value.end:
                    midpoint = (end_of_segment - start_of_segment) / 2
                    if start_of_segment + midpoint <= value.end:
                        return value
                if start_of_segment <= value.start and end_of_segment < value.end:
                    midpoint = (end_of_segment - start_of_segment) / 2
                    if start_of_segment + midpoint >= value.start:
                        return value
                if start_of_segment <= value.start and end_of_segment >= value.end:
                    return value

    def get_chromosome_non_overlapping_segments(self, chromosome_name, segments_length):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)

        segments_sequences = []
        segments_information = []
        step_length = np.floor(len(chromosome_sequence) / segments_length)
        for i in tqdm(range(int(step_length) - 1)):
            start_of_segment = i * segments_length
            end_of_segment = (i + 1) * segments_length
            segment_sequence = chromosome_sequence[start_of_segment:end_of_segment]
            if self.reverse_complement[chromosome_name]:
                segment_sequence = self.get_reverse_complement(segment_sequence)
            segments_sequences.append(segment_sequence)

            segment_information = self.get_annotation_of_segment(chromosome_name, start_of_segment, segments_length,
                                                                 group="cytoband")
            segments_information.append(segment_information)
        return {'segments_sequences': segments_sequences, 'segments_information': segments_information}

    def plot_fcgr(self, chromosome_name, start_of_segment=None, segment_length=None, k_mer=9):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        if start_of_segment is None:
            start_of_segment = 0
        if segment_length is None:
            segment_length = len(chromosome_sequence)
        segment_sequence = chromosome_sequence[start_of_segment:start_of_segment + segment_length]
        fcgr = CGR(segment_sequence, k_mer).get_fcgr()
        fcgr_image = CGR.array2img(fcgr, bits=8)
        fcgr_pil = Image.fromarray(fcgr_image, 'L')

        plt.imshow(fcgr_pil, cmap="gray")
        plt.xticks([])  # Remove x ticks
        plt.yticks([])  # Remove y ticks

        # Coordinates of points to label (example points)
        points = [(-12, 520), (-15, -6), (522, -6), (525, 520)]  # List of (x, y) tuples
        labels = ["A", "C", "G", "T"]  # Labels for each point

        # Annotating each point
        for (x, y), label in zip(points, labels):
            plt.text(x, y, label, color='black', fontsize=14, ha='center', va='center')
        plt.title(f"Chromosome {chromosome_name}", fontsize=14)

        save_path = os.path.join('Figures', 'FCGRs', self.species)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        plt.savefig(
            f"{save_path}/chr_{chromosome_name}_range_{start_of_segment}_{start_of_segment + segment_length}.png",
            dpi=300, bbox_inches='tight', transparent=True)
        # plt.show()

    def clear_cache(self):
        self._chromosome_sequence_cache = {}

    '''Run this method to create chromosome files from a whole genome file (Run only once)'''

    def create_chromosomes_files(self, whole_genome_path):
        file_path = os.path.join(self.root_path, self.species, 'extra', whole_genome_path)
        save_path = os.path.join(self.root_path, self.species, 'chromosomes')
        with open(file_path) as fasta_file:
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
                    path = os.path.join(save_path, filename)
                    f = open(path, "w")
                    f.write(line + "\n")
                else:
                    line = line.upper()
                    f.write(line)
            f.close()

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
        # TODO: Add centromere annotations
        for chromosome_name in self.get_all_chromosomes_name():
            self.cytobands[chromosome_name] = {}
        if self.species == "Human":
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
    def get_reverse_complement(sequence):
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
    genome = ChromosomesHolder("Human")

    # genome.get_random_segment(1000, remove_outlier=True)
    # genome.get_chromosome_non_overlapping_segments("1", 1000)
    genome.plot_fcgr("Y")
    # genome.create_chromosomes_files("GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic.fna")
    print("")
