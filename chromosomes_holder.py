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
from constants import RESOLUTION_DICT, GENOME_LENGTH, BITS_DICT
import pandas as pd
import plotly.graph_objects as go

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
                      "pu": "Purple", "R": "RS"}

        if self.color is not None:
            self.color_name = color_name[self.color]

    def __repr__(self):
        return f"{self.chr_name}_{self.name}, {self.start}, {self.end}, {self.group}"

    @staticmethod
    def get_color_display_dict():
        return {"White": "#C0C0C0", "Light Gray": "#808080", "Gray": "#696969", "Dark Gray": "#505050",
                "Black": "#000000", "Pink": "#ffb3d9", "Purple": "#c89efc", "RS": "#ff0000"}


class ChromosomesHolder:
    def __init__(self, species, root_path='Data', use_cache=True):
        self.species = species
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        self.root_path = os.path.join(project_root, root_path)
        self.use_cache = use_cache

        self._chromosomes_path = {}
        self._fill_chromosomes_path()

        self.reverse_complement = {}
        self._fill_reverse_complement_info()

        self._chromosome_sequence_cache = {}
        # if species is not in the genome_length list, find its genome_length
        if species not in GENOME_LENGTH.keys():
            self.genome_length = self._get_genome_length()
        else:
            self.genome_length = GENOME_LENGTH[species] * 1_000_000  # self._get_genome_length()

        self.cytobands = {}
        self._fill_cytobands_info()

    def get_all_chromosomes_name(self, include_whole_genome=False):
        # if self.genome_length is None:
        #     self.genome_length = self._get_genome_length()

        chr_name_list = natsorted(list(self._chromosomes_path.keys()))
        if include_whole_genome:
            if self.genome_length < 1e9:
                chr_name_list.append("Whole Genome")
        return chr_name_list

    def _get_genome_length(self):
        return sum([len(self.get_chromosome_sequence(chromosome_name)) for chromosome_name in
                    natsorted(list(self._chromosomes_path.keys()))])

    def get_chromosome_sequence(self, chromosome_name, apply_reverse_complement=False):
        if chromosome_name == "Whole Genome":
            all_chromosomes_sequence = ""
            for chromosome_name in self.get_all_chromosomes_name():
                all_chromosomes_sequence += self.get_chromosome_sequence(chromosome_name) + "N"
            return all_chromosomes_sequence

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
        if apply_reverse_complement and self.reverse_complement[chromosome_name]:
            sequence = self.get_reverse_complement(sequence)
        if self.use_cache:
            self._chromosome_sequence_cache[chromosome_name] = sequence
        return sequence

    def get_all_chromosome_length(self, include_whole_genome=False):
        return {chromosome_name: len(self.get_chromosome_sequence(chromosome_name)) for chromosome_name in
                self.get_all_chromosomes_name(include_whole_genome=include_whole_genome)}

    def get_largest_chromosome_length(self, include_whole_genome=False):
        return max(self.get_all_chromosome_length(include_whole_genome).values())

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
        assert start_of_segment < len(chromosome_sequence), "Start of segment is out of range."
        assert start_of_segment + sequence_length <= len(chromosome_sequence), "End of segment is out of range."
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

        if return_dict:
            return {'chromosome_name': chosen_chromosome_name, 'sequence': random_sequence, 'start': random_start}
        else:
            return random_sequence

    @staticmethod
    def get_random_segment_from_any(length, sequence=None, return_dict=False):
        # This condition only happens when the sequence is too short to fit the desired length
        if len(sequence) < (length // 2):
            raise ValueError("Sequence is too short to fit the desired length.")

        if len(sequence) < length:
            sequence += 'N' * (length - len(sequence))

        # Find the random start position where it doesn't return a sequence with all N
        random_start = -1

        processed_chosen_seq = ""
        while len(processed_chosen_seq) < (length // 2):
            random_start = random.randint(0, len(sequence) - length)
            random_sequence = sequence[random_start:random_start + length]
            processed_chosen_seq = random_sequence.replace("N", "")

        random_sequence = sequence[random_start:random_start + length]

        if return_dict:
            return {'chromosome_name': 'Custom', 'sequence': random_sequence, 'start': random_start}
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
                if overlap:
                    print("Number of segments times length is greater than chromosome")
                    break
                else:
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
        return chromosome_sequence[cytoband.start:cytoband.end]

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
            segments_sequences.append(segment_sequence)

            segment_information = self.get_annotation_of_segment(chromosome_name, start_of_segment, segments_length,
                                                                 group="cytoband")
            segments_information.append(segment_information)
        return {'segments_sequences': segments_sequences, 'segments_information': segments_information}

    def get_chromosome_fcgr_min_max(self, chromosome_name, start_of_segment=None, segment_length=None, k_mer=9):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        if start_of_segment is None:
            start_of_segment = 0
        else:
            if (start_of_segment < 0) or (start_of_segment >= len(chromosome_sequence)):
                raise ValueError(f"Start of segment {start_of_segment} is out of range of the chromosome length "
                                 f"{len(self.get_chromosome_sequence(chromosome_name))}")
        if segment_length is None:
            segment_length = len(chromosome_sequence)
        else:
            if (segment_length <= 0) or (start_of_segment + segment_length > len(chromosome_sequence)):
                raise ValueError(f"Segment length {segment_length} is out of range of the chromosome length "
                                 f"{len(self.get_chromosome_sequence(chromosome_name))}")
        segment_sequence = chromosome_sequence[start_of_segment:start_of_segment + segment_length]

        fcgr = CGR(segment_sequence, k_mer).get_fcgr()
        m, M = fcgr.min(), fcgr.max()
        return m, M

    def get_global_min_max(self, start_of_segment=None, segment_length=None, k_mer=9):
        global_min = np.inf
        global_max = 0
        for chr_name in tqdm(self.get_all_chromosomes_name()):
            m, M = self.get_chromosome_fcgr_min_max(chr_name, start_of_segment, segment_length, k_mer)
            if m < global_min:
                global_min = m
            if M > global_max:
                global_max = M
        return global_min, global_max

    def get_max_kmer_values(self, chromosome_name, start_of_segment=None, segment_length=None, k_mer=6, top_n=50,
                            bottom_n=50):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        if start_of_segment is None:
            start_of_segment = 0
        else:
            if (start_of_segment < 0) or (start_of_segment >= len(chromosome_sequence)):
                raise ValueError(f"Start of segment {start_of_segment} is out of range of the chromosome length "
                                 f"{len(self.get_chromosome_sequence(chromosome_name))}")
        if segment_length is None:
            segment_length = len(chromosome_sequence)
        else:
            if (segment_length <= 0) or (start_of_segment + segment_length > len(chromosome_sequence)):
                raise ValueError(f"Segment length {segment_length} is out of range of the chromosome length "
                                 f"{len(self.get_chromosome_sequence(chromosome_name))}")
        segment_sequence = chromosome_sequence[start_of_segment:start_of_segment + segment_length]

        fcgr = CGR(segment_sequence, k_mer).get_fcgr()
        kmer_matrix = CGR(segment_sequence, k_mer).get_kmer_matrix()

        # flat_indices = np.argsort(fcgr, axis=None)[-top_n:][::-1]
        # top_kmers = []
        # for idx in flat_indices:
        #     y, x = np.unravel_index(idx, fcgr.shape)
        #     kmer = kmer_matrix[y, x]
        #     count = fcgr[y, x]
        #     top_kmers.append((kmer, count))
        flat_fcgr = fcgr.flatten()

        zero_count = np.sum(flat_fcgr == 0)
        if zero_count > 0:
            print(f"The FCGR contains {zero_count} k-mers with zero frequency.")

        # Only consider non-zero entries
        non_zero_indices = np.nonzero(flat_fcgr)[0]
        sorted_indices = non_zero_indices[np.argsort(flat_fcgr[non_zero_indices])]

        bottom_indices = sorted_indices[:bottom_n]
        top_indices = sorted_indices[-top_n:][::-1]

        def idx_to_info(idx):
            y, x = np.unravel_index(idx, fcgr.shape)
            return kmer_matrix[y, x], fcgr[y, x]

        bottom_kmers = [idx_to_info(idx) for idx in bottom_indices]
        top_kmers = [idx_to_info(idx) for idx in top_indices]

        return fcgr, top_kmers, bottom_kmers

    def plot_fcgr(self, chromosome_name, start_of_segment=None, segment_length=None, k_mer=9, fcgr_cgr='fcgr',
                  label=True, global_min=None, global_max=None, _3d=False, resolution=None, bits=None):
        chromosome_sequence = self.get_chromosome_sequence(chromosome_name)
        if start_of_segment is None:
            start_of_segment = 0
        else:
            if (start_of_segment < 0) or (start_of_segment >= len(chromosome_sequence)):
                raise ValueError(f"Start of segment {start_of_segment} is out of range of the chromosome length "
                                 f"{len(self.get_chromosome_sequence(chromosome_name))}")
        if segment_length is None:
            segment_length = len(chromosome_sequence)
        else:
            if (segment_length <= 0) or (start_of_segment + segment_length > len(chromosome_sequence)):
                raise ValueError(f"Segment length {segment_length} is out of range of the chromosome length "
                                 f"{len(self.get_chromosome_sequence(chromosome_name))}")
        segment_sequence = chromosome_sequence[start_of_segment:start_of_segment + segment_length]

        if resolution is None:
            resolution = RESOLUTION_DICT[k_mer]
        if bits is None:
            bits = BITS_DICT[self.species]

        if fcgr_cgr == 'fcgr':
            fcgr = CGR(segment_sequence, k_mer).get_fcgr()
            rescale_fcgr, _ = CGR.array2img(fcgr, bits=bits, resolution=resolution, m=global_min, M=global_max,
                                            return_array=True)
            rescale_2, fcgr_image = CGR.array2img(rescale_fcgr, bits=bits, resolution=resolution, return_array=True)
            fcgr_pil = Image.fromarray(fcgr_image, 'L')
        else:
            fcgr = CGR(segment_sequence, k_mer).get_cgr()
            fcgr_image = CGR.array2img(fcgr, bits=bits, resolution=resolution)
            fcgr_pil = Image.fromarray(fcgr_image, 'L')

        if _3d:
            fig = go.Figure()
            fig.add_trace(go.Surface(
                z=fcgr,  # FCGR values as the height map
                colorscale="Blues",  # Choose a color scheme
                showscale=True  # Show color bar
            ))
            # Set layout for better visualization
            fig.update_layout(
                title=f"Interactive 3D FCGR Plot for chromosome {chromosome_name}",
                scene=dict(
                    xaxis_title="X Axis",
                    yaxis_title="Y Axis",
                    zaxis_title="FCGR Value",
                )
            )
            # Show the interactive plot
            fig.show()

        plt.imshow(fcgr_pil, cmap="gray")
        plt.xticks([])  # Remove x ticks
        plt.yticks([])  # Remove y ticks

        # Coordinates of points to label (example points)
        if k_mer == 9:
            points = [(-12, 520), (-15, -6), (522, -6), (525, 520)]  # List of (x, y) tuples
        elif k_mer == 6:
            points = [(-2, 65), (-2, -1), (65, -1), (65, 65)]

        labels = ["A", "C", "G", "T"]  # Labels for each point

        if label:
            if k_mer == 6 or k_mer == 9:
                for (x, y), label in zip(points, labels):
                    plt.text(x, y, label, color='black', fontsize=14, ha='center', va='center')
            plt.title(f"Chromosome {chromosome_name}", fontsize=14)

        # Get the absolute path of the project root directory
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        save_path = os.path.join(project_root, 'Figures', 'FCGRs', self.species)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        plt.savefig(
            f"{save_path}/chr-{chromosome_name}_range-{start_of_segment}-{start_of_segment + segment_length}_"
            f"kmer-{k_mer}_{fcgr_cgr}.png",
            dpi=300, bbox_inches='tight', transparent=True)
        # plt.show()
        plt.close()

    def find_n_counts(self):
        all_len = 0
        sequence_removed_n_len = 0
        for chromosome in self.get_all_chromosomes_name():
            sequence = self.get_chromosome_sequence(chromosome)
            all_len += len(sequence)
            sequence_removed_n = sequence.replace("N", "")
            sequence_removed_n_len += len(sequence_removed_n)
            if len(sequence) != len(sequence_removed_n):
                print(chromosome, len(sequence), len(sequence_removed_n))
        print(all_len, sequence_removed_n_len)
        print(all_len - sequence_removed_n_len)
        print((1 - (sequence_removed_n_len / all_len)) * 100)

    def clear_cache(self):
        self._chromosome_sequence_cache = {}

    def create_chromosomes_files(self, whole_genome_path):
        """
            Run this method to create chromosome files from a whole genome file (Run only once)
        """
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
        for chromosome_name in self.get_all_chromosomes_name(include_whole_genome=True):
            self.cytobands[chromosome_name] = {}
        if self.species == "Human":
            file_path = os.path.join(self.root_path, self.species, 'bedfiles', 'telomere.bed')
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

            file_path = os.path.join(self.root_path, self.species, 'bedfiles', 'centromere.bed')
            chr_name = None
            last_start = None
            last_end = None
            with open(file_path) as file:
                for index, line in enumerate(file):
                    if index == 0:
                        continue
                    line = line.strip()
                    parts = line.split("\t")

                    if chr_name != parts[0]:
                        if last_end is not None:
                            self.cytobands[chr_name]["centromere"] = AnnotationRecord(chr_name=chr_name,
                                                                                      start=last_start,
                                                                                      end=last_end,
                                                                                      name="cent",
                                                                                      group="centromere",
                                                                                      color=None)

                        last_start = int(parts[1])
                        chr_name = parts[0][3:]

                    last_end = int(parts[2])

            file_path = os.path.join(self.root_path, self.species, 'bedfiles', 'cytobands.bed')
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
        match = re.search(r'(chromosome\s*|chr|scaffold_|contig_)(\d+|[A-Za-z]+)', string, re.IGNORECASE)
        if match:
            return match.group(2)
        if match is None:
            return "N"

    @staticmethod
    def get_reverse_complement(sequence):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        bases = [complement[base] for base in sequence]
        bases = reversed(bases)
        return ''.join(bases)

    @staticmethod
    def read_fasta(file_path):
        sequence = ""
        with open(file_path) as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    file_name = line.split(",")[0].split(" ")[-1]
                else:
                    sequence += line
        return file_name, sequence


if __name__ == '__main__':
    # specie = "Human"
    # # ChromosomesHolder(specie).create_chromosomes_files("GCA_000011425.1_ASM1142v1_genomic.fna")
    # genome = ChromosomesHolder(specie)
    # genome.plot_fcgr("1", start_of_segment=None, segment_length=None, k_mer=9, fcgr_cgr='fcgr',
    #                  label=True, global_min=None, global_max=None, _3d=True)
    # genome.plot_fcgr("2", start_of_segment=None, segment_length=None, k_mer=9, fcgr_cgr='fcgr',
    #                  label=True, global_min=None, global_max=None, _3d=True)

    genome = ChromosomesHolder("Maize")
    for chr_name in genome.get_all_chromosomes_name():
        print(f"chromosome length is {len(genome.get_chromosome_sequence(chr_name))}")
        fcgr, top_kmers, _ = genome.get_max_kmer_values(chr_name, start_of_segment=None, segment_length=None,
                                                        k_mer=9, top_n=3, bottom_n=0)
        print(f"Top 3 most frequent k-mers in chromosome {chr_name} of Maize:")
        for kmer, count in top_kmers:
            print(f"{kmer}: {count}")
            print(f"{count / np.sum(fcgr) * 100}")

        # print(f"Last 10 most frequent k-mers in chromosome {chr_name} of Human:")
        # for kmer, count in bottom_kmers:
        #     print(f"{kmer}: {count}")

    # genome.find_n_counts()
    # genome.plot_fcgr("21", k_mer=9, fcgr_cgr='fcgr', label=True)

    '''FCGR global min and max test'''
    # global_m = np.inf
    # global_M = 0
    # for chr_name in tqdm(genome.get_all_chromosomes_name()):
    #     seq = genome.get_chromosome_sequence(chr_name)
    #     fcgr = CGR(seq, k_mer=9).get_fcgr()
    #     m, M = fcgr.min(), fcgr.max()
    #     global_m = min(global_m, m)
    #     global_M = max(global_M, M)
    # print(global_m, global_M)  # 0 1134132 human / 21 246147 maize

    # global_m = 21
    # global_M = 246147
    # for chr_name in tqdm(genome.get_all_chromosomes_name()):
    # genome.plot_fcgr("21", k_mer=9, fcgr_cgr='fcgr', label=True, global_min=global_m, global_max=global_M)

    # '''Genome length test'''
    # l_so_far = 0
    # lengths = []
    # for chr_name in tqdm(genome.get_all_chromosomes_name()):
    #     seq = genome.get_chromosome_sequence(chr_name)
    #     l = len(seq)
    #     lengths.append(l)
    #     print(f"chr {chr_name} length = ", l)
    #     l_so_far += l
    #     # print("Length so far is: ", l_so_far)
    # print("Total length: ", sum(lengths))
    # print("Max length: ", max(lengths))
    # print("Min length: ", min(lengths))
    # print("Mean length: ", np.mean(lengths))
