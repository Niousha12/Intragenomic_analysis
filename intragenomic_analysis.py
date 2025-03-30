import os
import random
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from tqdm import tqdm

from chromosomes_holder import ChromosomesHolder
from chaos_game_representation import CGR
from constants import DISTANCE_METRICS_LIST, HETERO_HETERO_DICT, HETERO_EU_DICT, TANDEM_REPEAT_DICT
from distances.distance_metrics import get_dist

random.seed(46)
np.random.seed(46)


class IntraGenomicAnalysis:
    def __init__(self, specie, kmer):
        self.chromosomes_holder = ChromosomesHolder(specie)
        self.kmer = kmer

    def telomere_vs_telomere(self, distance_metrics: list):
        first_chromosome_name = self.chromosomes_holder.get_all_chromosomes_name()[0]
        reference_telomere = self.chromosomes_holder.get_cytoband_segment(first_chromosome_name, "tel_1")
        reference_cgr = CGR(reference_telomere, self.kmer).get_fcgr()

        distance_values = {}
        for distance_metric in distance_metrics:
            distance_values[distance_metric] = []

        for chromosome_name in tqdm(self.chromosomes_holder.get_all_chromosomes_name()):
            if chromosome_name == first_chromosome_name:
                continue
            target_telomere = self.chromosomes_holder.get_cytoband_segment(chromosome_name, "tel_1")
            target_cgr = CGR(target_telomere, self.kmer).get_fcgr()
            for distance_metric in distance_metrics:
                distance = get_dist(reference_cgr, target_cgr, distance_metric)
                distance_values[distance_metric].append(distance)

        for distance_metric in distance_metrics:
            distance_values[distance_metric] = np.mean(np.asarray(distance_values[distance_metric]))
        return distance_values

    def cytoband_vs_cytoband(self, distance_metrics, cytoband_dictionary, exclude_chromosome="Y"):
        distance_values = {}
        for distance_metric in distance_metrics:
            distance_values[distance_metric] = []

        for chromosome_name, chromosome_cytobands in tqdm(cytoband_dictionary.items()):
            if chromosome_name == exclude_chromosome:
                continue
            reference_cytoband = self.chromosomes_holder.get_cytoband_segment(chromosome_name, chromosome_cytobands[0])
            reference_cgr = CGR(reference_cytoband, self.kmer).get_fcgr()

            for target_cytoband in chromosome_cytobands[1:]:
                target_cytoband = self.chromosomes_holder.get_cytoband_segment(chromosome_name, target_cytoband)
                target_cgr = CGR(target_cytoband, self.kmer).get_fcgr()
                for distance_metric in distance_metrics:
                    distance = get_dist(reference_cgr, target_cgr, distance_metric)
                    distance_values[distance_metric].append(distance)

        for distance_metric in distance_metrics:
            distance_values[distance_metric] = np.mean(np.asarray(distance_values[distance_metric]))
        return distance_values

    def tandem_repeat_vs_tandem_repeat(self, distance_metrics, cytoband_dictionary, exclude_chromosome=None):
        if exclude_chromosome is not None and exclude_chromosome in cytoband_dictionary:
            del cytoband_dictionary[exclude_chromosome]
        first_chromosome_name = list(cytoband_dictionary.keys())[0]
        reference_sequence = self.chromosomes_holder.get_cytoband_segment(first_chromosome_name,
                                                                          cytoband_dictionary[first_chromosome_name][0])
        reference_cgr = CGR(reference_sequence, self.kmer).get_fcgr()

        distance_values = {}
        for distance_metric in distance_metrics:
            distance_values[distance_metric] = []

        for chromosome_name, chromosome_cytobands in tqdm(cytoband_dictionary.items()):
            if chromosome_name == first_chromosome_name:
                continue
            for target_cytoband in chromosome_cytobands:
                target_cytoband = self.chromosomes_holder.get_cytoband_segment(chromosome_name, target_cytoband)
                target_cgr = CGR(target_cytoband, self.kmer).get_fcgr()
                for distance_metric in distance_metrics:
                    distance = get_dist(reference_cgr, target_cgr, distance_metric)
                    distance_values[distance_metric].append(distance)

        for distance_metric in distance_metrics:
            distance_values[distance_metric] = np.mean(np.asarray(distance_values[distance_metric]))
        return distance_values

    def p_vs_q(self, chromosome_name, distance_metrics, end_p, sequence_length=500000, num_samples=100):
        chromosome_sequence = self.chromosomes_holder.get_chromosome_sequence(chromosome_name)

        distance_values = {}
        for distance_metric in distance_metrics:
            distance_values[distance_metric] = []

        for _ in tqdm(range(num_samples)):
            random_start_p = random.randint(0, end_p - sequence_length)
            random_start_q = random.randint(end_p, len(chromosome_sequence) - sequence_length)

            cgr_p = CGR(chromosome_sequence[random_start_p:random_start_p + sequence_length], self.kmer).get_fcgr()
            cgr_q = CGR(chromosome_sequence[random_start_q:random_start_q + sequence_length], self.kmer).get_fcgr()

            for distance_metric in distance_metrics:
                distance = get_dist(cgr_p, cgr_q, distance_metric)
                distance_values[distance_metric].append(distance)

        for distance_metric in distance_metrics:
            distance_values[distance_metric] = np.mean(np.asarray(distance_values[distance_metric]))
        return distance_values

    def arbitrary_vs_arbitrary(self, distance_metrics, sequence_length=500000, num_samples=100):
        distance_values = {}
        for distance_metric in distance_metrics:
            distance_values[distance_metric] = []

        for _ in tqdm(range(num_samples)):
            segments_info = self.chromosomes_holder.get_random_segments_list(sequence_length, num_segments=2,
                                                                             overlap=False)
            segment_1_sequence = segments_info['segments_sequences'][0]
            segment_2_sequence = segments_info['segments_sequences'][1]

            segment_1_cgr = CGR(segment_1_sequence, self.kmer).get_fcgr()
            segment_2_cgr = CGR(segment_2_sequence, self.kmer).get_fcgr()
            for distance_metric in distance_metrics:
                distance = get_dist(segment_1_cgr, segment_2_cgr, distance_metric)
                distance_values[distance_metric].append(distance)

        for distance_metric in distance_metrics:
            distance_values[distance_metric] = np.mean(np.asarray(distance_values[distance_metric]))
        return distance_values

    def run_experiment(self, new_run=False):
        experiment_path = os.path.join('outputs', 'intragenome_analysis.csv')
        if not os.path.exists(experiment_path):
            new_run = True
        if new_run:
            experiments_list = ['Telomere vs. Telomere', 'Heterochromatin vs. Heterochromatin',
                                'Heterochromatin vs. Euchromatin', 'Y (p-arm vs. q-arm)', '13 (p-arm vs. q-arm)',
                                '14 (p-arm vs. q-arm)', '15 (p-arm vs. q-arm)', '21 (p-arm vs. q-arm)',
                                '22 (p-arm vs. q-arm)', 'Large Tandem Repeat Arrays-P1',
                                'Large Tandem Repeat Arrays-P2', 'Arbitrary Sequences']

            df = pd.DataFrame({'Experiment': experiments_list})

            data = [
                self.telomere_vs_telomere(DISTANCE_METRICS_LIST),
                self.cytoband_vs_cytoband(DISTANCE_METRICS_LIST, HETERO_HETERO_DICT),
                self.cytoband_vs_cytoband(DISTANCE_METRICS_LIST, HETERO_EU_DICT, exclude_chromosome="Y"),
                self.p_vs_q("Y", DISTANCE_METRICS_LIST, 10724418),
                self.p_vs_q("13", DISTANCE_METRICS_LIST, 16522942),
                self.p_vs_q("14", DISTANCE_METRICS_LIST, 11400261),
                self.p_vs_q("15", DISTANCE_METRICS_LIST, 17186630),
                self.p_vs_q("21", DISTANCE_METRICS_LIST, 11134529),
                self.p_vs_q("22", DISTANCE_METRICS_LIST, 14249622),
                self.tandem_repeat_vs_tandem_repeat(DISTANCE_METRICS_LIST, TANDEM_REPEAT_DICT),
                self.tandem_repeat_vs_tandem_repeat(DISTANCE_METRICS_LIST, TANDEM_REPEAT_DICT, exclude_chromosome="Y"),
                self.arbitrary_vs_arbitrary(DISTANCE_METRICS_LIST)
            ]
            data_df = pd.DataFrame(data, columns=DISTANCE_METRICS_LIST)
            df = pd.concat([df, data_df], axis=1)
            df.to_csv(experiment_path, index=False, sep='\t')

        df = pd.read_csv(experiment_path, sep='\t', header=0, index_col=0)
        return df

    @staticmethod
    def plot_intragenomic_analysis(df):
        save_path = os.path.join('Figures', 'intragenomic_experiment')
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        # colors = plt.cm.tab20.colors
        # colors = ['#C53A33',
        #           '#F19D99',
        #           '#F5BE82',
        #
        #           '#3A75AE',
        #           '#B3C6E5',
        #           '#519E3E',
        #           '#A7DD93',
        #           '#8D69B8',
        #           '#C1B1D2',
        #
        #           '#EF8636',
        #           '#BE9E96',
        #           '#83584D']

        colors = ['#ADD8E6',
                  '#6495ED',
                  '#1E3A8A',

                  '#FF4500',
                  '#FF6347',
                  '#FF7F50',
                  '#FF8C00',
                  '#FFA07A',
                  '#FFD700',

                  '#9370DB',
                  '#6A0DAD',
                  '#A89F91']

        exclude_color = '#6A0DAD'  # Color to remove
        exclude_index = colors.index(exclude_color) if exclude_color in colors else None

        # Remove corresponding row in df
        if exclude_index is not None and exclude_index < len(df):
            df = df.drop(df.index[exclude_index])  # Drop experiment by index
            colors.pop(exclude_index)  # Remove color to keep alignment

        metrics = list(df.columns)
        experiments = list(df.index)
        data = df.transpose().values  # Shape: (num_metrics, num_experiments)
        bar_width = 0.8
        x = np.arange(len(experiments))

        fig, axes = plt.subplots(1, len(data), figsize=(16, 5))

        for i, ax in enumerate(axes):
            for j in range(len(experiments)):
                if experiments[j] == "Large Tandem Repeat Arrays-P1":
                    experiments[j] = "Large Tandem Repeat Arrays"
                ax.bar(j, data[i][j], width=bar_width, color=colors[j], label=experiments[j] if i == 0 else "")
            ax.set_title(metrics[i])
            ax.grid(True, axis='y', color='lightgrey', linestyle='-', linewidth=0.5)  # Adjust grid properties
            ax.set_xticks(x)
            ax.set_xticklabels([])
            ax.tick_params(axis='y', pad=1)  # labelsize=8,

        handles, labels = axes[0].get_legend_handles_labels()

        fig.legend(handles, labels, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 0.1),
                   prop=FontProperties(style='italic', size=12))

        plt.tight_layout(rect=[0, 0.1, 1, 1])
        plt.subplots_adjust(wspace=0.4)

        plt.savefig(f"{save_path}/intragenomic_analysis.png", dpi=300, bbox_inches='tight', transparent=True)
        # plt.show()


if __name__ == '__main__':
    intragenome = IntraGenomicAnalysis('Human', kmer=6)
    dataframe = intragenome.run_experiment(new_run=False)
    intragenome.plot_intragenomic_analysis(dataframe)
