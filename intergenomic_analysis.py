import ast
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from tqdm import tqdm

from chaos_game_representation import CGR
from chromosomes_holder import ChromosomesHolder
from constants import DISTANCE_METRICS_LIST, SCIENTIFIC_NAMES
from distances.distance_metrics import get_dist


class InterGenomicAnalysis:
    def __init__(self, base_specie, target_species_list, length=500_000, num_samples=100, run=True):
        self.base_specie = base_specie
        self.target_species_list = target_species_list
        self.num_samples = num_samples

        if run:
            self.base_segments = []
            chromosomes_holder_base = ChromosomesHolder(self.base_specie)
            for _ in tqdm(range(self.num_samples)):
                if self.base_specie == "Human":
                    remove_outlier = True
                else:
                    remove_outlier = False
                random_segment = chromosomes_holder_base.get_random_segment(length, remove_outlier=remove_outlier,
                                                                            return_dict=False)
                self.base_segments.append(random_segment)

            self.species_segments_dict = {}
            for target_species in self.target_species_list:
                self.species_segments_dict[target_species] = []
                chromosomes_holder_target = ChromosomesHolder(target_species)
                for _ in tqdm(range(100)):
                    if target_species == "Human":
                        remove_outlier = True
                    else:
                        remove_outlier = False
                    random_segment = chromosomes_holder_target.get_random_segment(length, remove_outlier=remove_outlier,
                                                                                  return_dict=False)
                    self.species_segments_dict[target_species].append(random_segment)

    def run_experiment(self, trim=True, new_run=False):
        experiment_path = os.path.join('outputs', 'intergenome_analysis.csv')
        if not os.path.exists(experiment_path):
            new_run = True
        if new_run:
            df = pd.DataFrame(columns=['Category', 'Distance', 'Mean', 'Variance', 'Max', 'Distance_Values'])
            for target_species in self.target_species_list:
                print(f"Running experiment for {target_species}...")

                for distance_metric in DISTANCE_METRICS_LIST:
                    distances = []
                    for i in tqdm(range(self.num_samples)):
                        base_fcgr = CGR(self.base_segments[i], k_mer=6).get_fcgr()
                        target_fcgr = CGR(self.species_segments_dict[target_species][i], k_mer=6).get_fcgr()
                        distances.append(get_dist(base_fcgr, target_fcgr, distance_metric))

                    distances = np.asarray(distances)
                    if trim:
                        percent_to_trim = 5  # Total percentage to trim (2.5% from each side)
                        distances = self.trim_outliers(distances, percent_to_trim)

                    df.loc[len(df)] = [target_species, distance_metric, np.mean(distances), np.var(distances),
                                       np.max(distances), distances]
                    df.to_csv(experiment_path, index=False, sep='\t')

        df = pd.read_csv(experiment_path, sep='\t')
        return df

    @staticmethod
    def plot_means_variances(df, plot_type='barplot'):
        save_path = os.path.join('Figures', 'intergenomic_experiment')
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        # Iterate over each distance metric
        for distance_m in df['Distance'].unique():
            subset_df = df[df['Distance'] == distance_m]
            categories = subset_df['Category'].replace(SCIENTIFIC_NAMES)

            plt.figure(figsize=(8, 6))

            # Now you can create the box plot
            if plot_type == 'boxplot':
                subset_df.loc[:, 'Distance_Values'] = subset_df['Distance_Values'].apply(ast.literal_eval)
                exploded_df = subset_df.explode('Distance_Values')
                exploded_df.loc[:, 'Distance_Values'] = exploded_df['Distance_Values'].astype(float)

                distance_data_dict = {}
                for cat in categories:
                    group = exploded_df[exploded_df['Category'].replace(SCIENTIFIC_NAMES) == cat]
                    distance_data_dict[cat] = list(group['Distance_Values'])

                # Prepare the data for the boxplot
                distance_data = [np.array(distance_data_dict[cat]) for cat in categories if cat in distance_data_dict]

                plt.boxplot(distance_data)

                # Set x-tick labels as categories, preserving the order
                plt.xticks(ticks=np.arange(1, len(distance_data) + 1), labels=categories)

            elif plot_type == 'barplot':
                # Create the bar plot based on the mean and variance
                means = subset_df['Mean']
                stds = np.sqrt(subset_df['Variance'])
                bars = plt.bar(categories, means, yerr=stds, capsize=5, color='skyblue', ecolor='red',
                               edgecolor='black')

                # Annotate the bars with the value
                for i, bar in enumerate(bars):
                    plt.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height(),  # + stds.iloc[i] + 0.1,
                        f'{means.iloc[i]:.2f}',
                        ha='center',
                        va='bottom',
                        color='black',  # Change the text color
                        fontsize=10,
                        bbox=dict(facecolor='white', edgecolor='none', pad=0.05, alpha=0.6)
                    )

                # Set the y-axis label
                plt.ylabel(distance_m)

            # Rotate category labels for better readability
            plt.xticks(rotation=15, ha='right', fontstyle='italic')

            # Save the plot
            plt.savefig(f"{save_path}/{distance_m}_{plot_type}.png", bbox_inches='tight', transparent=False)
            plt.close()

    @staticmethod
    def trim_outliers(data, percent):
        # Calculate the lower and upper percentiles
        lower_percentile = percent / 2
        upper_percentile = 100 - lower_percentile

        # Calculate the lower and upper threshold values
        lower_threshold = np.percentile(data, lower_percentile)
        upper_threshold = np.percentile(data, upper_percentile)

        # Filter out the outliers
        trimmed_data = [x for x in data if lower_threshold <= x <= upper_threshold]

        return trimmed_data


if __name__ == '__main__':
    target_list = ["Human", "Chimp", "Mouse", "Insect", "Fungus", "Plant", "Protist", "Archaea", "Bacteria"]
    intergenome = InterGenomicAnalysis(base_specie="Human", target_species_list=target_list, run=False)
    data_frame = intergenome.run_experiment(trim=True, new_run=False)
    intergenome.plot_means_variances(data_frame, plot_type='boxplot')
