import math
import os
import pickle

import numpy as np
import pandas as pd
from sklearn.manifold import MDS
from matplotlib import pyplot as plt
from tqdm import tqdm

from chaos_game_representation import CGR
from chromosomes_holder import ChromosomesHolder, AnnotationRecord
from distances.distance_metrics import get_dist
import plotly.express as px


class ChromosomeRepresentativeSelection:
    def __init__(self, specie, kmer, distance_metric, segment_length=None):
        self.specie = specie
        self.chromosomes_holder = ChromosomesHolder(specie)
        self.kmer = kmer
        if segment_length is not None:
            self.length = segment_length
        else:
            self.length = self.chromosomes_holder.get_appropriate_segment_length(scale='genome')
        print(f"Segment length for {self.specie} is {self.length}")
        self.distance_metric = distance_metric

        self.pickle_path_root = f"./outputs/cache_pickles/{self.specie}/"
        if not os.path.exists(self.pickle_path_root):
            os.makedirs(self.pickle_path_root)

    def get_representative(self, chromosome_name):
        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name)
        distance_matrix = self.get_distance_matrix(chromosome_name)

        centroid = self.find_centroid(distance_matrix, exclude_indices=None)

        return {"sequence": segments['segments_sequences'][centroid],
                "index": centroid,
                "fcgr": fcgrs[centroid],
                "type": "non-approximative"}

    def get_random_representative(self, chromosome_name, outlier=True):
        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name)

        random_centroid_index = self.chromosomes_holder.choose_random_fragment_index(chromosome_name, self.length,
                                                                                     outlier=outlier)

        return {"sequence": segments['segments_sequences'][random_centroid_index],
                "index": random_centroid_index,
                "fcgr": fcgrs[random_centroid_index],
                "type": f"random_outlier_{outlier}"}

    def get_approximate_representative(self, chromosome_name, random_sequences_number=30,
                                       remove_outliers_function="IQR", verbose=False):
        if verbose:
            print(f"Finding the approximate representative for chromosome {chromosome_name}")
        random_sequences_list = []
        count, outlier_indices_number_total = 0, 0
        avgs = None
        while len(random_sequences_list) < random_sequences_number:
            count += 1
            for _ in range(random_sequences_number - len(random_sequences_list)):
                random_sequence_dict = self.chromosomes_holder.get_random_segment(self.length, chromosome_name,
                                                                                  return_dict=True)
                fcgr = CGR(random_sequence_dict['sequence'], self.kmer).get_fcgr()
                random_sequence_dict['fcgr'] = fcgr

                random_sequences_list.append(random_sequence_dict)

            # Get the distance matrix between these choices, then get the average distance for each choice
            distance_matrix = np.zeros((len(random_sequences_list), len(random_sequences_list)))
            for i in tqdm(range(distance_matrix.shape[0])):
                for j in range(i + 1):
                    distance_matrix[i, j] = get_dist(random_sequences_list[i]['fcgr'], random_sequences_list[j]['fcgr'],
                                                     dist_m=self.distance_metric)
                    distance_matrix[j, i] = distance_matrix[i, j]
            avgs = np.mean(distance_matrix, axis=1)
            if verbose:
                print(f"Averages at count {count} are : {avgs}, "
                      f"Mean of these averages : {np.mean(avgs)}")

            # Find the outlier indices
            if remove_outliers_function == "ZSCORE":
                outlier_indices = self.get_outliers_index_zscore(avgs)
            elif remove_outliers_function == "IQR":
                outlier_indices = self.get_outliers_index_iqr(avgs)
            else:
                raise ValueError("Invalid outlier removal function")
            outlier_indices_number_total += len(outlier_indices)
            if verbose:
                print(f"indices dropped at count {count} are {outlier_indices}")

            # Remove the outliers from the list of dictionaries
            random_sequences_list = [item for idx, item in enumerate(random_sequences_list) if
                                     idx not in outlier_indices]

        # Choose the minimum average distance from the remaining as the representative
        representative_dict = random_sequences_list[np.argmin(avgs)]
        if verbose:
            print(f"The process ran for {count} times, "
                  f"Total number of dropped indices {outlier_indices_number_total}")

        return {"sequence": representative_dict['sequence'],
                "start": representative_dict['start'],
                "fcgr": representative_dict['fcgr'],
                "type": "approximative"}

    def get_distance_from_representative(self, chromosome_name, representative_dict):
        if representative_dict['type'] == "approximative":
            fcgrs = self.get_fcgrs_of_segments(chromosome_name)
            distance_from_representative = np.zeros(len(fcgrs))

            for index, fcgr in enumerate(fcgrs):
                distance_from_representative[index] = get_dist(fcgr, representative_dict['fcgr'],
                                                               dist_m=self.distance_metric)
        else:
            distance_matrix = self.get_distance_matrix(chromosome_name)
            distance_from_representative = distance_matrix[representative_dict['index'], :]

        return distance_from_representative

    def get_approximation_error(self, chromosome_name, random_outliers=True):
        representative = self.get_representative(chromosome_name)
        distances_from_representative = self.get_distance_from_representative(chromosome_name, representative)

        approximative_representative = self.get_approximate_representative(chromosome_name)
        distances_from_approximative_representative = \
            self.get_distance_from_representative(chromosome_name, approximative_representative)

        approximation_error = np.mean(
            np.abs((distances_from_approximative_representative - distances_from_representative)))
        # np.sqrt(np.mean((distances_from_approximative_representative - distances_from_representative) ** 2))

        print(f"Approximation error: {approximation_error}")

        random_error = None
        if random_outliers:
            random_representative = self.get_random_representative(chromosome_name)
            distances_from_random_representative = \
                self.get_distance_from_representative(chromosome_name, random_representative)
            random_error = np.mean(
                np.abs((distances_from_random_representative - distances_from_representative)))
            # np.sqrt(np.mean((distances_from_random_representative - distances_from_representative) ** 2))
            print(f"Random error: {random_error}")

        return distances_from_representative, approximation_error, random_error

    def plot_distance_variations(self, chromosome_name, plot_random_outliers=True, plot_approximate=True, x_range=None):
        figure_path = os.path.join('Figures', 'Representative', self.specie, 'Different_centroids')
        if not os.path.exists(figure_path):
            os.makedirs(figure_path)

        plt.figure(figsize=(10, 5))

        pipeline_representative_dict = self.get_representative(chromosome_name)
        distances_from_centroid = self.get_distance_from_representative(chromosome_name, pipeline_representative_dict)

        if plot_random_outliers:
            random_outlier_representative_dict = self.get_random_representative(chromosome_name, True)
            distance_from_outlier_representative = \
                self.get_distance_from_representative(chromosome_name, random_outlier_representative_dict)
            plt.plot(distance_from_outlier_representative, marker='o', linestyle='-', markersize=4, color='purple')

        if plot_approximate:
            approximate_representative_dict = self.get_approximate_representative(chromosome_name, verbose=True)
            distance_from_approximate_representative = \
                self.get_distance_from_representative(chromosome_name, approximate_representative_dict)

            plt.plot(distance_from_approximate_representative, marker='o', linestyle='-', markersize=4, color='blue')

        plt.plot(distances_from_centroid, marker='o', linestyle='-', markersize=4, color='red')
        plt.grid(True)

        if chromosome_name == "Whole Genome":
            include_whole_genome = True
        else:
            include_whole_genome = False

        if x_range is None:
            x_range = math.ceil(
                self.chromosomes_holder.get_largest_chromosome_length(include_whole_genome) / self.length / 20)
        else:
            x_range = x_range
        x_ticks = []
        for i in range(0, int(x_range) + 1):
            x_ticks.append(i * 20)
        y_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        plt.xticks(x_ticks)
        plt.yticks(y_ticks)
        plt.savefig(f'{figure_path}/distance_from_centroid_{chromosome_name}_k{self.kmer}_size{self.length}.png')

    def get_non_overlapping_segments(self, chromosome_name):
        pickle_path = os.path.join(self.pickle_path_root, f"chr_{chromosome_name}_len_{self.length}_segments.pickle")
        if os.path.exists(pickle_path):
            return self.load_pickle(pickle_path)
        else:
            segments = self.chromosomes_holder.get_chromosome_non_overlapping_segments(chromosome_name, self.length)
            self.create_pickle(segments, pickle_path)
            return segments

    def get_fcgrs_of_segments(self, chromosome_name):
        pickle_path = os.path.join(self.pickle_path_root,
                                   f"chr_{chromosome_name}_len_{self.length}_kmer_{self.kmer}_fcgrs.pickle")
        if os.path.exists(pickle_path):
            return self.load_pickle(pickle_path)
        else:
            segments = self.get_non_overlapping_segments(chromosome_name)

            fcgrs_list = []
            for segment in segments['segments_sequences']:
                fcgrs_list.append(CGR(segment, self.kmer).get_fcgr())
            self.create_pickle(fcgrs_list, pickle_path)
            return fcgrs_list

    def get_distance_matrix(self, chromosome_name):
        pickle_path = os.path.join(self.pickle_path_root,
                                   f"chr_{chromosome_name}_len_{self.length}_kmer_{self.kmer}"
                                   f"_dist_{self.distance_metric}_distance_matrix.pickle")
        if os.path.exists(pickle_path):
            return self.load_pickle(pickle_path)
        else:
            fcgrs = self.get_fcgrs_of_segments(chromosome_name)

            distance_matrix = np.zeros((len(fcgrs), len(fcgrs)))
            for i in tqdm(range(distance_matrix.shape[0])):
                for j in range(i + 1):
                    distance_matrix[i, j] = get_dist(fcgrs[i], fcgrs[j], dist_m=self.distance_metric)
                    distance_matrix[j, i] = distance_matrix[i, j]
            self.create_pickle(distance_matrix, pickle_path)
            return distance_matrix

    @staticmethod
    def create_pickle(results, pickle_path):
        if not os.path.exists(pickle_path):
            with open(pickle_path, 'wb') as handle:
                pickle.dump(results, handle)
        else:
            print("Pickle file already exists")

    @staticmethod
    def load_pickle(pickle_path):
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as handle:
                results = pickle.load(handle)
            return results
        else:
            raise FileNotFoundError("Pickle file not found")

    @staticmethod
    def find_centroid(distance_matrix, median=False, exclude_indices=None):
        if exclude_indices is None:
            exclude_indices = []

        mask = np.ones(distance_matrix.shape[0], dtype=bool)
        mask[exclude_indices] = False
        # Apply the mask to rows and columns
        masked_matrix = distance_matrix[mask][:, mask]

        # Sum distances from each remaining point to all others
        sums = np.mean(masked_matrix, axis=1)

        # The index of the minimum sum corresponds to the centroid among the included indices, or is the median
        if median:
            centroid_index = np.argsort(sums)[len(sums) // 2]
        else:
            centroid_index = np.argmin(sums)

        # Convert the index back to the original matrix indices
        original_indices = np.arange(distance_matrix.shape[0])[mask]
        centroid_original_index = original_indices[centroid_index]

        return centroid_original_index

    @staticmethod
    def get_outliers_index_zscore(data, threshold=0.1):
        mean = np.mean(data)
        std = np.std(data)

        outlier_indices = []
        for index, x in enumerate(data):
            if not x < (mean + threshold):
                outlier_indices.append(index)

        return outlier_indices

    @staticmethod
    def get_outliers_index_iqr(data, multiplier=1.5):
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        upper_bound = q3 + multiplier * iqr

        outlier_indices = []
        for index, x in enumerate(data):
            if not x <= upper_bound:
                outlier_indices.append(index)
        return outlier_indices

    @staticmethod
    def plot_multi_dimensional_scaling(distance_matrix, info: list, representative_index=None, coloring_type='NCBI'):
        mds = MDS(n_components=3, dissimilarity='precomputed', random_state=42)
        mds_result = mds.fit_transform(distance_matrix)

        # Add description for each data
        data_description = np.asarray(
            [f"<br>Start: {info[i].start}<br>band: {info[i].name}" for i in range(len(mds_result))])
        labels = np.asarray([f"{info[i].color_name}" for i in range(len(mds_result))])
        chrs = np.asarray([f"{info[i].chr_name}" for i in range(len(mds_result))])
        data = {'x': mds_result[:, 0], 'y': mds_result[:, 1], 'z': mds_result[:, 2], 'description': data_description,
                'labels': labels, 'chromosome': chrs}

        if representative_index is not None:
            data['labels'][representative_index] = "RS"

        df = pd.DataFrame(data)

        if coloring_type == 'NCBI':
            # Draw based on NCBI colors for each chromosome
            fig = px.scatter_3d(df, x='x', y='y', z='z', color='labels', symbol='chromosome',
                                color_discrete_map=AnnotationRecord.get_color_display_dict(),
                                hover_data={'x': False, 'y': False, 'z': False, 'labels': False, 'chromosome': True,
                                            'description': True})
        elif coloring_type == 'Chromosome':
            # Draw each Chromosome with different color
            fig = px.scatter_3d(df, x='x', y='y', z='z', color='chromosome', symbol='chromosome',
                                hover_data={'x': False, 'y': False, 'z': False, 'labels': False, 'chromosome': True,
                                            'description': True})
        else:
            raise ValueError("Invalid coloring type")

        axis_value = 0.15
        fig.update_layout(
            scene=dict(xaxis=dict(range=[-axis_value, axis_value]), yaxis=dict(range=[-axis_value, axis_value]),
                       zaxis=dict(range=[-axis_value, axis_value])))
        fig.show()


if __name__ == '__main__':
    representative = ChromosomeRepresentativeSelection('Aspergillus nidulans', 6, 'DSSIM', segment_length=500_000)
    x_r = 1  # math.ceil(ChromosomesHolder('Maize').get_largest_chromosome_length() / 500_000 / 20)
    # representative.plot_distance_variations("Whole Genome", plot_random_outliers=False)
    segment_length, threshold_list, apx_list, rand_list = [], [], [], []
    for chr_name in ["Whole Genome"]:  # representative.chromosomes_holder.get_all_chromosomes_name():
        representative.plot_distance_variations(chr_name, x_range=x_r, plot_random_outliers=False)
        dist_list, apx, rand = representative.get_approximation_error(chr_name, random_outliers=False)
        threshold_list.append(np.sum(dist_list < 0.24))
        segment_length.append(len(dist_list))
        apx_list.append(apx)
        # rand_list.append(rand)

    print(sum(threshold_list) / sum(segment_length))

    print(f"Average approximation error: {np.mean(apx_list)}")
    print(f"Average random error: {np.mean(rand_list)}")

    # chr_n = "21"
    # human_representative = ChromosomeRepresentativeSelection('Human', 9, 'DSSIM', segment_length=500_000)
    # segments_info = human_representative.get_non_overlapping_segments(chr_n)['segments_information']
    # RSSP = human_representative.get_representative(chr_n)['index']
    # ChromosomeRepresentativeSelection.plot_multi_dimensional_scaling(human_representative.get_distance_matrix(chr_n),
    #                                                                  segments_info, RSSP, coloring_type='NCBI')
