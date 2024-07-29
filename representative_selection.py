import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

from chaos_game_representation import CGR
from chromosomes_holder import ChromosomesHolder
from distances.distance_metrics import get_dist


# TODO: Test the class and complete


class ChromosomeRepresentativeSelection:
    def __init__(self, specie, kmer, representative_length, distance_metric):
        self.specie = specie
        self.chromosomes_holder = ChromosomesHolder(specie)
        self.kmer = kmer
        self.length = representative_length
        self.distance_metric = distance_metric

        self.pickle_path_root = f"./outputs/cache_pickles/{self.specie}/"
        if not os.path.exists(self.pickle_path_root):
            os.makedirs(self.pickle_path_root)

    def get_representative(self, chromosome_name):
        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name, segments)
        distance_matrix = self.get_distance_matrix(chromosome_name, fcgrs)

        centroid = self.find_centroid(distance_matrix, exclude_indices=None)

        return {"representative_sequence": segments['segments_sequences'][centroid],
                "representative_index": centroid,
                "representative_fcgr": fcgrs[centroid]}

    def get_random_representative(self, chromosome_name, outlier=True):
        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name, segments)

        random_centroid_index = self.chromosomes_holder.choose_random_fragment_index(chromosome_name, self.length,
                                                                                     outlier=outlier)

        return {"representative_sequence": segments['segments_sequences'][random_centroid_index],
                "representative_index": random_centroid_index,
                "representative_fcgr": fcgrs[random_centroid_index]}

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

        return {"representative_sequence": representative_dict['sequence'],
                "representative_start": representative_dict['start'],
                "representative_fcgr": representative_dict['fcgr']}

    def plot_multi_dimensional_scaling(self):
        pass

    def plot_distance_variations(self, chromosome_name, plot_random_outliers=True, plot_approximate=True):
        figure_path = os.path.join('Figures', 'Representative', self.specie, 'Different_centroids')
        if not os.path.exists(figure_path):
            os.makedirs(figure_path)

        plt.figure(figsize=(10, 5))

        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name, segments)
        distance_matrix = self.get_distance_matrix(chromosome_name, fcgrs)

        pipeline_representative_dict = self.get_representative(chromosome_name)

        distances_from_centroid = distance_matrix[pipeline_representative_dict['representative_index'], :]

        if plot_random_outliers:
            random_outlier_representative_dict = self.get_random_representative(chromosome_name, True)
            distance_from_outlier_representative = distance_matrix[
                                                   random_outlier_representative_dict['representative_index'], :]
            plt.plot(distance_from_outlier_representative, marker='o', linestyle='-', markersize=4, color='purple')
        if plot_approximate:
            distance = np.zeros(len(fcgrs))
            apx_representative_fcgr = self.get_approximate_representative(chromosome_name)['representative_fcgr']
            for index, fcgr_matrix in enumerate(fcgrs):
                distance[index] = get_dist(fcgr_matrix, apx_representative_fcgr, dist_m=self.distance_metric)

            # Calculate the mean of the squared differences (MSE)
            # approximation_error = np.sqrt(np.mean((distance - distances_from_centroid) ** 2))

            plt.plot(distance, marker='o', linestyle='-', markersize=4, color='green')

        plt.plot(distances_from_centroid, marker='o', linestyle='-', markersize=4, color='red')
        plt.grid(True)

        x_ticks = []
        for i in range(0, 26):
            x_ticks.append(i * 20)
        y_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        plt.xticks(x_ticks)
        plt.yticks(y_ticks)
        plt.savefig(f'{figure_path}/distance_from_centroid_{chromosome_name}.png')

    def get_non_overlapping_segments(self, chromosome_name):
        pickle_path = os.path.join(self.pickle_path_root, f"chr_{chromosome_name}_len_{self.length}_segments.pickle")
        if os.path.exists(pickle_path):
            return self.load_pickle(pickle_path)
        else:
            segments = self.chromosomes_holder.get_chromosome_non_overlapping_segments(chromosome_name, self.length)
            self.create_pickle(segments, pickle_path)
            return segments

    def get_fcgrs_of_segments(self, chromosome_name, segments):
        pickle_path = os.path.join(self.pickle_path_root,
                                   f"chr_{chromosome_name}_len_{self.length}_kmer_{self.kmer}_fcgrs.pickle")
        if os.path.exists(pickle_path):
            return self.load_pickle(pickle_path)
        else:
            fcgrs_list = []
            for segment in segments['segments_sequences']:
                fcgrs_list.append(CGR(segment, self.kmer).get_fcgr())
            self.create_pickle(fcgrs_list, pickle_path)
            return fcgrs_list

    def get_distance_matrix(self, chromosome_name, fcgrs):
        pickle_path = os.path.join(self.pickle_path_root,
                                   f"chr_{chromosome_name}_len_{self.length}_kmer_{self.kmer}"
                                   f"_dist_{self.distance_metric}_distance_matrix.pickle")
        if os.path.exists(pickle_path):
            return self.load_pickle(pickle_path)
        else:
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
        sums = np.sum(masked_matrix, axis=1)

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


if __name__ == '__main__':
    human_representative = ChromosomeRepresentativeSelection('human', 6, 500_000, 'DSSIM')
    human_representative.plot_distance_variations('1')
