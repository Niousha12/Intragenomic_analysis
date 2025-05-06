import math
import os
import pickle
from datetime import datetime

import numpy as np
import pandas as pd
from sklearn.manifold import MDS
from sklearn_extra.cluster import KMedoids
from matplotlib import pyplot as plt
from tqdm import tqdm

from chaos_game_representation import CGR
from chromosomes_holder import ChromosomesHolder, AnnotationRecord
from distances.distance_metrics import get_dist
import plotly.express as px

np.random.seed(24)


class ChromosomeRepresentativeSelection:
    def __init__(self, specie, kmer, distance_metric, segment_length=None, root_path='Data'):
        self.specie = specie
        self.chromosomes_holder = ChromosomesHolder(specie, root_path=root_path)
        self.kmer = kmer
        self.length = segment_length
        self.distance_metric = distance_metric
        self.x_range = None

        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        self.pickle_path_root = os.path.join(project_root, 'outputs', 'cache_pickles', specie)
        if not os.path.exists(self.pickle_path_root):
            os.makedirs(self.pickle_path_root)

    def get_representative_of_representatives(self, pipeline="RSSP"):
        representative_dict_list = []
        for chromosome_name in self.chromosomes_holder.get_all_chromosomes_name():
            if pipeline == "RSSP":
                representative_dict_list.append(self.get_representative(chromosome_name))
            elif pipeline == "ARSSP":
                representative_dict_list.append(self.get_approximate_representative(chromosome_name))
        representative_fcgrs = [rep['fcgr'] for rep in representative_dict_list]

        distance_matrix = np.zeros((len(representative_fcgrs), len(representative_fcgrs)))
        for i in range(distance_matrix.shape[0]):
            for j in range(i + 1):
                distance_matrix[i, j] = get_dist(representative_fcgrs[i], representative_fcgrs[j],
                                                 dist_m=self.distance_metric)
                distance_matrix[j, i] = distance_matrix[i, j]

        representative_of_representatives_index = self.find_centroid(distance_matrix, exclude_indices=None)
        representative_of_representatives = representative_dict_list[representative_of_representatives_index]

        return representative_of_representatives

    def get_representative(self, chromosome_name):
        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name)
        distance_matrix = self.get_distance_matrix(chromosome_name)

        centroid = self.find_centroid(distance_matrix, exclude_indices=None)

        return {"sequence": segments['segments_sequences'][centroid],
                "chromosome": chromosome_name,
                "index": centroid,
                "fcgr": fcgrs[centroid],
                "type": "non-approximative"}

    def get_random_representative_outlier_condition(self, chromosome_name, outlier=True):
        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name)

        random_centroid_index = self.chromosomes_holder.choose_random_fragment_index(chromosome_name, self.length,
                                                                                     outlier=outlier)

        return {"sequence": segments['segments_sequences'][random_centroid_index],
                "index": random_centroid_index,
                "fcgr": fcgrs[random_centroid_index],
                "type": f"random_outlier_{outlier}"}

    def get_random_representative(self, chromosome_name):
        segments = self.get_non_overlapping_segments(chromosome_name)
        fcgrs = self.get_fcgrs_of_segments(chromosome_name)

        random_centroid_index = np.random.randint(0, len(segments['segments_sequences']))

        return {"sequence": segments['segments_sequences'][random_centroid_index],
                "index": random_centroid_index,
                "fcgr": fcgrs[random_centroid_index],
                "type": f"random"}

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
            for i in range(distance_matrix.shape[0]):
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

    def get_approximation_error(self, chromosome_name, num_samples=30, random_outliers=True):
        # start_time = datetime.now()
        representative = self.get_representative(chromosome_name)
        # print(f"Time to find the actual representative: {datetime.now() - start_time}")
        distances_from_representative = self.get_distance_from_representative(chromosome_name, representative)

        start_time = datetime.now()
        approximative_representative = self.get_approximate_representative(chromosome_name, num_samples)
        end_time = datetime.now() - start_time
        distances_from_approximative_representative = \
            self.get_distance_from_representative(chromosome_name, approximative_representative)

        approximation_error = np.mean(
            np.abs((distances_from_approximative_representative - distances_from_representative)))
        # np.sqrt(np.mean((distances_from_approximative_representative - distances_from_representative) ** 2))

        # print(f"Approximation error for {num_samples} samples: {approximation_error}")

        random_error = None
        if random_outliers:
            random_representative = self.get_random_representative_outlier_condition(chromosome_name)
            distances_from_random_representative = \
                self.get_distance_from_representative(chromosome_name, random_representative)
            random_error = np.mean(
                np.abs((distances_from_random_representative - distances_from_representative)))
            # np.sqrt(np.mean((distances_from_random_representative - distances_from_representative) ** 2))
            # print(f"Random error: {random_error}")

        return distances_from_representative, approximation_error, random_error, end_time

    def plot_distance_variations(self, chromosome_name, plot_random_outliers=True, plot_approximate=True,
                                 random_sequences_number=30, x_range=None):
        prefix = "RepSeg"
        if plot_approximate:
            prefix = "RepSeg_and_aRepSeg"
            if plot_random_outliers:
                prefix = "RepSeg_and_aRepSeg_and_random"
        elif plot_random_outliers:
            prefix = "RepSeg_and_random"

        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        figure_path = os.path.join(project_root, 'Figures', 'Representative', self.specie, prefix)
        if not os.path.exists(figure_path):
            os.makedirs(figure_path)

        plt.figure(figsize=(10, 5))

        pipeline_representative_dict = self.get_representative(chromosome_name)
        distances_from_centroid = self.get_distance_from_representative(chromosome_name, pipeline_representative_dict)

        if plot_random_outliers:
            random_outlier_representative_dict = self.get_random_representative_outlier_condition(chromosome_name, True)
            distance_from_outlier_representative = \
                self.get_distance_from_representative(chromosome_name, random_outlier_representative_dict)
            plt.plot(distance_from_outlier_representative, marker='o', linestyle='-', markersize=4, color='black')

        if plot_approximate:
            approximate_representative_dict = self.get_approximate_representative(chromosome_name,
                                                                                  random_sequences_number)
            distance_from_approximate_representative = \
                self.get_distance_from_representative(chromosome_name, approximate_representative_dict)

            plt.plot(distance_from_approximate_representative, marker='o', linestyle='-', markersize=4, color='blue')

        plt.plot(distances_from_centroid, marker='o', linestyle='-', markersize=4, color='red')
        plt.grid(True)

        if chromosome_name == "Whole Genome":
            include_whole_genome = True
        else:
            include_whole_genome = False

        if self.x_range is None:
            if x_range is None:
                self.x_range = math.ceil(
                    self.chromosomes_holder.get_largest_chromosome_length(include_whole_genome) / self.length / 20)
            else:
                self.x_range = x_range

        x_ticks = []
        for i in range(0, int(self.x_range) + 1):
            x_ticks.append(i * 20)
        y_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        plt.xticks(x_ticks)
        plt.yticks(y_ticks)
        plt.savefig(f'{figure_path}/chr-{chromosome_name}_kmer-{self.kmer}_length-{self.length}.png')
        # plt.show()
        plt.close()

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
            for segment in tqdm(segments['segments_sequences']):
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
    def find_centroid(distance_matrix, mode="mean", exclude_indices=None):
        if exclude_indices is None:
            exclude_indices = []

        # Mask the excluded indices
        mask = np.ones(distance_matrix.shape[0], dtype=bool)
        mask[exclude_indices] = False
        # Apply the mask to rows and columns
        masked_matrix = distance_matrix[mask][:, mask]

        # Initialize scores with infinities to handle fully NaN rows
        scores = np.full(masked_matrix.shape[0], np.inf)
        # Identify rows that are not fully NaN
        valid_rows = ~np.isnan(masked_matrix).all(axis=1)
        # Compute relevant scores
        if mode == "mean" or mode == "median":
            scores[valid_rows] = np.nanmean(masked_matrix[valid_rows], axis=1)  # Compute mean
            if mode == "median":
                centroid_index = np.argsort(scores)[len(scores) // 2]  # Choose index
            else:
                centroid_index = np.argmin(scores)  # Choose index
        elif mode == "medoid":
            scores[valid_rows] = np.nansum(masked_matrix[valid_rows], axis=1)
            centroid_index = np.argmin(scores)  # Choose index
        elif mode == "kmedoid":
            # Replace NaNs with a large number so KMedoids can run
            sanitized_matrix = np.where(np.isnan(masked_matrix),
                                        np.nanmax(masked_matrix[np.isfinite(masked_matrix)]) * 10, masked_matrix)

            kmedoids = KMedoids(n_clusters=1, metric='precomputed', init='random', random_state=0)
            kmedoids.fit(sanitized_matrix)
            centroid_index = kmedoids.medoid_indices_[0]
        else:
            raise ValueError("mode must be 'mean', 'median', or 'medoid'")

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

        fig.update_traces(marker=dict(
            size=8,  # Increase marker size
            opacity=0.9,  # Adjust opacity for a shiny effect
            # line=dict(width=1, color='DarkSlateGrey'),  # Add a border around markers
            # color=df['labels'].apply(lambda x: AnnotationRecord.get_color_display_dict().get(x, 'grey'))
        ))

        axis_value = 1
        fig.update_layout(
            scene=dict(xaxis=dict(range=[-axis_value, axis_value], gridcolor="lightgrey",
                                  zerolinecolor="lightgrey", linecolor="lightgrey"),
                       yaxis=dict(range=[-axis_value, axis_value], gridcolor="lightgrey",
                                  zerolinecolor="lightgrey", linecolor="lightgrey"),
                       zaxis=dict(range=[-axis_value, axis_value], gridcolor="lightgrey",
                                  zerolinecolor="lightgrey", linecolor="lightgrey")))
        fig.show()


if __name__ == '__main__':
    representative = ChromosomeRepresentativeSelection('Human', 6, 'DSSIM', segment_length=500_000)

    x_r = math.ceil(ChromosomesHolder('Human').get_largest_chromosome_length() / 500_000 / 20)
    segment_length, threshold_list, apx_list, rand_list = [], [], [], []
    for chr_name in representative.chromosomes_holder.get_all_chromosomes_name():
        representative.plot_distance_variations(chr_name, x_range=x_r, plot_random_outliers=True)
        dist_list, _, _, _ = representative.get_approximation_error(chr_name, random_outliers=False)
        threshold_list.append(np.sum(dist_list < 0.24))
        segment_length.append(len(dist_list))
        # rand_list.append(rand)

    print(sum(threshold_list) / sum(segment_length))

    #
    # n_list = [1, 50, 40, 30, 20, 10]
    # for n in n_list:
    #     apx_list = []
    #     time_list = []
    #     for i in range(100):
    #         dist_list, apx, _, t = representative.get_approximation_error("1", num_samples=n, random_outliers=False)
    #         apx_list.append(apx)
    #         time_list.append(t)
    #     print(f"Average approximation error for {n} samples: {np.mean(apx_list)}")
    #     print(f"Average time for {n} samples: {np.mean(time_list)}")

    # x_r = 1  # math.ceil(ChromosomesHolder('Maize').get_largest_chromosome_length() / 500_000 / 20)
    # # representative.plot_distance_variations("Whole Genome", plot_random_outliers=False)
    # segment_length, threshold_list, apx_list, rand_list = [], [], [], []
    # for chr_name in ["Whole Genome"]:  # representative.chromosomes_holder.get_all_chromosomes_name():
    #     x_r = math.ceil(len(representative.chromosomes_holder.get_chromosome_sequence(chr_name)) / 500_000 / 20)
    #     representative.plot_distance_variations(chr_name, x_range=x_r, plot_random_outliers=False)
    #     dist_list, apx, rand = representative.get_approximation_error(chr_name, random_outliers=False)
    #     threshold_list.append(np.sum(dist_list < 0.24))
    #     segment_length.append(len(dist_list))
    #     apx_list.append(apx)
    #     # rand_list.append(rand)
    #
    # print(sum(threshold_list) / sum(segment_length))
    #
    # print(f"Average approximation error: {np.mean(apx_list)}")
    # print(f"Average random error: {np.mean(rand_list)}")

    # x_r = math.ceil(len(representative.chromosomes_holder.get_chromosome_sequence("1")) / 500_000 / 20)
    # # for chr_name in representative.chromosomes_holder.get_all_chromosomes_name():
    # representative.plot_distance_variations("1", x_range=x_r, plot_random_outliers=False, plot_approximate=False)

    # chr_n = "21"
    # human_representative = ChromosomeRepresentativeSelection('Human', 9, 'DSSIM', segment_length=500_000)
    # segments_info = human_representative.get_non_overlapping_segments(chr_n)['segments_information']
    # RSSP = human_representative.get_representative(chr_n)['index']
    # ChromosomeRepresentativeSelection.plot_multi_dimensional_scaling(human_representative.get_distance_matrix(chr_n),
    #                                                                  segments_info, RSSP, coloring_type='NCBI')
