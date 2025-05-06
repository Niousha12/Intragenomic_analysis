import argparse
import logging
import os

import numpy as np
from tqdm import tqdm

from chromosomes_holder import ChromosomesHolder
from constants import DISTANCE_METRICS_LIST
from representative_selection import ChromosomeRepresentativeSelection


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', type=str, required=True, help='Species name')
    parser.add_argument('--root_path', type=str, default='Data', help='Path to the datasets directory')

    parser.add_argument('--k_mer', type=int, default=6, help='K-mer size for FCGR')
    parser.add_argument('--segment_length', type=int, default=500_000, help='Length of the representative segment')
    parser.add_argument('--distance_metric', type=str, default='DSSIM', help='Distance metric for representative '
                                                                             'segment selection')
    parser.add_argument('--chromosome_name', type=str, default='all', help='Chromosome name to plot, if it is set to '
                                                                           'all, it will plot all the chromosomes. '
                                                                           'In case that the length of '
                                                                           'the whole genome is less than 100_000_000 '
                                                                           'the representative will be chosen from the '
                                                                           'whole genome and this parameter '
                                                                           'will be ignored.')
    parser.add_argument('--plot_approximate', type=bool, default=True, help='Plot Experiment 3.1 and 3.2')
    parser.add_argument('--plot_random_outliers', type=bool, default=False, help='Plot Experiment 3.1 and 3.3. This is '
                                                                                 'only available for human genome')
    parser.add_argument('--n_value', type=int, default=30, help='Length of the set of random sequences '
                                                                'in Experiment 3.2.')
    parser.add_argument('--get_MAE_aRepSeg', type=bool, default=False, help='Get the average approximation error '
                                                                            'for aRepSeg')
    parser.add_argument('--get_MAE_random_outliers', type=bool, default=False, help='Get the average approximation '
                                                                                    'error for random outliers. '
                                                                                    'This is only available for human')
    parser.add_argument('--get_threshold', type=bool, default=False, help='Get the average number of segments '
                                                                          'that are below the threshold value')
    parser.add_argument('--threshold_value', type=float, default=0.24, help='Threshold value for the distance '
                                                                            'between segments')
    parser.add_argument('--plot_MDS', type=bool, default=False, help='Plot the Multi dimensional scaling')

    # --species "Human" --root_path "Data" --k_mer 6 --segment_length 500000 --distance_metric "DSSIM" --chromosome_name "all" --plot_approximate True --plot_random_outliers True --n_value 30 --get_MAE_aRepSeg True --get_MAE_random_outliers True --get_threshold True --threshold_value 0.24 --plot_MDS True

    args = parser.parse_args()

    # Check if the genome file exists
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    chromosomes_path = os.path.join(project_root, args.root_path, args.species)
    if not os.path.exists(chromosomes_path):
        raise FileNotFoundError(f"Genome file not found at {chromosomes_path}. Please check the path.")
    chromosomes_holder = ChromosomesHolder(args.species, args.root_path)

    # Check if the chromosome name is valid
    if chromosomes_holder.genome_length > 100_000_000:
        if (args.chromosome_name != 'all') and (
                args.chromosome_name not in chromosomes_holder.get_all_chromosomes_name()):
            raise ValueError(f"Invalid chromosome name: {args.chromosome_name}. "
                             f"Available chromosomes are: {chromosomes_holder.get_all_chromosomes_name()}")
    else:
        args.chromosome_name = 'Whole Genome'
        logging.warning(f"The genome length for {args.species} species is less than 100_000_000, "
                        f"so the chromosome name is ignored and the whole genome will be used.")

    # Check if the k-mer size is valid
    if (args.k_mer < 1) or (args.k_mer > 9):
        raise ValueError(f"K-mer size {args.k_mer} is invalid. It should be between 1 and 9.")

    # Check if the segment length is valid
    if (args.chromosome_name is not None) and (args.segment_length is not None):
        if chromosomes_holder.genome_length > 100_000_000:
            if args.chromosome_name != 'all':
                if args.segment_length > len(chromosomes_holder.get_chromosome_sequence(args.chromosome_name)):
                    raise ValueError(f"Segment length {args.segment_length} is greater than chromosome length "
                                     f"{len(chromosomes_holder.get_chromosome_sequence(args.chromosome_name))}")
            else:
                for chr_name in chromosomes_holder.get_all_chromosomes_name():
                    if args.segment_length > len(chromosomes_holder.get_chromosome_sequence(chr_name)):
                        raise ValueError(f"Segment length {args.segment_length} is greater than chromosome length "
                                         f"{len(chromosomes_holder.get_chromosome_sequence(chr_name))} in "
                                         f"chromosome {chr_name}")
        else:
            if args.segment_length > chromosomes_holder.genome_length:
                raise ValueError(f"Segment length {args.segment_length} is greater than genome length "
                                 f"{chromosomes_holder.genome_length}")

    # Check if the distance metric is valid
    if args.distance_metric not in DISTANCE_METRICS_LIST:
        raise ValueError(f"Invalid distance metric: {args.distance_metric}. "
                         f"Available options are: {', '.join(DISTANCE_METRICS_LIST)}")

    # Check if the plot_approximate is valid
    if args.plot_approximate not in [True, False]:
        raise ValueError(f"plot_approximate {args.plot_approximate} is invalid. It should be either True or False.")

    # Check if the plot_random_outliers is valid
    if args.species == 'Human':
        if args.plot_random_outliers not in [True, False]:
            raise ValueError(
                f"plot_random_outliers {args.plot_random_outliers} is invalid. It should be either True or False.")
    else:
        if args.plot_random_outliers:
            raise ValueError(f"plot_random_outliers {args.plot_random_outliers} is invalid. "
                             f"It should be False for species other than Human.")

    # Check if the n_value is valid
    if args.n_value < 1:
        raise ValueError(f"n_value {args.n_value} is invalid. It should be greater than 0.")

    # Check if the get_MAE_aRepSeg is valid
    if args.get_MAE_aRepSeg not in [True, False]:
        raise ValueError(f"get_MAE_aRepSeg {args.get_MAE_aRepSeg} is invalid. It should be either True or False.")

    # Check if the get_MAE_random_outliers is valid
    if args.species == 'Human':
        if args.get_MAE_random_outliers not in [True, False]:
            raise ValueError(f"get_MAE_random_outliers {args.get_MAE_random_outliers} is invalid. "
                             f"It should be either True or False.")
    else:
        if args.get_MAE_random_outliers:
            raise ValueError(f"get_MAE_random_outliers {args.get_MAE_random_outliers} is invalid. "
                             f"It should be False for species other than Human.")

    # Check if the get_threshold is valid
    if args.get_threshold not in [True, False]:
        raise ValueError(f"get_threshold {args.get_threshold} is invalid. It should be either True or False.")

    # Check if the threshold value is valid
    if args.get_threshold and (args.threshold_value < 0):
        raise ValueError(f"threshold_value {args.threshold_value} is invalid. It should be greater than 0.")
    if args.get_threshold and (args.threshold_value > 1):
        raise ValueError(f"threshold_value {args.threshold_value} is invalid. It should be less than 1.")
    if args.get_threshold:
        logging.warning(f"Threshold value is {args.threshold_value}.")

    # Check if the plot_MDS is valid
    if args.plot_MDS not in [True, False]:
        raise ValueError(f"plot_MDS {args.plot_MDS} is invalid. It should be either True or False.")

    representative = ChromosomeRepresentativeSelection(args.species, args.k_mer, args.distance_metric,
                                                       args.segment_length, args.root_path)

    if args.chromosome_name == 'Whole Genome':
        chromosomes_list = ["Whole Genome"]
    elif args.chromosome_name == 'all':
        chromosomes_list = representative.chromosomes_holder.get_all_chromosomes_name()
    else:
        chromosomes_list = [args.chromosome_name]

    approximate_error_list, rand_error_list = [], []
    segment_length, threshold_list = [], []
    for chr_name in tqdm(chromosomes_list):
        # for each chromosome plot the distance
        representative.plot_distance_variations(chr_name, plot_random_outliers=args.plot_random_outliers,
                                                plot_approximate=args.plot_approximate,
                                                random_sequences_number=args.n_value)

        # Calculate the MAE error or count the number of segments below the threshold
        if args.get_MAE_aRepSeg or args.get_MAE_random_outliers or args.get_threshold:
            distances, approximation_error, random_error, approximation_time = \
                representative.get_approximation_error(chr_name, num_samples=args.n_value,
                                                       random_outliers=args.plot_random_outliers)
            approximate_error_list.append(approximation_error)
            rand_error_list.append(random_error)
            segment_length.append(len(distances))
            threshold_list.append(np.sum(distances < args.threshold_value))

        # plot MDS
        if args.plot_MDS:
            distance_matrix = representative.get_distance_matrix(chr_name)
            segment_informations = representative.get_non_overlapping_segments(chr_name)['segments_information']
            representative_index = representative.get_representative(chr_name)['index']

            representative.plot_multi_dimensional_scaling(distance_matrix, segment_informations, representative_index)

    if args.get_MAE_aRepSeg:
        print(f"Average approximation error: {np.mean(approximate_error_list)}")
    if args.get_MAE_random_outliers:
        print(f"Average random error: {np.mean(rand_error_list)}")
    if args.get_threshold:
        print(f"Average number of segments below the {args.threshold_value} threshold: "
              f"{sum(threshold_list) / sum(segment_length)}")


if __name__ == '__main__':
    main()
