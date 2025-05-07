import argparse
import os
import pickle
import random

from matplotlib import pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from tqdm import tqdm

from chaos_game_representation import CGR
from chromosomes_holder import ChromosomesHolder
from constants import SCIENTIFIC_NAMES, DISTANCE_METRICS_LIST
from distances.distance_metrics import get_dist
from representative_selection import ChromosomeRepresentativeSelection


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Experiment_type', type=str, default='all', help='Two options are available:'
                                                                                          'all or no_chimp.'
                                                                                          'If you want to run the '
                                                                                          'experiment with all the '
                                                                                          'species then set it to all. '
                                                                                          'If you want to run the '
                                                                                          'experiment without chimp '
                                                                                          'then set it to no_chimp.')
    parser.add_argument('--root_path', type=str, default='Data', help='Path to the datasets directory')

    parser.add_argument('--number_of_samples', type=int, default=100, help='Number of samples from each species for '
                                                                           'creating test set')
    parser.add_argument('--segment_length', type=int, default=200_000, help='Length of the representative segment'
                                                                            'as well as test segments')
    parser.add_argument('--k_mer', type=int, default=6, help='K-mer size for FCGR')
    parser.add_argument('--distance_metric', type=str, default='DSSIM', help='Distance metric comparison')
    parser.add_argument('--run_times', type=int, default=50, help='Number of times to run the random experiment')

    args = parser.parse_args()

    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    # Check if the genome file exists
    species_list = list(SCIENTIFIC_NAMES.keys())
    to_remove = ['Paramecium caudatum', 'Aspergillus terreus']
    if args.Experiment_type == 'all':
        species_list = [s for s in species_list if s not in to_remove]
    elif args.Experiment_type == 'no_chimp':
        to_remove.append('Chimp')
        species_list = [s for s in species_list if s not in to_remove]
    else:
        raise ValueError(f"Invalid species list: {args.Experiment_type}. It should be either 'all' or 'no_chimp'.")
    for species in species_list:
        chromosomes_path = os.path.join(project_root, args.root_path, species)
        if not os.path.exists(chromosomes_path):
            raise FileNotFoundError(f"Genome file not found at {chromosomes_path}. Please check the path.")

    # Check if the number of samples is valid
    if args.number_of_samples < 1:
        raise ValueError(f"Number of samples {args.number_of_samples} is invalid. It should be greater than 0.")

    # Check if the segment length is valid
    if args.segment_length < 1:
        raise ValueError(f"Segment length {args.segment_length} is invalid. It should be greater than 0.")

    # Check if the k-mer size is valid
    if (args.k_mer < 1) or (args.k_mer > 9):
        raise ValueError(f"K-mer size {args.k_mer} is invalid. It should be between 1 and 9.")

    # Check if the distance metric is valid
    if args.distance_metric not in DISTANCE_METRICS_LIST:
        raise ValueError(f"Invalid distance metric: {args.distance_metric}. "
                         f"Available options are: {', '.join(DISTANCE_METRICS_LIST)}")

    # # # # # # # # # # # #
    # Create the test set
    # # # # # # # # # # # #
    X_test = []
    y_test = []

    pickle_path = os.path.join(project_root, 'outputs', 'knn')
    X_pickle_path = os.path.join(pickle_path,
                                 f'X_test_{args.segment_length}_{args.number_of_samples}_{args.Experiment_type}.pkl')
    y_pickle_path = os.path.join(pickle_path,
                                 f'y_test_{args.segment_length}_{args.number_of_samples}_{args.Experiment_type}.pkl')
    figures_path = os.path.join(project_root, 'Figures', 'knn',
                                f'{args.segment_length}_{args.number_of_samples}_{args.Experiment_type}')

    if os.path.exists(X_pickle_path) and os.path.exists(y_pickle_path):
        with open(X_pickle_path, 'rb') as f:
            X_test = pickle.load(f)
        with open(y_pickle_path, 'rb') as f:
            y_test = pickle.load(f)
    else:
        for species in species_list:
            print(f"Create test dataset for {species}...")
            chromosomes_holder = ChromosomesHolder(species, args.root_path)
            if chromosomes_holder.genome_length < 100_000_000:
                seg_dict = chromosomes_holder.get_random_segments_list(args.segment_length, args.number_of_samples,
                                                                       chromosome_name="Whole Genome", overlap=True)
            else:
                seg_dict = {'segments_sequences': []}
                for i in range(args.number_of_samples):
                    random_segment = chromosomes_holder.get_random_segment(args.segment_length, remove_outlier=False,
                                                                           return_dict=False)
                    seg_dict['segments_sequences'].append(random_segment)

            seg_fcgr = [CGR(seg, args.k_mer).get_fcgr() for seg in tqdm(seg_dict['segments_sequences'])]
            X_test.extend(seg_fcgr)
            y_test.extend([species] * args.number_of_samples)

        if not os.path.exists(pickle_path):
            os.makedirs(pickle_path)
        with open(X_pickle_path, 'wb') as f:
            pickle.dump(X_test, f)
        with open(y_pickle_path, 'wb') as f:
            pickle.dump(y_test, f)

    # # # # # # # # # # # #
    # Running the experiment
    # # # # # # # # # # # #
    pipeline_accuracy = 0
    approximate_pipeline_accuracy = 0
    random_accuracy = []
    run_times = args.run_times + 2
    for t in tqdm(range(run_times)):
        # Create train dataset from representative
        X_train = []
        y_train = []
        for species in species_list:
            # print(f"Creating train dataset for {species}...")
            chromosomes_holder = ChromosomesHolder(species, args.root_path)
            representative_selection = ChromosomeRepresentativeSelection(species, args.k_mer, args.distance_metric,
                                                                         segment_length=args.segment_length,
                                                                         root_path=args.root_path)
            if t == 0:
                # This is using RepSeg
                if chromosomes_holder.genome_length < 100_000_000:
                    rep_dict = representative_selection.get_representative("Whole Genome")
                else:
                    rep_dict = representative_selection.get_representative_of_representatives()
            elif t == 1:
                # This is using aRepSeg
                if chromosomes_holder.genome_length < 100_000_000:
                    rep_dict = representative_selection.get_approximate_representative("Whole Genome")
                else:
                    rep_dict = representative_selection.get_representative_of_representatives(pipeline="RepSeg")
            else:
                # This is using Random
                if chromosomes_holder.genome_length < 100_000_000:
                    rep_chr = "Whole Genome"
                else:
                    rep_chr = random.choice(chromosomes_holder.get_all_chromosomes_name())
                # Get representative randomly
                rep_dict = representative_selection.get_random_representative(rep_chr)
            X_train.append(rep_dict['fcgr'])
            y_train.append(species)

        # Classifier
        y_pred = []
        for x in X_test:
            distances = []
            for x_train in X_train:
                distances.append(get_dist(x, x_train, args.distance_metric))
            # Find y_pred
            y_pred_ = y_train[distances.index(min(distances))]
            y_pred.append(y_pred_)

        # Compute accuracy
        correct = 0
        for i in range(len(y_test)):
            if y_test[i] == y_pred[i]:
                correct += 1
        accuracy = correct / len(y_test) * args.number_of_samples
        print(f'Accuracy at {t}: {accuracy}')
        if t == 0:
            pipeline_accuracy = accuracy
            pipeline_used = "RepSeg"
        elif t == 1:
            approximate_pipeline_accuracy = accuracy
            pipeline_used = "aRepSeg"
        else:
            random_accuracy.append(accuracy)
            pipeline_used = "Random"

        # Step 1: Get all unique classes
        all_labels = sorted(set(y_test + y_pred))

        # Step 2: Map class names to integers
        label_to_int = {label: idx for idx, label in enumerate(all_labels)}
        int_to_label = {idx: label for label, idx in label_to_int.items()}  # optional reverse map

        # Step 3: Convert to integer labels
        y_test_int = [label_to_int[label] for label in y_test]
        y_pred_int = [label_to_int[label] for label in y_pred]

        # Compute confusion matrix
        cm = confusion_matrix(y_test_int, y_pred_int)

        plt.figure(figsize=(7, 6))  # Adjust the size as needed
        # disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
        disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                      display_labels=[int_to_label[i] for i in range(len(all_labels))])
        disp.plot(cmap='plasma', ax=plt.gca())  # You can change the color map 'plasma', 'viridis' or 'Blues'
        plt.xticks(rotation=45, ha='right')  # Rotate x labels (Predicted)
        plt.tight_layout()

        plt.title(f'{pipeline_used} Representative - Acc {accuracy:.2f}')
        plt.xlabel('Predicted Label')
        plt.ylabel('True Label')
        if not os.path.exists(figures_path):
            os.makedirs(figures_path)
        plt.savefig(f'{figures_path}/{pipeline_used}_{t}.png', dpi=300, bbox_inches='tight', transparent=True)
        # plt.show()
        plt.close()

    print(f'RepSeg accuracy: {pipeline_accuracy}')
    print(f'aRepSeg accuracy: {approximate_pipeline_accuracy}')
    print(f'Average random accuracy over {args.run_times} runs: {sum(random_accuracy) / len(random_accuracy)}')


if __name__ == '__main__':
    main()
