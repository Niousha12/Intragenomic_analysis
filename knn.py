import os
import pickle
import random

from matplotlib import pyplot as plt
from tqdm import tqdm

from chaos_game_representation import CGR
from chromosomes_holder import ChromosomesHolder
from constants import SCIENTIFIC_NAMES
from distances.distance_metrics import get_dist
from representative_selection import ChromosomeRepresentativeSelection
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

if __name__ == '__main__':
    # Create test dataset
    X_test = []
    y_test = []
    n_samples = 100
    sequence_length = 200_000
    species_list = ["Human", "Saccharomyces cerevisiae", "Maize", "Dictyostelium discoideum", "Pyrococcus furiosus",
                    "Escherichia coli"]

    X_pickle_path = f'outputs/knn/X_test_{sequence_length}_{n_samples}.pkl'
    y_pickle_path = f'outputs/knn/y_test_{sequence_length}_{n_samples}.pkl'

    Figures_path = f'Figures/knn/{sequence_length}_{n_samples}'

    if os.path.exists(X_pickle_path) and os.path.exists(y_pickle_path):
        with open(X_pickle_path, 'rb') as f:
            X_test = pickle.load(f)
        with open(y_pickle_path, 'rb') as f:
            y_test = pickle.load(f)
    else:
        for species in SCIENTIFIC_NAMES.keys():
            if species == "Paramecium caudatum" or species == "Aspergillus terreus":  # or species == "Chimp":
                continue
            print(f"Create test dataset for {species}...")
            chromosomes_holder = ChromosomesHolder(species)
            if chromosomes_holder.genome_length < 100_000_000:
                # rep_chr = "Whole Genome"
                seg_dict = chromosomes_holder.get_random_segments_list(sequence_length, n_samples,
                                                                       chromosome_name="Whole Genome", overlap=True)
            else:
                # rep_chr = chromosomes_holder.get_all_chromosomes_name()[0]
                seg_dict = {'segments_sequences': []}
                for i in range(n_samples):
                    random_segment = chromosomes_holder.get_random_segment(sequence_length,  # chromosome_name=rep_chr,
                                                                           remove_outlier=False, return_dict=False)
                    seg_dict['segments_sequences'].append(random_segment)

            seg_fcgr = [CGR(seg, 6).get_fcgr() for seg in tqdm(seg_dict['segments_sequences'])]
            X_test.extend(seg_fcgr)
            y_test.extend([species] * n_samples)

        if not os.path.exists('outputs/knn'):
            os.makedirs('outputs/knn')
        with open(X_pickle_path, 'wb') as f:
            pickle.dump(X_test, f)
        with open(y_pickle_path, 'wb') as f:
            pickle.dump(y_test, f)
        print("Turn of the random seed from chromosome_holder and representative_selection and "
              "run again with the same parameters.")
        quit()

    # Run the random t times to see the accuracy
    pipeline_accuracy = 0
    approximate_pipeline_accuracy = 0
    random_accuracy = []
    run_times = 52
    for t in range(run_times):
        # Create train dataset from representative
        X_train = []
        y_train = []
        for species in SCIENTIFIC_NAMES.keys():
            if species == "Paramecium caudatum" or species == "Aspergillus terreus":  # or species == "Chimp":
                continue
            # print(f"Creating train dataset for {species}...")
            chromosomes_holder = ChromosomesHolder(species)

            if t == 0:
                # This is using RSSP
                if chromosomes_holder.genome_length < 100_000_000:
                    rep_dict = ChromosomeRepresentativeSelection(species, 6, 'DSSIM', segment_length=sequence_length). \
                        get_representative("Whole Genome")
                else:
                    rep_dict = ChromosomeRepresentativeSelection(species, 6, 'DSSIM', segment_length=sequence_length). \
                        get_representative_of_representatives()
            elif t == 1:
                # This is using ARSSP
                if chromosomes_holder.genome_length < 100_000_000:
                    rep_dict = ChromosomeRepresentativeSelection(species, 6, 'DSSIM', segment_length=sequence_length). \
                        get_approximate_representative("Whole Genome")
                else:
                    rep_dict = ChromosomeRepresentativeSelection(species, 6, 'DSSIM', segment_length=sequence_length). \
                        get_representative_of_representatives(pipeline="ARSSP")
            else:
                # print("Starting the random representative selection...")
                if chromosomes_holder.genome_length < 100_000_000:
                    rep_chr = "Whole Genome"
                else:
                    rep_chr = random.choice(chromosomes_holder.get_all_chromosomes_name())
                # Get representative randomly
                rep_dict = ChromosomeRepresentativeSelection(species, 6, 'DSSIM', segment_length=sequence_length). \
                    get_random_representative(rep_chr)
            X_train.append(rep_dict['fcgr'])
            y_train.append(species)

        # Classifier
        y_pred = []
        for x in X_test:
            distances = []
            for x_train in X_train:
                distances.append(get_dist(x, x_train, 'DSSIM'))
            # Find y_pred
            y_pred_ = y_train[distances.index(min(distances))]
            y_pred.append(y_pred_)

        # Compute accuracy
        correct = 0
        for i in range(len(y_test)):
            if y_test[i] == y_pred[i]:
                correct += 1
        accuracy = correct / len(y_test) * 100
        print(f'Accuracy at {t}: {accuracy}')
        if t == 0:
            pipeline_accuracy = accuracy
            pipeline_used = "RSSP"
        elif t == 1:
            approximate_pipeline_accuracy = accuracy
            pipeline_used = "ARSSP"
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
        # Adjust layout to fit everything
        plt.tight_layout()

        plt.title(f'{pipeline_used} Representative - Acc {accuracy:.2f}')
        plt.xlabel('Predicted Label')
        plt.ylabel('True Label')
        if not os.path.exists(Figures_path):
            os.makedirs(Figures_path)
        plt.savefig(f'{Figures_path}/{pipeline_used}_{t}.png')
        # plt.show()
        plt.close()

    print(f'Pipeline Accuracy: {pipeline_accuracy}')
    print(f'Approximate Pipeline Accuracy: {approximate_pipeline_accuracy}')
    for acc in random_accuracy:
        if acc > pipeline_accuracy:
            print(f'Random accuracy was bigger than the Pipeline accuracy: {acc}')
    print(f'Random Accuracy: {sum(random_accuracy) / len(random_accuracy)}')
