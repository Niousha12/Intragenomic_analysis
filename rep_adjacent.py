import os
import pickle

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
    n_samples = 100
    sequence_length = 500_000

    species = "Human"
    chromosomes_holder = ChromosomesHolder(species)
    rep_chr = chromosomes_holder.get_all_chromosomes_name()[0]

    # Create a list of random samples + the representative based on RSSP
    seg_dict = {'segments_sequences': []}
    rep_dict = ChromosomeRepresentativeSelection(species, 6, 'DSSIM', segment_length=sequence_length). \
        get_approximate_representative(rep_chr)
    for i in range(n_samples):
        random_segment = chromosomes_holder.get_random_segment(sequence_length, chromosome_name=rep_chr,
                                                               remove_outlier=False, return_dict=False)
        seg_dict['segments_sequences'].append(random_segment)

    seg_fcgr = [CGR(seg, 6).get_fcgr() for seg in tqdm(seg_dict['segments_sequences'])]

    # create the global FCGR
    print("Create the global FCGR...")
    global_FCGR = CGR(chromosomes_holder.get_chromosome_sequence(rep_chr), 6).get_fcgr()
    print("Global FCGR created.")

    # Calculate the distance between the representative and the random samples to the global FCGR
    dist_rep = get_dist(rep_dict['fcgr'], global_FCGR, 'DSSIM')
    print(f"Distance between the representative and the global FCGR: {dist_rep}")
    dist_samples = [get_dist(seg, global_FCGR, 'DSSIM') for seg in seg_fcgr]
    print(f"Max distance between the random samples and the global FCGR: {max(dist_samples)}")
    print(f"Min distance between the random samples and the global FCGR: {min(dist_samples)}")
    print(f"Mean distance between the random samples and the global FCGR: {sum(dist_samples) / len(dist_samples)}")
    # print(f"Distance between the random samples and the global FCGR: {dist_samples}")


