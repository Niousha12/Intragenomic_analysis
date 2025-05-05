import argparse
import os

from intergenomic_analysis import InterGenomicAnalysis
from intragenomic_analysis import IntraGenomicAnalysis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Experiment_type', type=str, required=True, help='There are two experiments: '
                                                                           'intragenomic and intergenomic')
    parser.add_argument('--root_path', type=str, default='../Data', help='Path to the datasets directory')
    parser.add_argument('--k_mer', type=int, default=6, help='K-mer size for FCGR')

    args = parser.parse_args()

    # Check if the experiment type is valid
    if args.Experiment_type not in ['intragenomic', 'intergenomic']:
        raise ValueError(f"Invalid experiment type: {args.Experiment_type}. "
                         f"Available options are: intragenomic, intergenomic")

    # Check if the k-mer size is valid
    if (args.k_mer < 1) or (args.k_mer > 9):
        raise ValueError(f"K-mer size {args.k_mer} is invalid. It should be between 1 and 9.")

    # # # # # # # # # # # #
    # Running the first experiment which is intragenomic
    # # # # # # # # # # # #
    if args.Experiment_type == 'intragenomic':
        species = 'Human'
        # Check if the genome file exists
        chromosomes_path = os.path.join(args.root_path, species)
        if not os.path.exists(chromosomes_path):
            raise FileNotFoundError(f"Genome file not found at {chromosomes_path}. Please check the path.")

        intragenome = IntraGenomicAnalysis(args.root_path, species, args.k_mer)
        dataframe = intragenome.run_experiment(new_run=True)
        intragenome.plot_intragenomic_analysis(dataframe)

    # # # # # # # # # # # #
    # Running the second experiment which is intergenomic
    # # # # # # # # # # # #
    elif args.Experiment_type == 'intergenomic':
        base_species = 'Human'
        target_list = ['Human', 'Chimp', 'Mouse', 'Drosophila melanogaster', 'Saccharomyces cerevisiae',
                       'Arabidopsis thailana', 'Paramecium caudatum', 'Pyrococcus furiosus', 'Escherichia coli']
        length = 500_000
        num_samples = 100

        # Check if the genome file exists
        chromosomes_path = os.path.join(args.root_path, base_species)
        if not os.path.exists(chromosomes_path):
            raise FileNotFoundError(f"Genome file not found at {chromosomes_path}. Please check the path.")

        # Check if the target species genome files exist
        for target in target_list:
            chromosomes_path = os.path.join(args.root_path, target)
            if not os.path.exists(chromosomes_path):
                raise FileNotFoundError(f"Genome file not found at {chromosomes_path}. Please check the path.")

        intergenome = InterGenomicAnalysis(root_path=args.root_path, base_specie=base_species,
                                           target_species_list=target_list, length=length, num_samples=num_samples,
                                           run=True)
        data_frame = intergenome.run_experiment(new_run=True)
        intergenome.plot_means_variances(data_frame, plot_type='boxplot')


if __name__ == '__main__':
    main()
