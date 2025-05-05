import argparse
import os

from chromosomes_holder import ChromosomesHolder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', type=str, required=True, help='Species name')
    parser.add_argument('--root_path', type=str, default='../Data', help='Path to the datasets directory')

    parser.add_argument('--chromosome_name', type=str, default='all', help='Set it to the name of the chromosome that '
                                                                           'you want to plot, if you want to '
                                                                           'plot all chromosomes then set all')
    parser.add_argument('--start_of_segment', type=int, default=None, help='Start of the segment to be plotted')
    parser.add_argument('--segment_length', type=int, default=None, help='Length of the segment to be plotted')
    parser.add_argument('--k_mer', type=int, default=9, help='K-mer size for FCGR')
    parser.add_argument('--fcgr_cgr', type=str, default='fcgr', help='Type of plot: FCGR or cgr')
    parser.add_argument('--label', type=bool, default=True, help='Label the plot with ACGT, '
                                                                 'currently it is only available for k=6 and k=9')
    parser.add_argument('--global_normalization', type=bool, default=False, help='Normalize the FCGR plot using '
                                                                                 'the min and max values of '
                                                                                 'all FCGRs from the chromosomes '
                                                                                 'of the genome species')
    parser.add_argument('--plot_3d', type=bool, default=False, help='plot in 3D')
    parser.add_argument('--resolution', type=int, default=None, help='Resolution of the plot')
    parser.add_argument('--bits', type=int, default=None, help='Bits for the plot')

    args = parser.parse_args()

    # Check if the genome file exists
    chromosomes_path = os.path.join(args.root_path, args.species)
    if not os.path.exists(chromosomes_path):
        raise FileNotFoundError(f"Genome file not found at {chromosomes_path}. Please check the path.")

    chromosomes_holder = ChromosomesHolder(args.species, args.root_path)

    # Check if the chromosome name is valid
    if (args.chromosome_name != 'all') and (args.chromosome_name not in chromosomes_holder.get_all_chromosomes_name()):
        raise ValueError(f"Invalid chromosome name: {args.chromosome_name}. "
                         f"Available chromosomes are: {chromosomes_holder.get_all_chromosomes_name()}")

    # Check if the segment length is valid
    if (args.chromosome_name is not None) and (args.segment_length is not None):
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

    # Check if the start of segment is valid
    if (args.chromosome_name is not None) and (args.start_of_segment is not None):
        if args.chromosome_name != 'all':
            if args.start_of_segment > len(chromosomes_holder.get_chromosome_sequence(args.chromosome_name)):
                raise ValueError(f"Start of segment {args.start_of_segment} is greater than chromosome length "
                                 f"{len(chromosomes_holder.get_chromosome_sequence(args.chromosome_name))}")
        else:
            for chr_name in chromosomes_holder.get_all_chromosomes_name():
                if args.start_of_segment > len(chromosomes_holder.get_chromosome_sequence(chr_name)):
                    raise ValueError(f"Start of segment {args.start_of_segment} is greater than chromosome length "
                                     f"{len(chromosomes_holder.get_chromosome_sequence(chr_name))} in "
                                     f"chromosome {chr_name}")

    # Check if the k-mer size is valid
    if (args.k_mer < 1) or (args.k_mer > 9):
        raise ValueError(f"K-mer size {args.k_mer} is invalid. It should be between 1 and 9.")

    # Check if the fcgr_cgr is valid
    if args.fcgr_cgr not in ['fcgr', 'cgr']:
        raise ValueError(f"fcgr_cgr {args.fcgr_cgr} is invalid. It should be either 'fcgr' or 'cgr'.")

    # Check if the label is valid
    if args.label not in [True, False]:
        raise ValueError(f"label {args.label} is invalid. It should be either True or False.")

    # Check if the plot_3d is valid
    if args.plot_3d not in [True, False]:
        raise ValueError(f"plot_3d {args.plot_3d} is invalid. It should be either True or False.")

    # Check if the resolution is valid
    if (args.resolution is not None) and (args.resolution < 1):
        raise ValueError(f"Resolution {args.resolution} is invalid. It should be greater than 0.")

    # Check if the bits is valid
    if (args.bits is not None) and (args.bits < 1):
        raise ValueError(f"Bits {args.bits} is invalid. It should be greater than 0.")

    # Check if the global normalization is valid
    if args.global_normalization not in [True, False]:
        raise ValueError(f"global_normalization {args.global_normalization} is invalid. "
                         f"It should be either True or False.")

    global_min = None  # For human this is 0, for maize is 21
    global_max = None  # For human this is 1134132, for maize is 246147
    if args.global_normalization:
        # Get the global min and max values from the fcgrs of all chromosomes
        global_min, global_max = chromosomes_holder.get_global_min_max(args.start_of_segment, args.segment_length,
                                                                       args.k_mer)

    if args.chromosome_name == 'all':
        for chr_name in chromosomes_holder.get_all_chromosomes_name():
            print(f"Plotting {args.fcgr_cgr} for chromosome {chr_name}...")
            chromosomes_holder.plot_fcgr(chr_name, start_of_segment=args.start_of_segment,
                                         segment_length=args.segment_length, k_mer=args.k_mer, fcgr_cgr=args.fcgr_cgr,
                                         label=args.label, global_min=global_min, global_max=global_max,
                                         _3d=args.plot_3d, resolution=args.resolution, bits=args.bits)
    else:
        print(f"Plotting {args.fcgr_cgr} for chromosome {args.chromosome_name}...")
        chromosomes_holder.plot_fcgr(args.chromosome_name, start_of_segment=args.start_of_segment,
                                     segment_length=args.segment_length, k_mer=args.k_mer, fcgr_cgr=args.fcgr_cgr,
                                     label=args.label, global_min=global_min, global_max=global_max,
                                     _3d=args.plot_3d, resolution=args.resolution, bits=args.bits)


if __name__ == '__main__':
    main()
