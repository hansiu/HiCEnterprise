import argparse
import os
import sys


def exec_regions(args):
    """
    execute regions analysis
    """
    from .regions.extract import Extractor
    e = Extractor(os.path.basename(args.region_file).split('.')[0], args.hic_folders, args.chr, args.bin_res,
                  args.num_bins,
                  args.threshold, args.pickled_folder, args.figures_folder, args.stats_folder)
    e.run(args.region_file, args.num_regs, args.section, args.plotting, args.single_sig, args.remap,
          args.regs_algorithm, args.stat_formats)


def exec_domains(args):
    """
    execute domains analysis
    """
    from .domains.extract import Extractor
    e = Extractor(args.domain_file, args.chr, args.bin_res, args.sherpa_lvl, args.hic_folder, args.threshold, args.plot_title, args.ticks_separation, args.hic_color, args.interactions_color)
    e.run(args.stats_folder, args.plotting, args.figures_folder)


# Main argument parser with arguments used in both types of analysis
parser = argparse.ArgumentParser(prog='HiCEnterprise')
parent_parser = argparse.ArgumentParser(add_help=False)
subparsers = parser.add_subparsers(title='Types of contact analysis',
                                   description='With what should we establish contact?',
                                   help='Whether you want to find contacts between regions or domains')
parent_parser.add_argument('-c', '--chr', help='Chromosome for which to extract domain interactions', type=str,
                           required=True)
parent_parser.add_argument('-b', '--bin_res',
                           help='Resolution (size of the bins on the hicmaps) in bp i.e. 10000 for 10kb'
                                ' resolution', type=int, required=True)
parent_parser.add_argument('-t', '--threshold', help='Statistical cutoff threshold',
                           type=float,
                           default=0.01)
parent_parser.add_argument('-s', '--stats_folder', help="Folder to save the statistics & significant points in",
                           type=str,
                           default='../stats/')
parent_parser.add_argument('-f', '--figures_folder', help="Folder to save the figures (plots) in", type=str,
                           default='../figures/')

# Regions analysis subparser
parser_regions = subparsers.add_parser('regions',
                                       help='Extracts region positions from EnhancerAtlas files or BED files, finds '
                                            'significant point interactions and generates interaction profile plots '
                                            'from Hi-C maps', parents=[parent_parser])
parser_regions.add_argument('-r', '--region_file',
                            help='Any tab-delimited BED file or FASTA File from EnhancerAtlas with '
                                 'regions and their positions.', type=str, required=True)
parser_regions.add_argument('--hic_folders',
                            help='Folder/folders in which the HiC data is stored with file names in format '
                                 'mtx-N-N.npy, where N = chromosome', type=str, nargs='+', required=True)
parser_regions.add_argument('-n', '--num_bins', help='Number of bins left and right to be considered', type=int,
                            default=100)
parser_regions.add_argument('--num_regs', help='Number of regions to consider and plot', type=int)
parser_regions.add_argument('-a', '--section', help="Section of bp for the plot", type=str)
parser_regions.add_argument('-p', '--pickled_folder', help="Folder with pickled files (to load from/save in)", type=str,
                            default='../pickles/')
parser_regions.add_argument('--plotting',
                            help="If there should be plotting, and if so should it be with rpy2 or matplotlib."
                                 " Options: mpl, rpy2.", type=str)
parser_regions.add_argument('--single_sig', help="If the single map significant points should be saved or not",
                            action='store_true')
parser_regions.add_argument('--remap',
                            help="Assembly that maps & regions are currently in, and assembly in which to save the "
                                 "stats in format assembly:new_assembly i.e. hg_19:hg_38", type=str)
parser_regions.add_argument('--regs_algorithm',
                            help="Algorithm for sorting regions to bins. Default 'all', available 'one'.",
                            type=str, default="all")
parser_regions.add_argument('--stat_formats',
                            help="In which file formats to save statistics. It is possible to provide "
                                 "multiple. Available: txt, bed, gff", type=str, nargs='+', default=['txt'])
parser_regions.set_defaults(func=exec_regions)

# Domains analysis subparser
parser_domains = subparsers.add_parser('domains', help='Extracts domains (TADs) interactions from Hi-C maps with given'
                                                       ' domain positions', parents=[parent_parser])
parser_domains.set_defaults(func=exec_domains)

parser_domains.add_argument('-m', '--hic_folder',
                            help='Folder in which the HiC data is stored with file names in format '
                                 'mtx-N-N.npy, where N = chromosome', type=str, required=True)
parser_domains.add_argument('-d', '--domain_file',
                            help="Txt file with domain definition in the format: dom_id(integer) "
                                 "chromosome(integer) dom_start(bp) dom_end(bp) sherpa-lvl(OPTIONAL)", type=str,
                            required=True)
parser_domains.add_argument('--plotting', help="Plot results. You'll need matplotlib library",
                            action="store_true")
parser_domains.add_argument('--sherpa_lvl', help="If there are sherpa levels in the file and which one to use",
                            type=int)
parser_domains.add_argument('-l', '--plot_title', type=str, help="The title of the plot",
                    default='Interactions')
parser_domains.add_argument('-e', '--ticks_separation', type=int, help="Frequency of ticks on the plot", default=0)
parser.add_argument('-o', '--hic_color', type=str, help="The color of HiC map, use your favorite from https://matplotlib.org/api/pyplot_summary.html described as a Colormap option.  Recommended: Reds, Blues, YlOrBr, PuBu. Default is 'Greens'", 
                    default='Greens')
parser.add_argument('-r', '--interactions_color', type=str, help="The color of HiC map, use your favorite from https://matplotlib.org/api/pyplot_summary.html described as a Colormap option. Recommended: Reds, Blues, YlOrBr, PuBu. Default is 'YlOrBr'",
                    default='YlOrBr')



def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    args = parser.parse_args()
    args.func(args)
