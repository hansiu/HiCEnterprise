"""
Script for plotting domain interaction maps.
Based on code by Irina Tuszynska and Rafal Zaborowski.
"""
import argparse
import os
import numpy as np
import csv
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ..utils import load_hicmap, create_folders, clip_and_blur


class Plotter:
    """
    Plots domain interaction maps with matplotlib.
    """

    def __init__(self, hic_map, stats_folder, interactions, chrom, threshold, plot_title, ticks_separ, hic_color,
                 inter_color, bin_res, distribution, hic_name, all_dom, inter_indom):
        hic_map = os.path.abspath(hic_map)
        stats_folder = os.path.abspath(stats_folder)
        self.chr = chrom
        hic_folder = os.path.dirname(os.path.abspath(hic_map))
        filename = os.path.basename(os.path.abspath(hic_map))
        self.hicmap = clip_and_blur(load_hicmap(hic_folder, filename))
        self.hic_name = hic_name or os.path.basename(hic_map).split('.')[0]
        self.threshold = threshold
        self.plot_title = plot_title
        self.ticks_separ = ticks_separ
        self.distribution = distribution
        self.interactions = self._get_interactions(stats_folder, interactions)
        self.interactions_name = os.path.basename(os.path.abspath(interactions))
        self.corr_interactions = self._get_interactions(stats_folder, interactions, corr="corr_")
        self.hic_color = hic_color
        self.inter_color = inter_color
        self.bin_res = bin_res
        self.all_dom = all_dom
        self.inter_indom = inter_indom

    def _get_interactions(self, stats_folder, filename, corr=""):
        name = stats_folder + '/' + filename + '-' + self.hic_name + '-' + corr + 'stats' + self.chr + '-' + \
               "_".join(str(self.threshold).split('.')) + '-' + self.distribution + '.txt'
        try:
            interactions_file = open(os.path.abspath(name), 'r')
            interactions = csv.reader(interactions_file, delimiter='\t')
            next(interactions)  # skip headers
            return interactions        
        except IOError:
            logger.error(name + " file is not accessible")

    def prepare_interaction_matrix(self, interactions):
        """
        Prepares the matrix that visualizes significant interactions between domains
        :param interactions: list of significant interactions
        """
        i_m = np.zeros((self.hicmap.shape[0], self.hicmap.shape[0]))
        for l in interactions:
            start1 = int(int(l[1]) / self.bin_res)
            end1 = int(int(l[2]) / self.bin_res) - 1
            start2 = int(int(l[3]) / self.bin_res)
            end2 = int(int(l[4]) / self.bin_res) - 1
            dom1 = (start1, end1)
            dom2 = (start2, end2)
            if dom1 == dom2:
                i_m[dom1[0]:dom1[1] + 1, dom2[0]:dom2[1] + 1] = 0.0
            else:
                if float(l[5]) != 0.0:
                    i_m[dom1[0]:dom1[1] + 1, dom2[0]:dom2[1] + 1] = round(-np.log10(float(l[5])), 5)
                else:
                    i_m[dom1[0]:dom1[1] + 1, dom2[0]:dom2[1] + 1] = 500  # some big number here

        return np.triu(i_m)

    def plot(self, interaction_matrix, figures_folder, corr=""):
        """
        Plots the interaction map - Hi-C in one triangle and interaction matrix in the other
        """
        plt.imshow(np.tril(self.hicmap), origin='lower', norm=LogNorm(), cmap=self.hic_color, interpolation='nearest')
        plt.imshow(interaction_matrix, origin='lower', norm=LogNorm(), cmap=self.inter_color, interpolation='nearest')
        plt.colorbar()
        plt.axis([0, self.hicmap.shape[0], 0, self.hicmap.shape[0]])
        len_ma = self.hicmap.shape[0]
        if self.ticks_separ != 0:
            plt.xticks(np.arange(0, len_ma, self.ticks_separ))
            plt.yticks(np.arange(0, len_ma, self.ticks_separ))
        else:
            pass
        plt.title(self.plot_title, fontsize=7)
        if self.all_dom == True:
            output = figures_folder + '/' + self.hic_name + '-' + corr + self.interactions_name.split('.')[
                0] + '-' + self.distribution + ".png"
        else: 
            output = figures_folder + '/' + self.hic_name + '-' + corr + self.interactions_name.split('.')[
                0] + '-' + self.distribution + "-" + "pericent_multipl_" + str(self.inter_indom) + ".png"
        plt.savefig(output, dpi=1500, bbox_inches='tight')
        plt.close()


    def run(self, figures_folder):
        """
        Runs the plotter
        """
        figures_folder = create_folders([figures_folder])[0]
        interaction_matrix = self.prepare_interaction_matrix(self.interactions)
        self.plot(interaction_matrix, figures_folder)
        corr_interaction_matrix = self.prepare_interaction_matrix(self.corr_interactions)
        self.plot(corr_interaction_matrix, figures_folder, corr="corr_")


# Argument Parsing
parser = argparse.ArgumentParser(description='Script for visualizing hicMaps with significant domain interactions from '
                                             'extract.py script')
parser.add_argument('-m', '--hic_map',
                    help='Hi-C map in numpy format', type=str, required=True)
parser.add_argument('-i', '--interactions', type=str,
                    help="Name of the interactions provided while extracting", required=True)
parser.add_argument('-c', '--chr', type=str, help="Chromosome number", required=True)
parser.add_argument('-s', '--stats_folder', help="Folder to load the significant interactions from", type=str,
                    default='./stats/')
parser.add_argument('-f', '--figures_folder', help="Folder to save the plots in", type=str, default='./figures/')
parser.add_argument('-t', '--threshold', type=float, help="Threshold that was used for statistical analysis. Default = 0.01", default=0.01)
parser.add_argument('-l', '--plot_title', type=str,
                    help="The title of the plot. If it contains spaces, use quotation marks.",
                    default='Interactions')
parser.add_argument('-e', '--ticks_separation', type=int, help="Frequency of ticks on the plot", default=0)
parser.add_argument('-o', '--hic_color', type=str, help="The color of HiC map, use your favorite from "
                                                        "https://matplotlib.org/api/pyplot_summary.html described as a "
                                                        "Colormap option.  Recommended: Reds, Blues,YlOrBr, PuBu. "
                                                        "Default is 'Greens'",
                    default='Greens')
parser.add_argument('-r', '--interactions_color', type=str, help="The color of interactions, use your favorite from "
                                                                 "https://matplotlib.org/api/pyplot_summary.html "
                                                                 "described as a Colormap option. Recommended: Reds, "
                                                                 "Blues, YlOrBr, PuBu. Default is 'YlOrBr'",
                    default='YlOrBr')
parser.add_argument('-b', '--bin_res', help='Resolution (size of the bins on the hicmaps) in bp i.e. 10000 for 10kb'
                                            ' resolution', type=int, required=True)
parser.add_argument('--distribution', type=str,
                    help="The distribution on which identification of domain-domain interactions was based."
                         " Available: hypergeom, negbinom, poisson. Default: hypergeom",
                    default='poisson')
parser.add_argument('-n', '--hic_name', help="Name to use for Hi-C map. default is the name of the file.",
                    type=str)
parser.add_argument('--all_domains', help="Stop remove pericentromeric domains and domains with rare interdomains contacts, where a mean number of contacts in one row is less than n*numer_of_domains. n = 1 and can be changed by --interact_indomain parameter",
                            action="store_true")
parser.add_argument('-g','--interact_indomain', type =float, help="Multiplier of domains_number, that is a threeshold for neglecting pericentromeric domains. It is not needed if --all_domains option is on. Default = 1.0, higher number - more domains will be removed",
                            default = 1)

# Main
if __name__ == "__main__":
    args = parser.parse_args()
    p = Plotter(args.hic_map, args.stats_folder, args.interactions, args.chr, args.threshold, args.plot_title,
                args.ticks_separation, args.hic_color, args.interactions_color, args.bin_res, args.distribution,
                args.hic_name)
    p.run(args.figures_folder)

# TODO add TAD plotting
