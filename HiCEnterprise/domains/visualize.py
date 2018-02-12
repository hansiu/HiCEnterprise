"""
Script for plotting domain interaction maps.
Based on code by Irina Tuszynska and Rafal Zaborowski.
"""
import argparse, os, scipy.ndimage
import numpy as np
import matplotlib
import seaborn as sns
import csv

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ..utils import load_hicmap, create_folders


class Plotter:
    """
    Plots domain interaction maps with matplotlib.
    """

    def __init__(self, hic_folder, stats_folder, interactions, chrom, threshold, plot_title):
        hic_folder = os.path.abspath(hic_folder)
        stats_folder = os.path.abspath(stats_folder)
        self.chr = chrom
        self.hicmap = load_hicmap(hic_folder, 'mtx-' + self.chr + '-' + self.chr + '.npy')
        self.hic_name = os.path.basename(hic_folder)
        self.threshold = threshold
        self.plot_title = plot_title
        self.interactions = self._get_interactions(stats_folder, interactions)
        self.interactions_name = os.path.basename(os.path.abspath(interactions))
        self.corr_interactions = self._get_interactions(stats_folder, interactions, corr="corr_")

    def _get_interactions(self, stats_folder, filename, corr=""):
        name = stats_folder + '/' + filename + '-' + self.hic_name + '-' + corr + 'stats' + self.chr + '-' + \
               "_".join(str(self.threshold).split('.')) + '.txt'
        interactions_file = open(os.path.abspath(name), 'r')
        interactions = csv.reader(interactions_file, delimiter='\t')
        return (interactions)

    def prepare_interaction_matrix(self, interactions):
        """
        Prepares the matrix that visualizes significant interactions between domains
        :param interactions: list of significant interactions
        """
        i_m = np.zeros((self.hicmap.shape[0], self.hicmap.shape[0]))
        for l in interactions:
            start1 = float(l[1])
            end1 = (float(l[2])) - 1
            start2 = float(l[3])
            end2 = (float(l[4])) - 1
            dom1 = (start1, end1)
            dom2 = (start2, end2)
            if dom1 == dom2:
                i_m[int(dom1[0]):int(dom1[1]) + 1, int(dom2[0]):int(dom2[1]) + 1] = 0.0
            else:
                i_m[int(dom1[0]):int(dom1[1]) + 1, int(dom2[0]):int(dom2[1]) + 1] = -np.log10(float(l[5]))

        return np.triu(i_m)
        
    def clip_and_blur(self,arr, stddevs=5, blur=1):
        #print 'PRINT', arr
        arr = np.ma.masked_invalid(arr)
        mean = np.mean(arr)
        stddev = np.var(arr) ** 0.5
        np.clip(arr, 0, mean + stddevs * stddev, out=arr)
        arr = np.ma.filled(arr, 0)
        scipy.ndimage.gaussian_filter(arr, blur, output=arr)
        np.clip(arr, mean * 0.01, mean + stddevs * stddev, out=arr)
        return arr

    def plot(self, interaction_matrix, figures_folder, corr=""):
        """
        Plots the interaction map - Hi-C in one triangle and interaction matrix in the other
        """
        sns.set_style("ticks")
        sns.despine(right=True)
        hicmap= self.clip_and_blur(self.hicmap)
        plt.imshow(np.tril(hicmap), origin='lower', norm=LogNorm(), cmap="Blues", interpolation='nearest')
        plt.imshow(interaction_matrix, origin='upper', norm=LogNorm(), cmap="Reds", interpolation='nearest')
        plt.colorbar()
        plt.axis([0, self.hicmap.shape[0], 0, self.hicmap.shape[0]])
        len_ma = self.hicmap.shape[0]
        plt.xticks(np.arange(0, len_ma, 400))
        plt.title(self.plot_title, fontsize=7)
        #ax = sns.heatmap(np.tril(self.hicmap), cmap="Reds", cbar = True)
        output = figures_folder + '/' + self.hic_name + '-' + corr + self.interactions_name.split('.')[0] + ".png"
        plt.savefig(output, dpi=1500, bbox_inches='tight')
        plt.close()

    def run(self, figures_folder):
        """
        Runs the plotter
        """
        interaction_matrix = self.prepare_interaction_matrix(self.interactions)
        self.plot(interaction_matrix, figures_folder)
        corr_interaction_matrix = self.prepare_interaction_matrix(self.corr_interactions)
        self.plot(corr_interaction_matrix, figures_folder, corr="corr_")


# Argument Parsing
parser = argparse.ArgumentParser(description='Script for visualizing hicMaps with significant domain interactions from '
                                             'extract.py script')
parser.add_argument('--hic_folder', type=str, help='Folder in which the HiC data is stored with file names in format '
                                                   'mtx-N-N.npy, where N = chromosome', required=True)
parser.add_argument('-i', '--interactions', type=str,
                    help="Name of the interactions provided while extracting", required=True)
parser.add_argument('-c', '--chr', type=str, help="Chromosome number", required=True)
parser.add_argument('-s', '--stats_folder', help="Folder to load the significant interactions from", type=str,
                    default='../stats/')
parser.add_argument('-f', '--figures_folder', help="Folder to save the plots in", type=str,
                    default='../figures/')
parser.add_argument('-t', '--threshold', type=float, help="Threshold that was used for statistical analysis")
parser.add_argument('-p', '--plot_title', type=str, help="The title of the plot", type=str,
                    default='Interactions')
parser.add_argument('-t', '--threshold', type=float, help="Threshold that was used for statistical analysis")


# Main
if __name__ == "__main__":
    args = parser.parse_args()
    p = Plotter(args.hic_folder, args.stats_folder, args.interactions, args.chr, args.threshold, args.plot_title)
    p.run(args.figures_folder)


# TODO add TAD plotting