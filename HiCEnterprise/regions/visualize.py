"""
Script for plotting interaction profiles for given regions from analysis.
"""

import argparse
import re
from glob import glob

import numpy as np

try:
    # Python3
    import pickle
except ImportError:
    # Python2
    import cPickle as pickle
import os
import sys
import logging
from textwrap import wrap
from ..utils import create_folders
from .extract import Bin

logger = logging.getLogger('regions.visualize')
logging.basicConfig(level=logging.DEBUG)  # TODO change


class Plotter:
    """
    Plots interaction profiles with either rpy2 or matplotlib.
    """

    def __init__(self, regions_name, chrom, num_regs, num_bins, threshold, bin_res):
        self.regions_name = regions_name
        self.chr = chrom
        self.num_regs = num_regs
        self.num_bins = num_bins
        self.threshold = threshold
        self.bin_res = bin_res

    def plot(self, regions, figures_folder, section, type):
        """
        Plots interaction profiles for chosen regions to pdf file either with rpy2 or matplotlib.
        """

        filename = figures_folder + '/' + self.regions_name + '-' + self.chr + '-' + str(self.num_bins)
        if section:
            filename += '-' + str(section[0]) + '-' + str(section[1])
        if self.num_regs:
            filename += '-' + str(self.num_regs)
        if os.path.exists(filename + '.pdf'):
            os.remove(filename + '.pdf')
            logger.debug('Removed the previous file with the same name: ' + filename + '.pdf')

        logger.debug('Plotting type provided: ' + type)

        if type == 'rpy2':
            self._plot_with_rpy2(regions, filename)
        elif type == 'mpl':
            self._plot_with_mpl(regions, filename)
        else:
            logger.critical('Error: unknown type for plotting. Available: mpl (matplotlib), rpy2.')
            sys.exit(1)

    def _plot_with_rpy2(self, regions, filename):
        from rpy2 import robjects
        import rpy2.robjects.lib.ggplot2 as ggplot2
        from rpy2.robjects.lib import grid
        from rpy2.robjects.packages import importr
        grdevices = importr('grDevices')
        base = importr('base')
        grdevices.pdf(file=filename + '.pdf')

        t = [x for x in range(-self.num_bins, self.num_bins + 1)]
        for region in regions[:self.num_regs]:
            if not np.any(region.weighted):
                logger.warning(
                    "Warning: No data for region located on bin " + str(region.bin) + ". Not plotting this one.")
                continue
            middle = (len(region.weighted[0]) - 1) / 2
            if middle < self.num_bins:
                logger.error("Warning: There are less bins calculated for regions than you want to plot.")
                sys.exit(1)
            d = {'map': robjects.StrVector(
                [str(m) for sublist in [[x] * len(t) for x in range(len(region.weighted))] for m in sublist]),
                't': robjects.FloatVector(t * len(region.weighted)),
                'e': robjects.FloatVector([i for sublist in region.weighted for i in
                                           sublist[middle - self.num_bins:middle + self.num_bins + 1]]),
                'p': robjects.FloatVector([-np.log10(x) for sublist in region.pvalues for x in
                                           sublist[middle - self.num_bins:middle + self.num_bins + 1]]),
                'c': robjects.FloatVector([-np.log10(x) for sublist in region.corrected_pvalues for x in
                                           sublist[middle - self.num_bins:middle + self.num_bins + 1]])}
            dataf = robjects.DataFrame(d)
            gp = ggplot2.ggplot(dataf)  # first yellow second red
            p1 = gp + ggplot2.geom_line(mapping=ggplot2.aes_string(x='t', y='e', group='map', colour='map'),
                                        alpha=0.8) + ggplot2.scale_y_continuous(trans='log2') + ggplot2.ggtitle(
                "\n".join(wrap("Bin " + str(region.bin) + " : " + str(region.positions)))) + ggplot2.labs(
                y="log Intensity") + ggplot2.theme_classic() + ggplot2.theme(
                **{'axis.title.x': ggplot2.element_blank(), 'axis.text.y': ggplot2.element_text(angle=45),
                   'axis.text.x': ggplot2.element_blank(),
                   'legend.position': 'none'}) + ggplot2.scale_colour_brewer(palette="Set1")
            p2 = gp + ggplot2.geom_line(mapping=ggplot2.aes_string(x='t', y='p', group='map', colour='map'),
                                        alpha=0.8) + ggplot2.labs(
                y="-log10(p-value)") + ggplot2.theme_classic() + ggplot2.theme(
                **{'axis.title.x': ggplot2.element_blank(), 'axis.text.x': ggplot2.element_blank(),
                   'legend.position': 'none'}) + ggplot2.scale_colour_brewer(palette="Set1")
            p3 = gp + ggplot2.geom_line(mapping=ggplot2.aes_string(x='t', y='c', group='map', colour='map'),
                                        alpha=0.8) + ggplot2.labs(y="-log10(q-value)",
                                                                  x='bins (' + str(self.bin_res) + ' bp each)') + \
                 ggplot2.geom_hline(mapping=ggplot2.aes_string(yintercept=str(-np.log10(self.threshold))),
                                    colour='black', alpha=0.8, linetype='dashed') + ggplot2.theme_classic() + \
                 ggplot2.theme(**{'legend.position': 'none'}) + ggplot2.scale_colour_brewer(palette="Set1")
            g1 = ggplot2.ggplot2.ggplotGrob(p1)
            g2 = ggplot2.ggplot2.ggplotGrob(p2)
            g3 = ggplot2.ggplot2.ggplotGrob(p3)
            robjects.globalenv["g"] = base.rbind(g1, g2, g3, size='first')
            robjects.r("grid::grid.draw(g)")
            grid.newpage()
            logger.debug('Plotted region ' + str(region.bin))

        grdevices.dev_off()

    def _plot_with_mpl(self, regions, filename):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(filename + '.pdf')

        t = [x for x in range(-self.num_bins, self.num_bins + 1)]
        n = len(regions[0].weighted)
        colors = cm.Set1(np.linspace(0, 1, n))

        for region in regions[:self.num_regs]:
            if not np.any(region.weighted):
                logger.warning(
                    "Warning: No data for region located on bin " + str(region.bin) + ". Not plotting this one.")
                continue
            middle = (len(region.weighted[0]) - 1) / 2
            if middle < self.num_bins:
                logger.error("Warning: There are less bins calculated for regions than you want to plot.")
                sys.exit(1)
            fig = plt.figure()
            ax = fig.add_subplot(311)
            ax.set_yscale("log")
            ax.set_ylabel('log Intensity', fontsize=8)
            ax2 = fig.add_subplot(312)
            ax2.set_ylabel('-log10(p-value)', fontsize=8)
            ax3 = fig.add_subplot(313)
            ax3.set_ylabel('-log10(q-value)', fontsize=8)
            ax3.set_xlabel('bins (' + str(self.bin_res) + ' bp each)', fontsize=8)
            title = ax.set_title("\n".join(wrap("Bin " + str(region.bin) + " : " + str(region.positions))))
            ax.axes.xaxis.set_ticklabels([])
            ax2.axes.xaxis.set_ticklabels([])
            for a in [ax, ax2, ax3]:
                a.tick_params(axis='both', which='major', labelsize=7)
                a.spines['top'].set_visible(False)
            for i in range(n):
                ax.plot(t, region.weighted[i][middle - self.num_bins:middle + self.num_bins + 1], c=colors[i],
                        alpha=0.8, linewidth=0.8)
                ax2.plot(t,
                         [-np.log10(x) for x in region.pvalues[i][middle - self.num_bins:middle + self.num_bins + 1]],
                         c=colors[i], alpha=0.8, linewidth=0.8)
                ax3.plot(t, [-np.log10(x) for x in
                             region.corrected_pvalues[i][middle - self.num_bins:middle + self.num_bins + 1]],
                         c=colors[i], alpha=0.8, linewidth=0.8)
            ax3.axhline(y=-np.log10(self.threshold), color='black', linestyle='dashed', alpha=0.8, linewidth=0.7)
            fig.tight_layout()
            title.set_y(1.05)
            plt.draw()
            pp.savefig()
            plt.close()
            logger.debug('Plotted region ' + str(region.bin))
        pp.close()

    def run(self, section, pickled_folder, figures_folder, type):
        """
        Runs the visualization according to provided args
        """
        pickled_folder, figures_folder = create_folders([pickled_folder, figures_folder])
        pck_filename_base = pickled_folder + '/' + self.chr + '-' + str(self.num_bins) + 'x' + str(
            self.bin_res) + 'bp-'
        reg_files = glob(
            "-[0-9]*x".join(
                pck_filename_base.split('-' + str(self.num_bins) + 'x')) + 'regs-' + self.regions_name)
        logger.debug('Reg files found: ' + str(reg_files))
        best = 0
        if reg_files:
            p = re.compile(r'-(?P<n_b>[0-9]*)x')
            best = max([int(p.search(x).group('n_b')) for x in reg_files])
            logger.info('Loading region data from pickled file')
            with open(("-" + str(best) + "x").join(
                    pck_filename_base.split('-' + str(self.num_bins) + 'x')) + 'regs-' + self.regions_name,
                      'rb') as fp:
                regs = pickle.load(fp)
        else:
            logger.critical('No pickled region file '+pck_filename_base + 'regs' + self.regions_name+' in ' + pickled_folder)
            sys.exit(1)

        if section:
            section = section.split(':')
            section = tuple(int(x) for x in section)
            regs = [r for r in regs if
                    r.positions[0][0] <= section[1] and r.positions[-1][1] >= section[0]]
        if not self.num_regs:
            self.num_regs = len(regs)
        logger.info('Plotting for ' + str(self.num_regs) + ' first unique-bin regions, on span of ' + str(
            self.num_bins) + ' bins...')

        self.plot(regs, figures_folder, section, type)

        logger.info('Plotting finished.')


# Argument parsing
parser = argparse.ArgumentParser(
    description='Script for visualizing hicMaps with significant points (& tss enrichments) from extract.py '
                '(and enrichments.py) scripts')
parser.add_argument('-r', '--regions_name',
                    help='Name of file from which extract.py generated profiles, ie Fetal_brain', type=str)
parser.add_argument('-c', '--chr', help='Chromosome for which to extract bins', type=str,
                    default='1')  # it has to be str, because of 'X'
parser.add_argument('-b', '--bin_res', help='Resolution (size of the bins on the hicmaps) in bp', type=int)
parser.add_argument('-n', '--num_bins', help='Number of bins left and right to be plotted', type=int, default=100)
parser.add_argument('--num_regs', help='Number of regions to consider and plot', type=int)
parser.add_argument('-t', '--threshold', help='FDR threshold', type=float,
                    default=0.01)
parser.add_argument('-s', '--section', help="Section of bp for the plot", type=str)
parser.add_argument('-p', '--pickled_folder', help="Folder with pickled files (to load from)", type=str,
                    default='../pickles/')
parser.add_argument('-f', '--figures_folder', help="Folder to save the figures (plots) in", type=str,
                    default='../plots/')
parser.add_argument('--type', help="If the plotting should be with rpy2 or matplotlib. Options: mpl, rpy2. Default is"
                                   " mpl.", type=str, default='mpl')

# Main

if __name__ == "__main__":
    args = parser.parse_args()
    p = Plotter(args.regions_name, args.chr, args.num_regs, args.num_bins, args.threshold, args.bin_res)
    p.run(args.section, args.pickled_folder, args.figures_folder, args.type)
