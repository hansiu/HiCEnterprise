"""
Script for extracting significant domain-domain contact frequencies for given domains and HiC maps.
Idea by Irina Tuszynska. Based on code by Irina Tuszynska.
"""

import os
import numpy as np
import scipy.stats as ss
import csv
import logging
import sys
from operator import itemgetter
from HiCEnterprise.utils import create_folders, load_hicmap

from statsmodels.sandbox.stats.multicomp import fdrcorrection0

logger = logging.getLogger('domains.extract')
logging.basicConfig(level=logging.INFO)


class Extractor:
    """
    Performs all extracting of significant contact frequencies from provided maps, domains and arguments.
    """

    def __init__(self, domain_file, chrom, bin_res, sherpa_lvl, hic_folder, threshold):
        self.domains_name = str(os.path.basename(domain_file).split('.')[0])
        self.chr = chrom
        self.bin_res = bin_res
        self.sherpa_level = sherpa_lvl or None
        self.domains = self._load_domains(domain_file)
        self._check_overlap()
        self.hic_folder = os.path.abspath(hic_folder)
        self.hicmap = load_hicmap(hic_folder, 'mtx-' + self.chr + '-' + self.chr + '.npy')
        self.hic_name = os.path.basename(os.path.abspath(hic_folder))
        self.threshold = threshold

    def _load_domains(self, domain_file):
        logger.info('Loading Domains: ' + domain_file)
        domain_file = open(domain_file, 'r')
        reader = csv.reader(domain_file, delimiter=' ')
        domains = []
        for row in reader:
            if row[1] == self.chr:
                if (not self.sherpa_level) or (len(row) >= 5 and row[4] == str(self.sherpa_level)):
                    domains.append(
                        [int(int(row[2]) / self.bin_res), int(int(row[3]) / self.bin_res) - 1])
        domains = sorted(domains, key=itemgetter(0))
        return domains

    def _check_overlap(self):
        for i in range(len(self.domains) - 1):
            if self.domains[i + 1][0] <= self.domains[i][1]:
                logger.error('The domains ' + str(self.domains[i]) + str(',') + str(
                    self.domains[i + 1]) + ' are overlapping. Check input data')
                sys.exit(1)

    def create_domain_matrix(self):
        """
        Creates the domain matrix. Each cell contains interactions between two domains (by indices in self.domains)

        """
        logger.info('Creating the domain matrix')
        domain_matrix = np.zeros((len(self.domains), len(self.domains)))
        for i1, d1 in enumerate(self.domains):
            for i2, d2 in enumerate(self.domains):
                domain_matrix[i1, i2] = self.hicmap[d1[0]:d1[1] + 1, d2[0]:d2[1] + 1].sum()

        return domain_matrix

    def calc(self, domain_matrix):
        """
        Calculates statistical signigicance of inter-domain contacts using hypergeometric test.
        :param domain_matrix: Matrix with frequencies of contacts between domains
        """
        logger.info('Calculating the probability of inter-domain contacts (this may take some time...)')
        sigs = []
        n = np.sum(domain_matrix)
        pvalue_matrix = np.zeros(domain_matrix.shape)
        for i in range(domain_matrix.shape[0]):
            a = sum(domain_matrix[i][:])
            for j in range(domain_matrix.shape[1]):
                b = sum(domain_matrix[j][:])
                model = ss.hypergeom(n, a, b)
                expected = (a*b)/n
                k = domain_matrix[i][j]
                if expected and k:
                    pval = model.sf(k)
                else:
                    pval = 1.0
                pvalue_matrix[i, j] = pval

                if pval < self.threshold and pval != 0:
                    sigs.append((i, j, pval))

        logger.info('FDR correcting the pvalues')

        corrected = fdrcorrection0(pvalue_matrix.flatten(), self.threshold)[1].reshape(domain_matrix.shape)

        corr_sigs = []
        for i in range(domain_matrix.shape[0]):
            for j in range(domain_matrix.shape[1]):
                pval = corrected[i,j]
                if pval < self.threshold and pval != 0:
                    corr_sigs.append((i, j, pval))

        return sigs,corr_sigs

    def save_sigs(self, sigs, stats_folder, corr=""):
        """
        Saves significant interactions between domains to givenn folder

        """
        stats_folder = create_folders([stats_folder])[0]
        logger.debug(
            'Saving stats file: ' + stats_folder + '/' + self.domains_name + '-' + self.hic_name + '-' + corr + 'stats' + self.chr + '-' + "_".join(
                str(
                    self.threshold).split('.')) + '.txt')
        stats = open(
            stats_folder + '/' + self.domains_name + '-' + self.hic_name + '-' + corr + 'stats' + self.chr + '-' + "_".join(
                str(
                self.threshold).split('.')) + '.txt',
                     'w')
        stats_writer = csv.writer(stats, delimiter='\t')
        rows = []
        for sig in sigs:
            rows.append([self.chr] + self.domains[sig[0]] + self.domains[sig[1]] + [sig[2]])
        stats_writer.writerows(rows)

    def run(self, stats_folder, plotting, figures_folder):
        """
        Runs the analysis
        """
        dom_matrix = self.create_domain_matrix()
        sigs, corr_sigs = self.calc(dom_matrix)
        self.save_sigs(sigs, stats_folder)
        self.save_sigs(corr_sigs, stats_folder, corr="corr_")

        if plotting is not False:
            logger.debug('Getting to plotter')
            from .visualize import Plotter
            p = Plotter(self.hic_folder, stats_folder, self.domains_name, self.chr, self.threshold)
            p.run(figures_folder)
