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

    def __init__(self, domain_file, chrom, bin_res, sherpa_lvl, hic_map, threshold, plot_title, ticks_separ,
                 hic_color, inter_color, distribution, hic_name, all_dom, inter_indom):
        self.domains_name = str(os.path.basename(domain_file).split('.')[0])
        self.chr = chrom
        self.bin_res = bin_res
        self.sherpa_level = sherpa_lvl or None
        self.domains = self._load_domains(domain_file)
        self._check_overlap()
        hic_folder = os.path.dirname(os.path.abspath(hic_map))
        filename = os.path.basename(os.path.abspath(hic_map))
        self.hicmap = self._symm(load_hicmap(hic_folder, filename))
        self.hic_name = hic_name or os.path.basename(hic_map).split('.')[0]
        self.hic_path = hic_map
        self.threshold = threshold
        self.plot_title = plot_title
        self.ticks_separ = ticks_separ
        self.hic_color = hic_color
        self.inter_color = inter_color
        self.distribution = distribution
        self.all_dom = all_dom
        self.inter_indom = inter_indom

    def _symm(self, hicmap):
        hicmap = np.nan_to_num(hicmap)
        return ((hicmap + np.transpose(hicmap)) / 2)


    def _load_domains(self, domain_file):
        logger.info('Loading Domains: ' + domain_file)
        domain_file = open(domain_file, 'r')
        reader = csv.reader(domain_file, delimiter=' ')
        domains = []
        warned = False
        for row in reader:
            if len(row) > 5 and not warned:
                logger.warning(
                    'There are more than 5 columns in the input domains. This may suggest that your input domains '
                    'are in the wrong format.')
                warned = True
            if len(row) < 4:
                logger.critical('There is something wrong with the input domains: not enough columns.')
                sys.exit(1)

            try:
                row = [int(row[0]), str(row[1]).lstrip('chr')] + list(map(int, row[2:5]))
            except ValueError:
                logger.critical(
                    'There is something wrong with the input domains: columns are in the wrong format (should be '
                    'integers).')
                sys.exit(1)
            if row[1] == self.chr:
                if (not self.sherpa_level) or (len(row) >= 5 and row[4] == self.sherpa_level):
                    domains.append(
                        [int(row[2] / self.bin_res), int(row[3] / self.bin_res) - 1])
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
        dom_size = np.zeros(len(self.domains)) #moja
        for i1, d1 in enumerate(self.domains):
            for i2, d2 in enumerate(self.domains):
                domain_matrix[i1, i2] = self.hicmap[d1[0]:d1[1] + 1, d2[0]:d2[1] + 1].sum()
            dom_size[i1] =  d1[1]+1-d1[0] #moja
        if self.all_dom == False:
            dom_sum = domain_matrix.sum(axis =0) # in order to eliminate near to centromere domains with small amount of contacts we calculate the mean number of contacts in one row inside each domain, and if it is less than  n*len_domains (default n = 1), all values for this domain = 0.0
            points_in_one_row_domains = np.divide(dom_sum, dom_size)
            suspicious_dom = np.nonzero(points_in_one_row_domains < self.inter_indom * len(self.domains))
            if suspicious_dom[0].shape != (0,):
                for d in suspicious_dom[0]:
                    domain_matrix[d] = domain_matrix[:,d] = 0.0
            else: pass
        else: pass
        if domain_matrix[np.nonzero(domain_matrix)].size == 0:
            logger.error("WARNING! Due to the sparsity of the matrix, please use --all_domains option")
            sys.exit(1)
        return domain_matrix

    def calc(self, domain_matrix):
        """
        Calculates statistical signigicance of inter-domain contacts using hypergeometric test.
        :param domain_matrix: Matrix with frequencies of contacts between domains
        """
        logger.info('Calculating the probability of inter-domain contacts with ' + self.distribution + ' distribution '
                                                                                                       '(this may take some time...)')

        if self.distribution == 'hypergeom':
            logger.info("Hypergeometric test is calculating")
            sigs, corr_sigs = self._calc_hypergeom(domain_matrix)
        elif self.distribution == 'negbinom':
            logger.info("Negativ binomial test is calculating")
            sigs, corr_sigs = self._calc_negbinom(domain_matrix)
        elif self.distribution == 'poisson':
            logger.info("Poisson test is calculating")
            sigs, corr_sigs = self._calc_poisson(domain_matrix)
        else:
            logger.error('The chosen distribution ' + self.distribution + ' is not available. Choose from: hypergeom, '
                                                                          'negbinom, poisson')
            sys.exit(1)

        return sigs, corr_sigs

    def _calc_hypergeom(self, domain_matrix):
        sigs = []
        n = np.sum(np.triu(domain_matrix))
        #n =  np.sum(domain_matrix)
        pvalue_matrix = np.ones(domain_matrix.shape)
        for i in range(domain_matrix.shape[0]):
            a = sum(domain_matrix[i][:])
            for j in range(domain_matrix.shape[1]):
                b = sum(domain_matrix[j][:])
                model = ss.hypergeom(n, a, b)
                expected = (a * b) / n
                k = domain_matrix[i][j]
                #print expected, k, a, b, n
                if expected and k:
                    pval = model.sf(k)
                else:
                    pval = 1.0
                pvalue_matrix[i, j] = pval

                if pval < self.threshold:
                    sigs.append((i, j, pval))

        return sigs, self._fdr_correct(pvalue_matrix, domain_matrix.shape)

    def _calc_negbinom(self, domain_matrix):
        sigs = []
        means = [np.mean(self.hicmap.diagonal(i)) for i in range(self.hicmap.shape[0])]
        lens = [len(self.hicmap.diagonal(i)) for i in range(self.hicmap.shape[0])]

        def sum_mean(i, j):
            """
            Counts the mean across several consecutive diagonals in the hicmap
            """
            s = sum([m * l for (m, l) in list(zip(means, lens))[i:j]])
            l = sum(lens[i:j])
            try:
                return s / l
            except ZeroDivisionError:
                return 0

        def sum_var(i, j, m):
            """
            Counts the variance in several consecutive diagonals given their mean
            """
            mses = [np.mean((self.hicmap.diagonal(i) - m) ** 2) for i in range(i, j)]
            s = sum([m * l for (m, l) in zip(mses, lens[i:j])])
            l = sum(lens[i:j])
            try:
                return s / l
            except:
                return 0

        pvalue_matrix = np.ones(domain_matrix.shape)

        for i in range(domain_matrix.shape[0]):
            li = self.domains[i][1] - self.domains[i][0] + 1
            for j in range(i + 1, domain_matrix.shape[1]):
                lj = self.domains[j][1] - self.domains[j][0] + 1
                dist = self.domains[j][0] - self.domains[i][1]
                span = self.domains[j][1] - self.domains[i][0] + 1
                expected = sum_mean(dist, span)
                var = sum_var(dist, span, expected)
                mean = expected * li * lj
                if var < mean:
                    var = mean + 1
                k = domain_matrix[i][j]
                r = mean ** 2 / (var - mean) if (var - mean) != 0 else np.nan
                p = (var - mean) / var if var != 0 else np.nan
                model = ss.nbinom(n=r, p=1 - p)
                if expected and k:
                    pval = model.sf(k)
                    pvalue_matrix[i, j] = pval
                    if pval < self.threshold:
                        sigs.append((i, j, pval))

        return sigs, self._fdr_correct(pvalue_matrix, domain_matrix.shape)

    def _calc_poisson(self, domain_matrix):
        sigs = []
        means = [np.mean(self.hicmap.diagonal(i)) for i in range(self.hicmap.shape[0])]
        lens = [len(self.hicmap.diagonal(i)) for i in range(self.hicmap.shape[0])]

        def sum_mean(i, j):
            """
            Counts the mean across several consecutive diagonals in the hicmap
            """
            s = sum([m * l for (m, l) in list(zip(means, lens))[i:j]])
            l = sum(lens[i:j])
            try:
                return s / l
            except:
                return 0

        pvalue_matrix = np.ones(domain_matrix.shape)

        for i in range(domain_matrix.shape[0]):
            li = self.domains[i][1] - self.domains[i][0] + 1
            for j in range(i + 1, domain_matrix.shape[1]):
                lj = self.domains[j][1] - self.domains[j][0] + 1
                dist = self.domains[j][0] - self.domains[i][1]
                span = self.domains[j][1] - self.domains[i][0] + 1
                expected = sum_mean(dist, span) * li * lj
                k = domain_matrix[i][j]
                model = ss.poisson(expected)
                if expected and k:
                    pval = model.sf(k)
                    pvalue_matrix[i, j] = pval
                    if pval < self.threshold:
                        sigs.append((i, j, pval))

        return sigs, self._fdr_correct(pvalue_matrix, domain_matrix.shape)

    def _fdr_correct(self, pvalue_matrix, dom_shape):
        logger.info('FDR correcting the pvalues')
        pvalue_matrix[np.isnan(pvalue_matrix)] = 1.0
        corrected = fdrcorrection0(pvalue_matrix.flatten(), self.threshold)[1].reshape(dom_shape)
        np.save(open('pval.npy','w'), pvalue_matrix)
        #for i,j in  zip(pvalue_matrix.flatten(),corrected.flatten()):
        #    if i<1.0:
        #        print 'p-qval', i,j

        corr_sigs = []
        for i in range(dom_shape[0]):
            for j in range(dom_shape[1]):
                pval = corrected[i, j]
                if pval < self.threshold:
                    corr_sigs.append((i, j, pval))

        return corr_sigs

    def save_sigs(self, sigs, stats_folder, corr=""):
        """
        Saves significant interactions between domains to givenn folder

        """
        stats_folder = create_folders([stats_folder])[0]
        logger.debug(
            'Saving stats file: ' + stats_folder + '/' + self.domains_name + '-' + self.hic_name + '-' + corr +
            'stats' + self.chr + '-' + "_".join(str(self.threshold).split('.')) + '-' + self.distribution + '.txt')
        stats = open(
            stats_folder + '/' + self.domains_name + '-' + self.hic_name + '-' + corr + 'stats' + self.chr + '-' +
            "_".join(str(self.threshold).split('.')) + '-' + self.distribution + '.txt', 'w')
        stats_writer = csv.writer(stats, delimiter='\t')
        stats_writer.writerow(['chr', 'd1_start', 'd1_end', 'd2_start', 'd2_end', corr + 'pval'])
        rows = []
        for sig in sigs:
            rows.append([self.chr] + self._get_coordinates(self.domains[sig[0]]) + self._get_coordinates(
                self.domains[sig[1]]) + [sig[2]])
        stats_writer.writerows(rows)

    def _get_coordinates(self, domain):
        return [domain[0] * self.bin_res, (domain[1] + 1) * self.bin_res]

    def run(self, stats_folder, plotting, figures_folder):
        """
        Runs the analysis
        """
        dom_matrix = self.create_domain_matrix()
        sigs, corr_sigs = self.calc(dom_matrix)
        print "sigs", sigs
        print 'corr', corr_sigs
        if sigs == []:
            logger.error('There is no p-val lower tan threeshold. Change threeshold value to bigger value and run again. Try 0.01, 0.05 or 0.1')
            sys.exit(1)
        else:
            self.save_sigs(sigs, stats_folder)
        if corr_sigs == []:
            logger.info('There is no appropriate q-value after FDR correction!')
        else:
            self.save_sigs(corr_sigs, stats_folder, corr="corr_")

        if plotting is not False:
            logger.debug('Getting to plotter')
            from .visualize import Plotter
            p = Plotter(self.hic_path, stats_folder, self.domains_name, self.chr, self.threshold, self.plot_title,
                        self.ticks_separ, self.hic_color, self.inter_color, self.bin_res, self.distribution,
                        self.hic_name, self.all_dom, self.inter_indom)
            p.run(figures_folder)
