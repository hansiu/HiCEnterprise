"""
Script for extracting interaction profiles and significant long-distance contact frequencies for given regions and HiC
maps.
"""

import csv
import logging
import numpy as np
import os
import random
import re
import scipy.stats as ss
import sys
from glob import glob
from ..utils import create_folders, load_hicmap

from statsmodels.sandbox.stats.multicomp import fdrcorrection0

try:
    # Python3
    import pickle
except ImportError:
    # Python2
    import cPickle as pickle

logger = logging.getLogger('regions.extract')
logging.basicConfig(level=logging.INFO)


class Bin:
    """
    Container for regions occurring in one bin.
    """

    def __init__(self, bin, positions):
        self.positions = positions
        self.bin = bin
        self.intensities = []
        self.weighted = []
        self.pvalues = []
        self.corrected_pvalues = []


class HiCmap:
    """
     class for HiC map operations
    """

    def __init__(self, hic_folder, chr):
        self.folder = os.path.abspath(hic_folder)
        self.filename = '/mtx-' + chr + '-' + chr + '.npy'
        self.loaded = False
        self.map = None
        self.fits = {}
        self.means = {}

    def get_name(self):
        """
        Returns the name of the highest folder that map is kept in.
        :return:
        """
        return "_".join(os.path.basename(self.folder).split())

    def load(self):
        """
        Loads HiC Map for given chromosome from the given format and checks the symmetry.
        Maps names should be in 'mtx-N-N.npy' form, where N is the chromosome name.
        Files should contain one HiC map in numpy format.
        """
        if not self.loaded:
            self.map = load_hicmap(self.folder, self.filename)
            self.loaded = True

    def get_means(self, num_bins):
        """
        Calculates means from lines parallel to the diagonal that will be used as weights for intensities.
        This is necessary to account for the lower number of contacts on bigger distances.

        """
        self.load()

        logger.info('Calculating the missing mean intensities for diagonals for Map ' + self.get_name())

        for i in range(-num_bins, num_bins + 1):
            if i not in self.means.keys():
                self.means[i] = sum(self.map.diagonal(i)) / len(self.map.diagonal(i)) if self.map.diagonal(i).size and \
                                                                                         sum(self.map.diagonal(
                                                                                             i)) else 1.0
                logger.debug('Calculated missing mean for ' + str(i) + ': ' + str(self.means[i]))

    def get_Weibull_fits(self, num_bins):
        """
        Calculates Weibull fits parameters from lines parallel to the diagonal that will be used to calculate the
        significance of contact frequencies. Idea by (Won, 2016)

        """
        self.load()

        logger.info('Getting missing Weibull fits for Map ' + self.get_name())

        for i in range(-num_bins, num_bins + 1):
            if i not in self.fits.keys():
                diag = self.map.diagonal(i)
                try:
                    p95 = np.percentile(diag[diag.nonzero()], 95)
                except IndexError:
                    p95 = 0.0
                diag = sorted([x for x in diag[diag.nonzero()] if x < p95])
                if diag:
                    self.fits[i] = ss.exponweib.fit(diag, floc=0)
                else:
                    self.fits[i] = None
                logger.debug('Calculated missing fit for ' + str(i) + ': ' + str(self.fits[i]))


class Extractor:
    """
    Performs all the extracting of interaction profiles and significant contact frequencies from provided
    maps, regions and arguments.
    """

    def __init__(self, regions_name, hic_folders, chrom, bin_res, num_bins, threshold, pickled_folder, figures_folder,
                 stats_folder):
        self.regions_name = regions_name
        self.maps = []
        for folder in sorted(hic_folders):
            self.maps.append(HiCmap(folder, chrom))
            logger.debug('Created a map for ' + folder + ', ' + chrom)
        self.chr = chrom
        self.bin_res = bin_res
        self.num_bins = num_bins
        self.threshold = threshold
        self.pickled_folder, self.figures_folder, self.stats_folder = create_folders([pickled_folder, figures_folder,
                                                                                      stats_folder])

    def get_regions(self, reg_file, algorithm='all'):
        """
        Reads regions from given file. Supported file formats: Tab-delimited BED or FASTA from EnhancerAtlas with
        positions in the header of the sequence.

        :param reg_file: path to the file with regions
        :param algorithm: algorithm type:
            - 'all' is set as a default - region is set to every bin it is even partially in
            - 'one' is optional - region is set to the bin in which its biggest part belongs OR in case there are few
               of equal length to random choice from them
        :return: list of loaded regions (list of Bins)
        """

        format = reg_file.split('.')[-1]
        if format.lower() == 'fasta':
            regs = self._get_reg_fasta(reg_file, algorithm)
        elif format.lower() == 'bed':
            regs = self._get_reg_bed(reg_file, algorithm)
        else:
            logger.critical('Unsupported file format.')
            sys.exit(1)
        return regs

    def _get_reg_fasta(self, reg_file, algorithm):
        logger.debug('Loading regions from FASTA file with alg ' + algorithm)
        regs = []
        regseq = open(reg_file, 'r')
        bins = {}
        for line in regseq:
            if line.startswith(">chr" + self.chr + ":"):
                position = [int(x) for x in line.split(":")[1].split("_")[0].split('-')]
                if algorithm == 'all':
                    bin = tuple(int(x / self.bin_res) for x in position)
                    for b in range(bin[0], bin[1] + 1):
                        if b in bins.keys():
                            bins[b].append(position)
                        else:
                            bins[b] = [position]
                elif algorithm == 'one':
                    b = self._get_reg_justone(position)
                    if b in bins.keys():
                        bins[b].append(position)
                    else:
                        bins[b] = [position]
        for b in sorted(bins.keys()):
            regs.append(Bin(b, bins[b]))
        logger.debug('Found ' + str(len(bins)) + ' bins and ' + str(len(regs)) + ' regions')
        return regs

    def _get_reg_bed(self, reg_file, algorithm):
        logger.debug('Loading regions from BED file with alg ' + algorithm)
        regs = []
        regseq = open(reg_file, 'r')
        bins = {}
        for line in regseq:
            if line.startswith("chr" + self.chr + '\t') or line.startswith(self.chr + '\t'):
                position = [int(x) for x in line.split()[1:3]]
                position[1] -= 1
                if algorithm == 'all':
                    bin = tuple(int(x / self.bin_res) for x in position)
                    for b in range(bin[0], bin[1] + 1):
                        if b in bins.keys():
                            bins[b].append(position)
                        else:
                            bins[b] = [position]
                elif algorithm == 'one':
                    b = self._get_reg_justone(position)
                    if b in bins.keys():
                        bins[b].append(position)
                    else:
                        bins[b] = [position]
        for b in sorted(bins.keys()):
            regs.append(Bin(b, bins[b]))
        logger.debug('Found ' + str(len(bins)) + ' bins and ' + str(len(regs)) + ' regions')
        return regs

    def _get_reg_justone(self, position):
        bin = tuple(int(x / self.bin_res) for x in position)
        rem = tuple(x % self.bin_res for x in position)
        if bin[1] - bin[0] == 1:
            if self.bin_res - rem[0] > rem[1]:
                b = bin[0]
            elif self.bin_res - rem[0] < rem[1]:
                b = bin[1]
            else:
                b = random.choice(bin)
        elif bin[1] - bin[0] >= 2:
            b = random.choice(range(bin[0] + 1, bin[1]))
        else:
            b = bin[0]
        return b

    def set_regs_intensities(self, regs, hicmap):
        """
        Sets intensity profiles for regions (Bins with regions) from a given HiC map.

        :param regs: list ot Bins (with regions)
        :param hicmap: HiC map loaded as numpy array
        """

        hicmap.load()
        np.nan_to_num(hicmap)

        for region in regs:
            region.intensities.append([])
            for i in range(-self.num_bins, self.num_bins + 1):
                if i < 0:
                    h = region.bin + i
                else:
                    h = region.bin
                if hicmap.map.diagonal(i).size <= abs(h) or h < 0:
                    region.intensities[-1].append(0.0)
                else:
                    region.intensities[-1].append(hicmap.map.diagonal(i)[h])
            logger.debug('Added intensities to region ' + str(region.bin))

    def calc(self, regs, hicmap):
        """
        Calculates weighted intensities, p-values and q-values (FDR-corrected p-values) for region intensity profiles.

        :param hicmap: HiC map to run calculations for
        :param regs: list ot Bins (with regions)
        """

        pvals = []
        ints = []
        for region in regs:
            pvals.append([])
            for i in range(-self.num_bins, self.num_bins + 1):
                if hicmap.fits[i] is not None:
                    pvals[-1].append(ss.exponweib.sf(region.intensities[len(region.pvalues)][len(pvals[-1])],
                                               *hicmap.fits[i]))
                else:
                    pvals[-1].append(1.0)
            ints.append(region.intensities[len(region.pvalues)])
            p = np.array(pvals[-1])
            pvals[-1] = p
            p[p == 0.0] = 0.000000000000000000000000000001  # TODO delete that? useful for log representation
            region.pvalues.append(p)
        logger.debug('Calculated pvalues for map ' + hicmap.get_name())
        pvals, ints, corr_big = np.array(pvals), np.array(ints), np.array(ints)

        corrected = fdrcorrection0(np.array(pvals)[ints.nonzero()], self.threshold)[1]
        logger.debug('Calculated qvalues for map ' + hicmap.get_name())
        corrected[corrected == 0.0] = 0.000000000000000000000000000001 # TODO delete that? useful for log representation
        corr_big[corr_big.nonzero()] = corrected
        corr_big[np.nonzero(corr_big == 0.0)] = 1.0
        corr_big.reshape(pvals.shape)

        for r in range(len(regs)):
            x = len(regs[r].corrected_pvalues)
            regs[r].corrected_pvalues.append(corr_big[r])
            regs[r].weighted.append(
                [0.0 if np.isnan(hicmap.means[-(int(len(regs[r].intensities[x]) / 2.0) - regs[r].intensities[x].index(y))]) else y / hicmap.means[-(int(len(regs[r].intensities[x]) / 2.0) - regs[r].intensities[x].index(y))]  for y in
                 regs[r].intensities[x]])
            

    def save_stats_and_sigs(self, regions, single_sig, remap, formats, h_fs):
        """
        Saves statistics from the analysis and significant points in the custom txt, tab-delimited bed or gff format

        :param h_fs:
        :param regions: list ot Bins (with regions)
        :param single_sig: if statistics for the single maps should be saved or not
        :param remap: if remap and from which assembly to which remap the positions
        :param formats: list of formats, available: txt, bed, gff
        """
        formats = list({'bed', 'txt', 'gff'} & {f.lower() for f in formats})

        if remap and formats != ['gff']:
            logger.warning('For now remapping between assemblies is available only for GFF output format')

        if not formats:
            logger.warning('No accepted formats for significant contact file provided, no files will be generated ('
                        'available: txt, bed, gff)')

        rest_of_name = 'x' + str(self.bin_res) + 'bp-' + "_".join(str(self.threshold).split('.'))
        significant = {}
        single_sigs = {}
        t = [x for x in range(-self.num_bins, self.num_bins + 1)]
        r_pred = 0
        pred = 0

        if remap is not None and 'gff' in formats:
            logger.debug('Setting up pyliftover')
            from pyliftover import LiftOver
            lo = LiftOver(*remap.split(':'))
            rest_of_name += '_remapped' + remap.split(':')[1]
        else:
            lo = None

        for f in formats:
            sign_file = open(
                self.stats_folder + '/' + self.regions_name + '-' + "-".join(
                    h_fs) + '-significant' + self.chr + '-' + str(
                    self.num_bins) +
                rest_of_name + '.' + f, 'w')
            if f != 'txt':
                sign_file.writelines(['#track name="LongRanger predictions ' + self.regions_name + ' chr' +
                                      self.chr + ' FDR' + str(self.threshold) + '"\n',"#chr start    end chr_int  start_int end_int -log10(q-val)\n" ])
                #sign_file.writelines(['track name="LongRanger predictions ' + self.regions_name + ' chr' +
                #                      self.chr + ' FDR' + str(self.threshold) + '"\n'])
            significant[f] = csv.writer(sign_file, delimiter='\t')
            logger.debug('Created significant file writer for ' + f + ' format')
            if single_sig and len(regions[0].weighted) > 1:
                for i in range(len(regions[0].weighted)):
                    single_file = open(self.stats_folder + '/' + self.regions_name + '-' + h_fs[
                        i] + 'single_sig' + self.chr + '-' + str(self.num_bins) + rest_of_name + '.' + f, 'w')
                    if f != 'txt':
                        single_file.writelines(['track name="LongRanger predictions ' + self.regions_name + ' chr' +
                                                self.chr + ' FDR' + str(self.threshold) + ' map' + str(i) + '"\n'])
                    single_file.close()
                single_sigs[f] = [
                    csv.writer(open(self.stats_folder + '/' + self.regions_name + '-' + h_fs[i] + 'single_sig' +
                                    self.chr + '-' + str(self.num_bins) + rest_of_name + '.' + f, 'a'), delimiter='\t')
                    for i in range(len(regions[0].weighted))]
                logger.debug('Created significant single file writers for ' + f + ' format')

        for reg in regions:
            sig = []
            supp = set(
                [j for j, v in enumerate(reg.corrected_pvalues[0]) if v < self.threshold])

            if supp and single_sig and len(reg.weighted) > 1:
                s = []
                supp = sorted(list(supp))
                inds = [x + 1 for x, v in enumerate(supp) if
                        x == len(supp) - 1 or supp[x] + 1 != supp[x + 1]][:-1]
                pairs = [tuple([p[0], p[-1]]) for p in
                         np.split(np.array(supp), inds)]
                for pair in pairs:
                    s.append([tuple([reg.bin + t[pair[0]], reg.bin + t[pair[-1]]]),
                              reg.corrected_pvalues[0][pair[0]:pair[1] + 1]])
                if s:
                    if 'txt' in formats:
                        single_sigs['txt'][0].writerow([str(reg.bin) + ';' + str(reg.positions),
                                                        ";".join([",".join([str(el) for el in line]) for line in s])])
                    if 'bed' in formats:
                        for x in s:
                            for pos in reg.positions:
                                for b in range((x[0][1] - x[0][0]) + 1):
                                    single_sigs['bed'][0].writerow(
                                        ["chr" + self.chr, str(pos[0]), str(pos[1] + 1), "chr" + self.chr,  (x[0][0] + b) * self.bin_res,
                                         (x[0][0] + b + 1) * self.bin_res, -np.log10(x[1][b])])
                                    #single_sigs['bed'][0].writerow(
                                    #    ["chr" + self.chr, (x[0][0] + b) * self.bin_res,
                                    #     (x[0][0] + b + 1) * self.bin_res, "chr" + self.chr + "_" + str(pos[0]) +
                                    #     "_" + str(pos[1] + 1), -np.log10(x[1][b])])
                    if 'gff' in formats:
                        if not self._save_sig_gff(single_sigs['gff'][0], s, True, reg.positions, reg.bin, remap, lo):
                            continue

            supp_in_all_ind = set(supp)

            for i in range(1, len(reg.corrected_pvalues)):
                supp = set([j for j, v in enumerate(reg.corrected_pvalues[i]) if v < self.threshold])
                if supp and single_sig:
                    s = []
                    supp = sorted(list(supp))
                    inds = [x + 1 for x, v in enumerate(supp) if
                            x == len(supp) - 1 or supp[x] + 1 != supp[x + 1]][:-1]
                    pairs = [tuple([p[0], p[-1]]) for p in
                             np.split(np.array(supp), inds)]
                    for pair in pairs:
                        s.append([tuple([reg.bin + t[pair[0]], reg.bin + t[pair[-1]]]),
                                  reg.corrected_pvalues[i][pair[0]:pair[1] + 1]])
                    if s:
                        if 'txt' in formats:
                            single_sigs['txt'][i].writerow(
                                [str(reg.bin) + ';' + str(reg.positions),
                                 ";".join([",".join([str(el) for el in line]) for line in s])])
                        if 'bed' in formats:
                            for x in s:
                                for pos in reg.positions:
                                    for b in range((x[0][1] - x[0][0]) + 1):
                                        single_sigs['bed'][i].writerow(
                                            ["chr" + self.chr, str(pos[0]), str(pos[1] + 1), "chr" + self.chr,  (x[0][0] + b) * self.bin_res,
                                             (x[0][0] + b + 1) * self.bin_res, -np.log10(x[1][b])])
                                        #single_sigs['bed'][i].writerow(
                                        #    ["chr" + self.chr, (x[0][0] + b) * self.bin_res,
                                        #     (x[0][0] + b + 1) * self.bin_res, "chr" + self.chr + "_" +
                                        #     str(pos[0]) + "_" + str(pos[1] + 1), -np.log10(x[1][b])])
                        if 'gff' in formats:
                            if not self._save_sig_gff(single_sigs['gff'][i], s, True, reg.positions, reg.bin, remap,
                                                      lo):
                                continue

                supp_in_all_ind &= set(supp)
            supp_in_all_ind = sorted(list(supp_in_all_ind))

            if supp_in_all_ind:
                inds = [x + 1 for x, v in enumerate(supp_in_all_ind) if
                        x == len(supp_in_all_ind) - 1 or supp_in_all_ind[x] + 1 != supp_in_all_ind[x + 1]][:-1]
                pairs = [tuple([p[0], p[-1]]) for p in np.split(np.array(supp_in_all_ind), inds)]
                for pair in pairs:
                    sig.append([tuple([reg.bin + t[pair[0]], reg.bin + t[pair[-1]]]), list(zip(
                        *[reg.corrected_pvalues[i][pair[0]:pair[1] + 1] for i in
                          range(len(reg.corrected_pvalues))]))])

            if sig:
                r_pred += 1
                pred += len(sig)
                if 'txt' in formats:
                    significant['txt'].writerow([str(reg.bin) + ';' + str(reg.positions), ";".join(
                        [",".join([str(el) for el in line]) for line in sig])])
                if 'bed' in formats:
                    for x in sig:
                        for pos in reg.positions:
                            for b in range((x[0][1] - x[0][0]) + 1):
                                #significant['bed'].writerow(
                                #    ["chr" + self.chr, (x[0][0] + b) * self.bin_res,
                                #     (x[0][0] + b + 1) * self.bin_res, "chr" + self.chr + "_" + str(pos[0]) + "_"
                                #     + str(pos[1] + 1), str(min(list(-np.log10(x[1][b]))))])
                                significant['bed'].writerow(
                                    ["chr" + self.chr, str(pos[0]), str(pos[1] + 1), "chr" + self.chr, (x[0][0] + b) * self.bin_res,
                                     (x[0][0] + b + 1) * self.bin_res, str(min(list(-np.log10(x[1][b]))))])
                if 'gff' in formats:
                    if not self._save_sig_gff(significant['gff'], sig, False, reg.positions, reg.bin, remap, lo):
                        continue

        self._save_stats(rest_of_name, len(regions), r_pred, pred, h_fs)

    def _save_sig_gff(self, file, sig, single, positions, bin, remap, lo):
        if sig:
            rows = []
            region_strand = ''
            for pos in positions:
                if remap is not None:
                    rem = self._remap(pos[0], pos[1], lo)
                else:
                    rem = [[('chr' + self.chr, pos[0] + 1, '+')],
                           [('chr' + self.chr, pos[1] + 1, '+')]]
                if rem is not None:
                    rows.append(
                        [rem[0][0][0], 'LongRanger', 'region', rem[0][0][1], rem[1][0][1], 0, '.',
                         '.', 'Parent=bin_' + self.chr + '_' + str(
                            bin)])
                    region_strand = rem[0][0][2]
            if len(set([row[0] for row in rows])) != 1 or len(set([row[6] for row in rows])) != 1:
                return False
            mi = min([row[3] for row in rows])
            ma = max([row[4] for row in rows])
            for x in sig:
                for b in range((x[0][1] - x[0][0]) + 1):
                    if remap is not None:
                        rem = self._remap(((x[0][0] + b) * self.bin_res) + 1,
                                          ((x[0][0] + b + 1) * self.bin_res), lo)
                    else:
                        rem = [[('chr' + self.chr, ((x[0][0] + b) * self.bin_res) + 1, '+')],
                               [('chr' + self.chr, ((x[0][0] + b + 1) * self.bin_res), '+')]]
                    if rem is not None and rem[0][0][0] == rows[0][0] and rem[0][0][2] == region_strand:
                        if rem[0][0][1] < mi:
                            mi = rem[0][0][1]
                        if rem[1][0][1] > ma:
                            ma = rem[1][0][1]
                        if single:
                            rows.append(
                                rows[0][:2] + ['prediction', rem[0][0][1], rem[1][0][1], str(-np.log10(x[1][b]))] +
                                rows[0][6:])
                        else:
                            rows.append(
                                rows[0][:2] + ['prediction', rem[0][0][1], rem[1][0][1],
                                               str(min(list(-np.log10(x[1][b]))))] +
                                rows[0][6:])
            if ma - mi < self.num_bins * self.bin_res * 2 * 4:
                file.writerows([rows[0][:2] + ['group', mi, ma] + rows[0][5:-1] + [
                    'ID=bin_' + self.chr + '_' + str(bin)]] + rows)
        return True

    def _save_stats(self, rest_of_name, len_regions, r_pred, pred, h_fs):
        logger.debug('Saving stats file')
        stats = open(
            self.stats_folder + '/' + self.regions_name + '-' + "-".join(h_fs) + '-stats' + self.chr + '-' + str(
                self.num_bins) +
            rest_of_name + '.txt', 'w')
        stats.write('Number of regions: ' + str(len_regions) + '\n')
        stats.write('\twith predictions: ' + str(r_pred) + '\n')
        stats.write('\twithout predictions: ' + str(len_regions - r_pred) + '\n')
        stats.write('Number of all predictions: ' + str(pred) + '\n')
        try:
            stats.write('\tper region(all): ' + str(float(pred) / float(len_regions)) + '\n')
        except ZeroDivisionError:
            stats.write('\tper region(all): ' + str(0) + '\n')
        try:
            stats.write('\tper region(with predictions): ' + str(float(pred) / float(r_pred)) + '\n')
        except ZeroDivisionError:
            stats.write('\tper region(with predictions): ' + str(0) + '\n')
        stats.write('Threshold: ' + str(self.threshold) + '\n')

    def _remap(self, pos_s, pos_e, lo):
        new_s = lo.convert_coordinate('chr' + self.chr, pos_s)
        new_e = lo.convert_coordinate('chr' + self.chr, pos_e)
        if new_s and new_e:
            if new_s[0][2] != new_e[0][2] or new_s[0][0] != new_e[0][0]:
                logger.debug(
                    'Did not remap ' + str([self.chr, pos_s, pos_e]) + ', result was: ' + str(new_s) + ' ' + str(new_e))
                return None
            elif new_s[0][2] == '-':
                return new_e, new_s
            else:
                return new_s, new_e
        else:
            logger.debug(
                'Did not remap ' + str([self.chr, pos_s, pos_e]) + ', result was: ' + str(new_s) + ' ' + str(new_e))
            return None

    def run(self, region_file, num_regs, section, plotting, single_sig, remap, regs_algorithm,
            stat_formats):
        """
        Runs the visualization according to provided arguments

        :param region_file: Path to file with regions to extract significant contacts for
        :param num_regs: number of regions to plot for - optional
        :param section: section of chromosome to plot for - optional
        :param plotting: if there should be plots generated and with what - mpl or rpy2
        :param single_sig: if significant contact files should be generated for single maps too
        :param remap: assemblies to map from and to the significant contact positions in gff.
        :param regs_algorithm: algorithm for regions extraction
        :param stat_formats: formats to output statistics (significant contacts)
        """
        pck_filename_base = self.pickled_folder + '/' + self.chr + '-' + str(self.num_bins) + 'x' + str(
            self.bin_res) + 'bp-'

        if os.path.exists(pck_filename_base + 'regs-' + self.regions_name + '-' + '-'.join(
                [hicmap.get_name() for hicmap in self.maps])):
            logger.info('Loading region data from pickled file')
            with open(pck_filename_base + 'regs-' + self.regions_name + '-' + '-'.join(
                    [hicmap.get_name() for hicmap in self.maps]), 'rb') as fp:
                regs = pickle.load(fp)

        else:
            logger.info('Extracting regions '+ self.regions_name)
            regs = self.get_regions(region_file, regs_algorithm)

            for hicmap in self.maps:
                fit_files = glob(
                    "-[0-9]*x".join(
                        pck_filename_base.split('-' + str(self.num_bins) + 'x')) + 'fits-' + hicmap.get_name())
                logger.debug('Fit files found: ' + str(fit_files))
                best = 0
                if fit_files:
                    p = re.compile(r'-(?P<n_b>[0-9]*)x')
                    best = max([int(p.search(x).group('n_b')) for x in fit_files])
                    logger.info('Loading ' + str(
                        best * 2 + 1) + ' Weibull fits for Map ' + hicmap.get_name() + ' from pickled file')
                    with open(("-" + str(best) + "x").join(
                            pck_filename_base.split('-' + str(self.num_bins) + 'x')) + 'fits-' + hicmap.get_name(),
                              'rb') as fp:
                        hicmap.fits = pickle.load(fp)

                if self.num_bins > best:
                    hicmap.get_Weibull_fits(self.num_bins)
                    with open(pck_filename_base + 'fits-' + hicmap.get_name(), 'wb') as fp:
                        pickle.dump(hicmap.fits, fp)

                mean_files = glob(
                    "-[0-9]*x".join(
                        pck_filename_base.split('-' + str(self.num_bins) + 'x')) + 'means-' + hicmap.get_name())
                logger.debug('Mean files found: ' + str(mean_files))
                best = 0
                if mean_files:
                    p = re.compile(r'-(?P<n_b>[0-9]*)x')
                    best = max([int(p.search(x).group('n_b')) for x in mean_files])
                    logger.info(
                        'Loading ' + str(best * 2 + 1) + ' means for Map ' + hicmap.get_name() + ' from pickled file')
                    with open(("-" + str(best) + "x").join(
                            pck_filename_base.split('-' + str(self.num_bins) + 'x')) + 'means-' + hicmap.get_name(),
                              'rb') as fp:
                        hicmap.means = pickle.load(fp)

                if self.num_bins > best:
                    hicmap.get_means(self.num_bins)
                    with open(pck_filename_base + 'means-' + hicmap.get_name(), 'wb') as fp:
                        pickle.dump(hicmap.means, fp)

                logger.info('Setting intensities for regions from Map ' + hicmap.get_name())
                self.set_regs_intensities(regs, hicmap)

                logger.info('Calculating the pvalues and weighted values for Map ' + hicmap.get_name())
                self.calc(regs, hicmap)

            with open(pck_filename_base + 'regs-' + self.regions_name + '-' + '-'.join(
                    [hicmap.get_name() for hicmap in self.maps]), 'wb') as fp:
                pickle.dump(regs, fp)

        logger.info('Saving predictions and statistics for regions...')
        self.save_stats_and_sigs(regs, single_sig, remap, stat_formats,
                                 [hicmap.get_name() for hicmap in self.maps])

        logger.info('Extracting finished.')

        if plotting is not None:
            logger.debug('Getting to plotter')
            from .visualize import Plotter
            p = Plotter(self.regions_name + '-' + '-'.join(
                [hicmap.get_name() for hicmap in self.maps]), self.chr, num_regs, self.num_bins, self.threshold,
                        self.bin_res)
            p.run(section, self.pickled_folder, self.figures_folder, plotting)
