import logging
import pytest
import os
import numpy as np

from HiCEnterprise.regions.extract import Bin, HiCmap, Extractor
from HiCEnterprise.regions.visualize import Plotter

logging.disable(logging.CRITICAL)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

class TestBin():
    def test_init(self):
        positions = [(10, 11), (11, 12)]

        b = Bin(1, [(10, 11), (11, 12)])
        assert b.bin == 1
        assert b.positions == positions


class TestHiCmap():
    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.tmpdir = tmpdir.strpath

    def test_init(self):
        folder = self.tmpdir + '/some/path'
        filename = '/mtx-1-1.npy'

        h = HiCmap(self.tmpdir + '/some/path/', '1')
        assert h.folder == folder
        assert h.filename == filename

    def test_get_name(self):
        folder_name = 'path'

        h = HiCmap(self.tmpdir + '/some/path', '1')
        assert h.get_name() == folder_name

    def test_bad_load(self):
        h = HiCmap(self.tmpdir + '/some/path/reallyrandom', '1')
        with pytest.raises(IOError):
            h.load()

    def test_good_load(self):
        h = HiCmap(TEST_DIR + '/test_files/maps/fbd', '12')
        assert not h.loaded
        assert h.map is None
        h.load()
        assert h.loaded
        assert h.map is not None

    def test_get_means(self):
        expected_means = {-2:70.018349,-1:1,0:1,1:1,2:70.018349}

        h = HiCmap(TEST_DIR + '/test_files/maps/fbd', '12')
        assert not h.means
        h.get_means(2)
        assert len(h.means)==len(expected_means)
        assert h.means == pytest.approx(expected_means)

    def test_get_Weibull_fits(self):
        expected_fits = {-2: (2.166370, 2.259363, 0, 66.175952), -1: None, 0: None, 1: None, 2: (2.166370, 2.259363, 0, 66.175952)}

        h = HiCmap(TEST_DIR + '/test_files/maps/fbd', '12')
        assert not h.fits
        h.get_Weibull_fits(2)
        assert len(h.fits) == len(expected_fits)
        for fit in h.fits.keys():
            assert h.fits[fit] == pytest.approx(expected_fits[fit])


class TestExtractor():
    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.tmpdir = tmpdir.strpath

    def create(self):
        region = 'Fetal_brain'
        hicmaps = [TEST_DIR + '/test_files/maps/fbd']
        chrom = '12'
        resolution = 10000
        bins = 5
        threshold = 0.02
        p_folder, f_folder, s_folder = self.tmpdir + '/pickles', self.tmpdir + '/figures', self.tmpdir + '/stats'

        e = Extractor(region, hicmaps, chrom, resolution, bins, threshold, p_folder,
                      f_folder, s_folder)

        return (e)

    def test_init(self):
        e = self.create()
        assert e.regions_name == 'Fetal_brain'
        assert len(e.maps) == 1
        assert isinstance(e.maps[0],HiCmap)
        assert e.chr == '12'
        assert e.bin_res == 10000
        assert e.num_bins == 5
        assert e.threshold == 0.02

    def test_get_regions_all(self):
        e = self.create()
        regs_all_fasta = e.get_regions(TEST_DIR + '/test_files/enhseq/Fetal_brain.fasta')
        regs_all_bed = e.get_regions(TEST_DIR + '/test_files/enhseq/Fetal_brain.bed')
        assert len(regs_all_fasta) == 10 == len(regs_all_bed)
        assert isinstance(regs_all_fasta[0], Bin) == isinstance(regs_all_bed[0], Bin)
        assert regs_all_fasta[5].bin == 67 == regs_all_bed[5].bin
        assert regs_all_fasta[5].positions[0] == regs_all_fasta[6].positions[0] == regs_all_bed[5].positions[0] == \
               regs_all_bed[6].positions[0]
        assert regs_all_fasta[6].positions[1] == [685290, 685980] == regs_all_bed[6].positions[1]
        assert len(regs_all_fasta[8].positions) == 2 == len(regs_all_bed[8].positions)

    def test_get_regions_fasta_one(self):
        e = self.create()
        regs_one_fasta = e.get_regions(TEST_DIR + '/test_files/enhseq/Fetal_brain.fasta', 'one')
        regs_one_bed = e.get_regions(TEST_DIR + '/test_files/enhseq/Fetal_brain.bed', 'one')
        assert len(regs_one_fasta) == 10 == len(regs_one_bed)
        assert isinstance(regs_one_fasta[0], Bin) == isinstance(regs_one_bed[0], Bin)
        assert regs_one_fasta[6].bin == 68 == regs_one_bed[6].bin
        assert regs_one_fasta[6].positions == [[685290, 685980]] == regs_one_bed[6].positions
        assert len(regs_one_fasta[8].positions) == 2 == len(regs_one_bed[8].positions)

    def test_set_regs_intensities(self):
        e = self.create()
        regs = e.get_regions(TEST_DIR + '/test_files/enhseq/Fetal_brain.bed')
        e.set_regs_intensities(regs, e.maps[0])
        for reg in regs:
            assert len(reg.intensities) == len(e.maps)
            assert len(reg.intensities[0]) == 11
        assert regs[0].intensities[0][2] == pytest.approx(49.47158)

    def test_calc(self):
        e = self.create()
        regs = e.get_regions(TEST_DIR + '/test_files/enhseq/Fetal_brain.bed')
        hicmap = e.maps[0]
        hicmap.get_Weibull_fits(e.num_bins)
        hicmap.get_means(e.num_bins)
        e.set_regs_intensities(regs, hicmap)
        e.calc(regs, hicmap)
        for reg in regs:
            assert len(reg.weighted) == len(reg.pvalues) == len(reg.corrected_pvalues) == 1
            assert len(reg.weighted[0]) == len(reg.pvalues[0]) == len(reg.corrected_pvalues[0]) == 11

        assert np.any(regs[2].corrected_pvalues[0][1] < e.threshold)

        assert list(regs[0].pvalues[0]) == [1.0, pytest.approx(0.7064177), pytest.approx(0.5278167),
                                            pytest.approx(0.9939315), 1.0, 1.0, 1.0, pytest.approx(0.0344187),
                                            pytest.approx(0.0028666697), pytest.approx(0.51194545),
                                            pytest.approx(0.0144503)]

    def test_save_stats_and_sigs(self):
        # txt, gff, bed
        e = self.create()
        regs = e.get_regions(TEST_DIR + '/test_files/enhseq/Fetal_brain.bed')
        hicmap = e.maps[0]
        hicmap.get_Weibull_fits(e.num_bins)
        hicmap.get_means(e.num_bins)
        e.set_regs_intensities(regs, hicmap)
        e.calc(regs, hicmap)
        single_sig = False
        remap = None
        stat_formats = ['txt', 'gff', 'bed']
        e.save_stats_and_sigs(regs, single_sig, remap, stat_formats, [hicmap.get_name()])
        for form in stat_formats:
            assert os.path.exists(e.stats_folder + '/Fetal_brain-' + "-".join(
                [hicmap.get_name()]) + '-significant12-5x10000bp-0_02.' + form)

    def test_plotting(self):
        e = self.create()
        region_file = TEST_DIR + '/test_files/enhseq/Fetal_brain.bed'
        num_regs = 3
        section = '340000:690000'
        plotting = 'mpl'
        single_sig = False
        remap = None
        regs_algorithm = 'all'
        stat_formats = ['txt']
        e.run(region_file, num_regs, section, plotting, single_sig, remap, regs_algorithm, stat_formats)

        assert os.path.exists(self.tmpdir + '/figures/Fetal_brain-fbd-12-340000-690000-3.pdf')
