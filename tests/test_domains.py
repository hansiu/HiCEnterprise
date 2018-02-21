import logging
import pytest
import os
import numpy as np

from HiCEnterprise.domains.extract import Extractor

logging.disable(logging.CRITICAL)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


class TestExtractor:
    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.tmpdir = tmpdir.strpath

    @staticmethod
    def create():
        domain_file = TEST_DIR + '/test_files/doms/sherpa-tads'
        chrom = '1'
        bin_res = 150000
        sherpa_lvl = 1
        hic_folder = TEST_DIR + '/test_files/maps/'
        threshold = 0.01
        plot_title = "Plot"
        ticks_separation = 0
        hic_color = 'Reds'
        interactions_color = "Blues"

        e = Extractor(domain_file, chrom, bin_res, sherpa_lvl, hic_folder, threshold, plot_title, ticks_separation,
                      hic_color, interactions_color)
        return e

    def test_init(self):
        e = self.create()
        assert e.domains_name == 'sherpa-tads'
        assert e.chr == '1'
        assert e.bin_res == 150000
        assert e.sherpa_level == 1
        assert isinstance(e.domains, list)
        assert e.domains[0] == [0, 0]
        assert isinstance(e.hicmap, np.ndarray)
        assert e.hic_folder == TEST_DIR + '/test_files/maps'
        assert e.hic_name == 'maps'
        assert e.threshold == 0.01
        assert e.plot_title == "Plot"
        assert e.ticks_separ == 0
        assert e.hic_color == 'Reds'
        assert e.inter_color == "Blues"
        assert e.hicmap.shape == (100, 100)

    def test_check_overlap_bad(self):
        e = self.create()
        e.domains = [[0, 3], [2, 5]]
        with pytest.raises(SystemExit) as pytest_e:
            e._check_overlap()
        assert pytest_e.type == SystemExit
        assert pytest_e.value.code == 1

    def test_create_domain_matrix(self):
        e = self.create()
        domain_matrix = e.create_domain_matrix()
        assert isinstance(domain_matrix, np.ndarray)
        print("domains", e.domains)
        assert domain_matrix.shape == (len(e.domains), len(e.domains))

    def test_calc(self):
        # TODO failing now - fix when the decision is made about the sum of contacts and symm
        e = self.create()
        sigs, corr_sigs = e.calc(e.create_domain_matrix())
        assert len(sigs) == 100
        assert len(corr_sigs) == 100
        assert sigs[0] == (5, 6, pytest.approx(0))
        assert corr_sigs[1] == (5, 7, pytest.approx(0))

    def test_plotting(self):
        e = self.create()
        e.run(self.tmpdir + '/stats/', True, self.tmpdir + '/figures/')
        assert os.path.exists(self.tmpdir + '/stats/sherpa-tads-maps-stats1-0_01.txt')
        assert os.path.exists(self.tmpdir + '/figures/maps-sherpa-tads.png')
