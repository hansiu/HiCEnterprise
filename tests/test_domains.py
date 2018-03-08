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
    def create(distribution='hypergeom'):
        domain_file = TEST_DIR + '/test_files/doms/sherpa-tads'
        chrom = '1'
        bin_res = 150000
        sherpa_lvl = 1
        hic_map = TEST_DIR + '/test_files/maps/mtx-1-1.npy'
        threshold = 0.01
        plot_title = "Plot"
        ticks_separation = 0
        hic_color = 'Reds'
        interactions_color = "Blues"
        hic_name = 'test_maps'

        e = Extractor(domain_file, chrom, bin_res, sherpa_lvl, hic_map, threshold, plot_title, ticks_separation,
                      hic_color, interactions_color, distribution, hic_name)
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
        assert e.hic_path == TEST_DIR + '/test_files/maps/mtx-1-1.npy'
        assert e.threshold == 0.01
        assert e.plot_title == "Plot"
        assert e.ticks_separ == 0
        assert e.hic_color == 'Reds'
        assert e.inter_color == "Blues"
        assert e.hicmap.shape == (100, 100)
        assert e.distribution == 'hypergeom'
        assert e.hic_name == 'test_maps'

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

    def test_calc_hypergeom(self):
        e = self.create()
        sigs, corr_sigs = e.calc(e.create_domain_matrix())
        assert len(sigs) == 100
        assert len(corr_sigs) == 100
        assert sigs[0] == (5, 6, pytest.approx(0))
        assert corr_sigs[1] == (5, 7, pytest.approx(0))

    def test_calc_negbinom(self):
        e = self.create('negbinom')
        sigs, corr_sigs = e.calc(e.create_domain_matrix())
        assert len(sigs) == 70
        assert len(corr_sigs) == 59
        assert sigs[0] == (5, 17, pytest.approx(0, abs=1e-8))
        assert corr_sigs[1] == (5, 31, pytest.approx(0.0066, abs=1e-3))

    def test_calc_poisson(self):
        e = self.create('poisson')
        sigs, corr_sigs = e.calc(e.create_domain_matrix())
        assert len(sigs) == 106
        assert len(corr_sigs) == 94
        assert sigs[0] == (5, 7, pytest.approx(0))
        assert corr_sigs[1] == (5, 17, pytest.approx(0))


    def test_plotting(self):
        e = self.create()
        e.run(self.tmpdir + '/stats/', True, self.tmpdir + '/figures/')
        assert os.path.exists(self.tmpdir + '/stats/sherpa-tads-test_maps-stats1-0_01-hypergeom.txt')
        assert os.path.exists(self.tmpdir + '/figures/test_maps-sherpa-tads-hypergeom.png')
