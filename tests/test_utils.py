import pytest
import os
import logging

from HiCEnterprise.utils import create_folders, load_hicmap

logging.disable(logging.CRITICAL)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

class Test_create_folders():
    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.tmpdir = tmpdir.strpath

    def test_created(self):
        paths = create_folders([self.tmpdir + '/stats/'])
        assert paths == [self.tmpdir + '/stats']
        assert os.path.exists(paths[0])


class Test_load_hicmap():
    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.tmpdir = tmpdir.strpath

    def test_bad_load(self):
        with pytest.raises(IOError):
            h = load_hicmap(self.tmpdir + '/some/path/reallyrandom', 'mtx-blah-blah.npy')

    def test_good_load(self):
        h = load_hicmap(TEST_DIR + '/test_files/maps/fbd', 'mtx-12-12.npy')
        assert h.shape == (100, 100)
