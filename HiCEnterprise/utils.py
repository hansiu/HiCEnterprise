"""
Utilities for both domains and regions.
"""

import os
import numpy as np
import logging

logger = logging.getLogger('utils')
logging.basicConfig(level=logging.INFO)


def create_folders(folders):
    """
    Creates necessary folders if they are not there already and returns the full paths.
    """
    folder_paths = []
    for folder in folders:
        folder_path = os.path.abspath(folder)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        folder_paths.append(folder_path)
    return folder_paths


def load_hicmap(hic_folder, filename):
    """
    Loads HiC Map from given folder and filename and checks the symmetry.
    Maps filenames should be in 'mtx-N-N.npy' form, where N is the chromosome name.
    Files should contain one HiC map in numpy format.
    """
    name = hic_folder + '/' + filename
    logger.info('Loading Hic-Map: ' + name)
    hicmap = np.load(os.path.abspath(name))
    if not np.allclose(hicmap.transpose(), hicmap):
        logger.warn('This map is not too symmetric')
    return hicmap
