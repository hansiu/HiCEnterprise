"""
Utilities for both domains and regions.
"""

import os
import numpy as np
import logging
import scipy.ndimage

logger = logging.getLogger('utils')
logging.basicConfig(level=logging.INFO)


def clip_and_blur(arr, stddevs=5, blur=1):
    """
    Clips and blurs the matrix as needed. By Krzysiek Krolak.
    """
    arr = np.ma.masked_invalid(arr)
    mean = np.mean(arr)
    stddev = np.var(arr) ** 0.5
    np.clip(arr, 0, mean + stddevs * stddev, out=arr)
    arr = np.ma.filled(arr, 0)
    scipy.ndimage.gaussian_filter(arr, blur, output=arr)
    np.clip(arr, mean * 0.01, mean + stddevs * stddev, out=arr)
    return arr


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
        logger.warning('This map is not too symmetric')
    return hicmap
