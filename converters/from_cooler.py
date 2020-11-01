##### Save cooler files as numpy arrays
#wget ftp://cooler.csail.mit.edu/coolers/hg19/Rao2014-HUVEC-MboI-allreps-filtered.1000kb.cool
#cd HiCEnterprise/
#python3.6 converters/from_cooler.py -o your/path/here/Rao1000kb/ -i your/path/to/file/Rao2014-HUVEC-MboI-allreps-filtered.1000kb.cool
# OR
#python3.6 converters/from_cooler.py -o your/path/here/Rao1000kb/ -i your/path/to/file/Rao2014-HUVEC-MboI-allreps-filtered.1000kb.cool -c 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X

import argparse
chrs = [str(x) for x in range(1, 23)] + ['X']

parser = argparse.ArgumentParser(prog='HiCE convert from hdf5')
parser.add_argument('-o', '--output_directory', help='Directory to which to output the mtx-C-C.npy files', type=str,
                           required=True)
parser.add_argument('-i','--input_file', help='Input hdf5 file to convert to numpy format', type=str, required=True)
parser.add_argument('-c','--chromosome_names', help='The chromosome numerings (in the format of 1 2 3 ... X, or I II III...'
                                                    ' without the "chr" preceeding) for the analysis. Default is human,'
                                                    ' from 1 to 22 including X, without Y. Provide as a list.', nargs='+', default=chrs)

if __name__ == '__main__':
    import cooler
    import numpy as np
    import os

    args = parser.parse_args()
    save_folder = os.path.abspath(args.output_directory)
    os.makedirs(save_folder,exist_ok=True)
    input_file = os.path.abspath(args.input_file)
    chrs = args.chromosome_names

    file_format = 'mtx-{}-{}.npy'

    c = cooler.Cooler(input_file)
    for chr in chrs:
        print('now preparing: {}'.format(file_format.format(chr,chr)))
        save_path = '{}/{}'.format(save_folder, file_format.format(chr,chr))
        data = c.matrix().fetch('chr'+chr)
        data = np.nan_to_num(data)
        np.save(save_path, data)
