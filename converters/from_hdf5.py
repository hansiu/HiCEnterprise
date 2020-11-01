#for the data in the format like:
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77565/suppl/GSE77565_FBD_IC-heatmap-res-100k.hdf5.gz
#if its gzipped, gunzip
#gunzip your/path/to/file/GSE77565_FBD_IC-heatmap-res-100k.hdf5.gz
#python3.6 converters/from_hdf5.py -o your/path/here/fbd100 -i your/path/to/file/GSE77565_FBD_IC-heatmap-res-100k.hdf5
#OR
#python3.6 converters/from_hdf5.py -o your/path/here/fbd100 -i your/path/to/file/GSE77565_FBD_IC-heatmap-res-100k.hdf5 -z -m "{23:'X'}"

import argparse
import ast

parser = argparse.ArgumentParser(prog='HiCE convert from hdf5')
parser.add_argument('-o', '--output_directory', help='Directory to which to output the mtx-C-C.npy files', type=str,
                           required=True)
parser.add_argument('-i','--input_file', help='Input hdf5 file to convert to numpy format', type=str, required=True)
parser.add_argument('-z','--start_at_zero', help='If in hdf5 file chromosomes are enumerated from 0, not 1, provide this'
                                                 ' flag to convert to 1-based numering.', action='store_true')
parser.add_argument('-m','--chromosome_mapping', help='Mapping of chromosome numbers to other, written as a string: '
                                                      '"{23:\'X\'}"', type = str, default='"{23:\'X\'}"')

if __name__ == '__main__':
    import h5py
    import os
    import numpy as np

    args = parser.parse_args()
    save_folder = os.path.abspath(args.output_directory)
    os.makedirs(save_folder, exist_ok=True)

    input_file = os.path.abspath(args.input_file)
    start_at_0 = args.start_at_zero
    file_format = 'mtx-{}-{}.npy'

    chr_mapping = ast.literal_eval(ast.literal_eval(args.chromosome_mapping))

    with h5py.File(input_file,'r') as f:
        chrms_inx = f['chromosomeIndex'][...]
        chrms_starts = list(f['chromosomeStarts'][...])
        heatmap = f['heatmap']
        size = heatmap.shape[0]

        for start,end in zip(chrms_starts,chrms_starts[1:]+[size]):
            current_chrom = chrms_inx[start]
            if start_at_0:
                current_chrom+=1
            if current_chrom in chr_mapping.keys():
                current_chrom = chr_mapping[current_chrom]
            heatmap = f['heatmap'][start:end,start:end]
            save_path = '{}/{}'.format(save_folder, file_format.format(current_chrom,current_chrom))
            np.save(save_path, heatmap)
