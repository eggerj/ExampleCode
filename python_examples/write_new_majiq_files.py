#!/usr/bin/env python3

import argparse
import numpy as np
import pickle


def write_sample_data(majiq_file, lsv_dict, junc_dict, out_dir):

    new_data = {}
    with np.load(majiq_file) as np_fh:
        sample_name = np_fh['meta'][0][0].decode('utf-8').split('.')[0]
        meta = np_fh['meta']
        junc_info = np_fh['junc_info']
        coverage = np_fh['coverage']
        lsv_types = np_fh['lsv_types']
        # Create new dictionaries for storing filtered data
    new_data['meta'] = meta
    new_data['junc_info'] = []
    new_data['coverage'] = []
    new_data['lsv_types'] = []
    print('Writing new majiq file for', sample_name)

    # Add junctions and coverage info to new array dict if junction remains after filtering
    for j,cov in zip(junc_info, coverage):
        if j[0] in lsv_dict:
            j_key = (j[0], j[1], j[2]) # (LSV ID, pos1, pos2)
            if j_key in junc_dict:
                new_data['junc_info'] = new_data['junc_info'] + [j,]
                new_data['coverage'] =  new_data['coverage'] + [cov,]
    
    # Add LSV type info if LSV remains after filtering
    for lsv in lsv_types:
        if lsv[0] in lsv_dict:
            new_data['lsv_types'] =  new_data['lsv_types'] + [lsv,]
    
    # Cast new_data listings to numpy arrays
    new_data['junc_info'] = np.array(new_data['junc_info'])
    new_data['coverage'] = np.array(new_data['coverage'])
    new_data['lsv_types'] = np.array(new_data['lsv_types'])
    fn = out_dir + '/' + sample_name
    np.savez(fn, **new_data)


def load_numpy_arrays(out_dir):

    fn = out_dir + '/tmp/lsv_dict.pkl'
    with open(fn, 'rb') as f:
        lsv_dict = pickle.load(f)

    fn = out_dir + '/tmp/junc_dict.pkl'
    with open(fn, 'rb') as f:
        junc_dict = pickle.load(f)

    return lsv_dict, junc_dict

if __name__ == "__main__":

    # Program arguments
    parser = argparse.ArgumentParser(description='Write new MAJIQ splicing graphs across sample set')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--majiq_file', type=str, help='.majiq file from MAJIQ builder.', required=True)
    requiredNamed.add_argument('--out_dir', type=str, help='Full path to output directory to store simplified majiq files.', required=True)
    args = parser.parse_args()

    # Get arguments
    majiq_file = args.majiq_file
    out_dir = args.out_dir

    lsv_dict, junc_dict = load_numpy_arrays(out_dir)

    write_sample_data(majiq_file, lsv_dict, junc_dict, out_dir)

