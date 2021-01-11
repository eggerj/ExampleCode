#!/usr/bin/env python3

import numpy as np
import argparse

def get_sample_counts(files):

    # Create dictionary to store LSV sample counts and list of sample names
    lsv_dict = {}
    samples = []

    with open(files, 'r') as fh:
        sample_files = [l.strip() for l in fh.readlines()]

    for sf in sample_files:
        with np.load(sf.strip()) as np_fh:
            junc_info = np_fh['junc_info']
            sample_name = str(np_fh['meta'][0][0].decode('utf-8').split('.')[0].split('_')[1])
            samples.append(sample_name)
        print('Loading Junctions For', sample_name)        
        for j in junc_info:
            lsv = str(j[0].decode('utf-8'))
            read_count = j[3]
            if lsv not in lsv_dict:
                lsv_dict[lsv] = {}
            if sample_name not in lsv_dict[lsv]:
                lsv_dict[lsv][sample_name] = float(read_count)
            else:
                lsv_dict[lsv][sample_name] = lsv_dict[lsv][sample_name] + float(read_count) 

    return lsv_dict,samples        
    
# Output LSV x sample matrix of read counts
def write_counts_matrix(lsv_dict, samples, fn):

    with open(fn, 'w') as fh:
        fh.write(',' + ','.join(samples) + '\n')
        n = 0 # Counter for debugging
        for lsv,counts in lsv_dict.items():
            fh.write(lsv)
            for s in samples:
                fh.write(',' + str(counts[s]))
            fh.write('\n')
            n += 1
            #if n > 10:
            #    break

if __name__ == "__main__":

    # Arguments list
    parser = argparse.ArgumentParser(description='Program to get distribution of total reads covering each LSV from each sample.')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--majiq_files', type=str, help='File containing all (simplified) majiq files with full paths', required=True)
    requiredNamed.add_argument('--outfile', type=str, help='Name (and path) of outfile to contain sample by LSV read count matrix (csv file).', required=True)
    args = parser.parse_args()

    # Get arguments
    majiq_files = args.majiq_files
    outfile = args.outfile

    # Example of use from AML dataset
    #majiq_files = '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/beatAML_splicing/majiq_simplify_v5/majiq_simplified_files_v5.txt'
    #outfile = '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/beatAML_splicing/lsv_data_v5/lsv_read_counts.x484.csv'

    lsv_dict, samples = get_sample_counts(majiq_files)

    write_counts_matrix(lsv_dict, samples, outfile)


