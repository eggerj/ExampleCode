#!/usr/bin/env python3

'''
This program is used to simplify alternative splicing graphs created by the MAJIQ 
 framework prior to LSV quantification.

This is part 1 of a two part series of scripts ran using a Slurm script on Exacloud.
 Part 1 loads in .majiq files from MAJIQ (they are actually .npz files from numpy) and
  evaluates the coverage of each junction across the samples. If junctions do not
  have meaninful read counts across a reasonable number of samples they are removed.
'''


import argparse
import numpy as np
from collections import Counter
import pickle

# Function to create dictionaries of sample counts for all LSV junctions
def make_junction_dictionaries(files):

    # {Junction:[sample read counts]} --> key: (LSV, pos1, pos2)
    junc_dict = {}

    with open(files, 'r') as fh:
        for sample in fh.readlines():
            with np.load(sample.strip()) as np_fh:
                junc_info = np_fh['junc_info']
                sample_name = np_fh['meta'][0][0].decode('utf-8').split('.')[0]
            print('Loading Junctions For', sample_name)
            for j in junc_info:
                j_key = (j[0], j[1], j[2]) # (LSV ID, pos1, pos2)
                read_count = j[3]
                if j_key not in junc_dict:
                    junc_dict[j_key] = [read_count,]
                else:
                    junc_dict[j_key] = junc_dict[j_key] + [read_count,]

    print('Total Junctions Before Filtering:', str(len(junc_dict.keys())))

    return junc_dict

# Function to remove junctions not meeting read thresholds
def filter_junctions(j_dict, min_reads, min_exp, ss_min_reads):

    new_j_dict = {}

    for j,counts in j_dict.items():
        read_counts = np.array(counts)
        if (np.quantile(read_counts, min_exp) >= min_reads) or (np.max(read_counts) >= ss_min_reads):
            new_j_dict[j] = read_counts

    print('Total Junctions After Filtering:', str(len(new_j_dict.keys())))

    return new_j_dict

# Function to remove LSVs not meeting read thresholds or if less than two junctions 
#  remain after filtering
def filter_LSVs(j_dict, min_reads, min_exp):

    lsv_dict = {} 

    # Fill initial LSV dictionary with list of lists of read counts (keep separate for each junction)
    for j,counts in j_dict.items():
        lsv = j[0]
        if lsv not in lsv_dict:
            lsv_dict[lsv] = [np.array(counts),]
        else:
            lsv_dict[lsv] = lsv_dict[lsv] + [np.array(counts),]  

    # Create new LSV dictionary containing LSVs having > 1 junction and total read counts meeting filtering threshold
    filtered_lsv_dict = {}
    for lsv,counts in lsv_dict.items():
        if len(counts) > 1:
            read_totals = np.array(counts).sum(axis=0)
            if np.quantile(read_totals, min_exp) >= min_reads:
                filtered_lsv_dict[lsv] = counts

    return filtered_lsv_dict  

# OLD FUNCTION --> NO LONGER USED (USE WRITE SCRIPT ON SLURM)
# Create new .npz files mimicking .majiq files  
def write_sample_data(majiq_files, lsv_dict, junc_dict, out_dir):

    # Parse list of majiq files and create new files based on filtered LSV set                      
    with open(majiq_files, 'r') as fh:
        for sample in fh.readlines():
            new_data = {}
            with np.load(sample.strip()) as np_fh:                
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

# NEW FUNCTION
#  Store all data in pkl files and write new files in parallel using Slurm
def save_numpy_array(lsv_dict, junc_dict, out_dir):

    fn = out_dir + '/tmp/lsv_dict.pkl'
    with open(fn, 'wb') as f:
        pickle.dump(lsv_dict, f, pickle.HIGHEST_PROTOCOL)

    fn = out_dir + '/tmp/junc_dict.pkl'
    with open(fn, 'wb') as f:
        pickle.dump(junc_dict, f, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":

    # Parameter arguments
    parser = argparse.ArgumentParser(description='Simplify MAJIQ splicing graphs across sample set')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--majiq_files', type=str, help='Line by line listing of .majiq files from MAJIQ builder (must include full paths).  \
                                All samples used during build step must be included.', required=True)
    requiredNamed.add_argument('--junction_min_reads', type=int, help='Minimum number of reads covering junction to keep in LSV.', required=True)
    requiredNamed.add_argument('--junction_min_experiments', type=float, help='Minimum number (ratio) of samples for applying junction read filter.', required=True)
    requiredNamed.add_argument('--lsv_min_reads', type=int, help='Minimum number of total reads from all junctions to keep LSV.', required=True)
    requiredNamed.add_argument('--lsv_min_experiments', type=float, help='Minimum number (ratio) of samples for applying LSV read total filter.', required=True)
    requiredNamed.add_argument('--single_sample_junction_min_reads', type=int, help='Keep junctions meeting threshold even if only in single sample. (High count may be of interest)', required=True)
    requiredNamed.add_argument('--out_dir', type=str, help='Name (and path) of output directory to store simplified majiq files.', required=True)
    args = parser.parse_args()


    # Get arguments
    majiq_files = args.majiq_files
    junc_min_reads =float(args.junction_min_reads)
    junc_min_experiments = args.junction_min_experiments
    lsv_min_reads = float(args.lsv_min_reads)
    lsv_min_experiments = args.lsv_min_experiments
    ss_junc_min_reads = float(args.single_sample_junction_min_reads)
    out_dir = args.out_dir

    # For log
    with open(out_dir + '/log_file.txt', 'w') as fh:
        for arg in vars(args):
            fh.write(str(arg) + ': ' + str(getattr(args, arg)) + '\n')


    # Make junction dictionaries 
    #    {Junction:[sample read counts]} --> key: (LSV, pos1, pos2)
    #    Can get remainder of sample arrays from remaining junctions after filtering 
    junction_read_counts = make_junction_dictionaries(majiq_files)

    # Filter junctions having read count below threshold for specified ratio of samples    
    filtered_junctions = filter_junctions(junction_read_counts, junc_min_reads, junc_min_experiments, ss_junc_min_reads)

    # Compute LSV total read counts & filter LSVs having less than 2 junctions and/or having read count below threshold for specified ratio
    filtered_LSVs = filter_LSVs(filtered_junctions, lsv_min_reads, lsv_min_experiments)
    print('Remaining LSVs:', len(filtered_LSVs.keys()))
    with open(out_dir + '/summary_file.txt', 'w') as fh:   
        fh.write('Total junctions before filtering: ' + str(len(junction_read_counts.keys())) + '\n')
        fh.write('Total junctions after filtering: ' + str(len(filtered_junctions.keys())) + '\n')
        fh.write('Remaining LSVs: ' + str(len(filtered_LSVs.keys())) + '\n')

    # Save lsv and junction dictionaries to numpy array for later
    save_numpy_array(filtered_LSVs, filtered_junctions, out_dir)


