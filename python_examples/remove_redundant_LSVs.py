#!/usr/bin/env python3

'''
This program is for identifying potential redundant LSVs that could be removed prior to SVR formulation.
  Redundant LSVs can be optionally removed, but are usually kept for reference (can be removed later).
'''


import argparse
import operator
import numpy as np

def load_sample_data(fn):

    with open(fn, 'r') as fh:
        sample_files = [l.strip() for l in fh.readlines()]

    # Try this: LSV dictionary of dictionaries
    #  {LSV_ID: {sample:line from tsv}}
    lsv_dict = {}

    sample_names = []
    for sf in sample_files:
        sample = sf.split('/')[-1].split('.')[0].split('_')[1]
        sample_names.append(sample)
        print(sample)
        with open(sf, 'r') as tsv_fh:
            lines = [l.strip().split('\t') for l in tsv_fh.readlines()]
        for rec in lines[1:]: # skip header
            lsv = rec[2] # get LSV ID
            if lsv not in lsv_dict:
                lsv_dict[lsv] = {}
            lsv_dict[lsv][sample] = rec

    return lsv_dict, sample_names

def get_shared_LSVs(lsv_dict, samples, min_sample_ratio):

    shared_lsvs = {}

    for lsv,data in lsv_dict.items():
        lsv_ratio = float(len(data.keys())) / float(len(samples))
        if lsv_ratio < min_sample_ratio:
            foo = 1 #TODO: write missing LSVs to file
        else:
            shared_lsvs[lsv] = data

    print('Total LSVs:', len(lsv_dict.keys()))
    print('Total shared LSVs:', len(shared_lsvs.keys()))

    return shared_lsvs

def get_lsv_pairs(lsv_set):

    # Dictionary of dictionaries: {geneID: {lsv_id:rec}}
    source_LSVs = {}
    target_LSVs = {}

    # Collect LSVs that do not contain alt splice sites or intron retentions
    for rec in lsv_set:
        lsv_type = rec[5]
        alt5 = rec[6]
        alt3 = rec[7]        
        if lsv_type.split('|')[-1] != 'i' and alt5 != 'True' and alt3 != 'True':
            gene_id = rec[1]
            lsv_id = rec[2]
            if lsv_id.split(':')[1] == 's':
                if gene_id not in source_LSVs:
                    source_LSVs[gene_id] = {}
                source_LSVs[gene_id][lsv_id] = rec
            else: 
                if gene_id not in target_LSVs:
                    target_LSVs[gene_id] = {}
                target_LSVs[gene_id][lsv_id] = rec

    return source_LSVs, target_LSVs       

def search_rec(source_rec, target_rec):

    # Look for the following:
    #  1) Share a junction
    #  2) Share first and last exons if 2 junctions exist
    #  3) Share all exons if 3 or more junctions exist 

    sourceID = source_rec[2]
    targetID = target_rec[2]
    # Pull junctions from records
    source_juncs = source_rec[14].split(';')
    target_juncs = target_rec[14].split(';')
    # Pull exons from records
    source_exons = sorted([sorted(pair.split('-')) for pair in source_rec[15].split(';')])
    target_exons = sorted([sorted(pair.split('-')) for pair in target_rec[15].split(';')])

    junction_shared = False
    for sj in source_juncs:
        if sj in target_juncs:
            junction_shared = True
            break

    if (junction_shared) and len(source_juncs) == len(target_juncs):
        # If binary LSVs, only first and last exons have to match
        if len(source_juncs) == 2:
            if source_exons[0] == target_exons[0] and source_exons[-1] == target_exons[-1]:
                return True, False
        # If 3+ junctions, all exons have to match
        else:
            if source_exons == target_exons:
                return True, False
    # Find mutually exclusive exons?
    else:
        if len(source_exons) == 3 and len(target_exons) == 3 and len(source_juncs) == 2 and len(target_juncs) == 2:
            if source_exons[1] == target_exons[0] and source_exons[2] == target_exons[1]:
                return True, True

    return False, False

def find_redundant_pairs(sources, targets):

    redundant_pairs = []
    mutually_exclusive = {}

    for gene_id, source_LSVs in sources.items():
        if gene_id in targets:
            for source_lsv, source_rec in source_LSVs.items(): # source LSVs for current gene
                for target_lsv, target_rec in targets[gene_id].items(): # target LSVs for current gene
                    pair_found, mxe = search_rec(source_rec, target_rec)
                    if (pair_found):
                        redundant_pairs.append((source_lsv, target_lsv, mxe))
                    if (mxe):
                        mutually_exclusive[source_lsv] = 'mutually_exclusive'
                        mutually_exclusive[target_lsv] = 'mutually_exclusive'

    print('Redundant LSV Pairs:', len(redundant_pairs))

    return redundant_pairs, mutually_exclusive

def get_variance(lsv, lsv_dict):

    junc_vars = [v[4].split(';') for v in lsv_dict[lsv].values()]
    return np.array([float(v) for var in junc_vars for v in var])

def filter_LSV_pairs(lsv_dict, pairs, mut_excl):

    redundant_lsvs = {}
    filtered_lsvs = {}
    removed_lsvs = {}
    hold_pairs = []

    # Try this: Keep LSV having smallest average variance of expected PSI values across samples (average of junctions, average of samples)
    for pair in pairs:
        source,target,mxe = pair[0],pair[1],pair[2]
        # Check that one LSV is part of a MXE set and if so save for later
        if (mxe == False) and (source in mut_excl or target in mut_excl):
            hold_pairs.append(pair)
        # If source or target already removed, keep it that way and keep the other LSV 
        elif target in removed_lsvs:
            filtered_lsvs[source] = lsv_dict[source]  
            if source in redundant_lsvs:
                redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
            else:
                redundant_lsvs[source] = target
        elif source in removed_lsvs:
            filtered_lsvs[target] = lsv_dict[target]
            if target in redundant_lsvs:
                redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
            else:
                redundant_lsvs[target] = source
        # Else, select LSV with lowest variance between the two across samples
        elif np.mean(get_variance(source, lsv_dict)) > np.mean(get_variance(target, lsv_dict)): # Keep target LSV
            filtered_lsvs[target] = lsv_dict[target]
            removed_lsvs[source] = source
            if target in redundant_lsvs:
                redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
            else:
                redundant_lsvs[target] = source
        else: 
            filtered_lsvs[source] = lsv_dict[source]
            removed_lsvs[target] = target
            if source in redundant_lsvs:
                redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
            else:
                redundant_lsvs[source] = target

    # Take care of redundant LSVs inside of mutually exclusive exon pair
    for pair in hold_pairs:
        source,target = pair[0],pair[1]
        # Keep LSV that was kept to represent MXE and remove other LSV
        if source in mut_excl and source in filtered_lsvs:
            removed_lsvs[target] = target
            redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
        elif source in mut_excl and source not in filtered_lsvs:
            filtered_lsvs[target] = lsv_dict[target]
            removed_lsvs[source] = source
            if target in redundant_lsvs:
                redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
            else:
                redundant_lsvs[target] = source 
        elif target in mut_excl and target in filtered_lsvs:
            removed_lsvs[source] = source
            redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
        else:
            filtered_lsvs[source] = lsv_dict[source]
            removed_lsvs[target] = target
            if source in redundant_lsvs:
                redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
            else:
                redundant_lsvs[source] = target

    # Store remaining LSVs redundant pairs as 'NaN' and store remaining LSVs without a redundant pairing 
    for lsv in lsv_dict.keys():
        if lsv not in redundant_lsvs and lsv not in removed_lsvs:
            redundant_lsvs[lsv] = 'NaN'
            filtered_lsvs[lsv] = lsv_dict[lsv] 

    print('Remaining LSVs:', len(filtered_lsvs.keys()))

    return filtered_lsvs, redundant_lsvs

def store_redundant_pairs(lsv_dict, pairs):

    redundant_lsvs = {}

    for pair in pairs:
        source, target = pair[0], pair[1]
        if target in redundant_lsvs:
            redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
        else:
            redundant_lsvs[target] = source
        if source in redundant_lsvs:
            redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
        else:
            redundant_lsvs[source] = target 

    # Store remaining LSVs as having no redundant pairing
    for lsv in lsv_dict.keys():
        if lsv not in redundant_lsvs:
            redundant_lsvs[lsv] = 'NaN'    

    return lsv_dict, redundant_lsvs       

def select_junctions(lsv_dict):

    psi_dict = {}
    junction_indices = {}

    samples_index = [sample for samples in lsv_dict.values() for sample in samples.keys()][0]

    # Get numpy arrays of sample PSI values for every junction set of every LSV
    for lsv,samples in lsv_dict.items():
        # Create dictionary to store sample PSIs for each junction of LSV
        psi_lists = {j:[] for j in range(0,len(lsv_dict[lsv][samples_index][3].split(';')))}
        for sample,rec in samples.items():
            junctionPSIs = rec[3].split(';')
            for j in range(len(junctionPSIs)):
                psi_lists[j] = psi_lists[j] + [float(junctionPSIs[j]),]
        # Determine which junction has greatest variance across samples
        psi_arrays = []
        for j in sorted(psi_lists.keys()):
            psi_arrays.append(np.array(psi_lists[j]))
        var_array = np.array([np.var(a) for a in psi_arrays])
        junc = np.argmax(var_array) # Index for which junction to use
        # Put LSV and sample PSIs in dictionary (dictionary of dictionaries --> {LSV:{sample:PSI}}
        # Store junction index with lsv as LSV_ID
        junction_indices[lsv] = str(junc + 1)
        #lsv_id = lsv + ':' + str(junc)
        psi_dict[lsv] = {}
        for sample,rec in samples.items():
            junctionPSIs = rec[3].split(';')
            psi_dict[lsv][sample] = float(junctionPSIs[junc])
            

    n = 0
    for lsv,samples in psi_dict.items():
        if n < 10:
            print(lsv)
            m = 0
            for s,psi in samples.items():
                if m < 10:
                    print('   ', s, psi)
                m = m + 1
            print() 
        n = n + 1

    return psi_dict, junction_indices

def output_psi_matrix(psi_dict, junc_indices, redundant_lsvs, out_dir, group_name):

    lsvs = sorted([lsv for lsv in psi_dict.keys()])
    samples = sorted([sample for sample in psi_dict[lsvs[0]].keys()])

    fn = out_dir + '/psi_matrix.' + group_name + '.csv'

    # Create CSV file for holding PSI expression matrix of remaining LSVs
    with open(fn, 'w') as oh:
        oh.write(','+','.join(lsvs)+'\n')
        for sample in samples:
            oh.write(sample + ',' + ','.join(str(psi_dict[lsv][sample]) for lsv in lsvs) + '\n')            

    # write file to keep track of which junction was selected to represent LSV
    # TODO: Add redundant LSV pair information and annotations
    fn = out_dir + '/junction_indices.' + group_name + '.csv'
    with open(fn, 'w') as oh:
        for lsv,j in junc_indices.items():
            oh.write(lsv + ',' + str(j) + '\n')

if __name__ == "__main__":

    # Program Arguments
    parser = argparse.ArgumentParser(description='Remove redundant LSVs for better network analysis')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--tsv_files', type=str, required = True, help='File containing file names (and paths) of TSV files from voila tsv.')
    requiredName.add_argument('--sample_ratio', type=float, required = True, help='Minimum ratio of samples in which LSV must exist (from majiq psi step).')
    requiredName.add_argument('--out_dir', type=str, required = True, help='Name of directory (and path) to write output files.')
    requiredName.add_argument('--group_name', type=str, required = True, help='Name of group (and version) for which LSVs were derived. (Will be file prefix)')
    args = parser.parse_args()

    # Get arguments
    fn = args.tsv_files
    sample_ratio = args.sample_ratio
    out_dir = args.out_dir
    group_name = args.group_name

    # Change to loading in LSVs from list of samples
    lsv_dict, sample_names = load_sample_data(fn)

    # Find common set of LSVs across all samples
    lsv_dict = get_shared_LSVs(lsv_dict, sample_names, sample_ratio)

    # Use arbitrary sample ID for selecting LSV information
    test_sample = [sample for samples in lsv_dict.values() for sample in samples.keys()][0] 
    test_data = [lsv_dict[lsv][test_sample] for lsv in lsv_dict.keys() if test_sample in lsv_dict[lsv]]

    # Find and remove redundant pairs of LSVs and create final LSV dictionary
    sources, targets = get_lsv_pairs(test_data)
    redundant_pairs, mutually_exclusive = find_redundant_pairs(sources, targets)
    # Filter out redundant LSVs
    lsv_dict, redundant_lsvs = filter_LSV_pairs(lsv_dict, redundant_pairs, mutually_exclusive)

    # Select which junction will be used for representing LSV in network (greatest variance across samples)
    psi_dict, junction_indices = select_junctions(lsv_dict)

    # Create expression matrix of LSVs (Rows = Samples, Columns = LSVs) 
    output_psi_matrix(psi_dict, junction_indices, redundant_lsvs, out_dir, group_name)

    
