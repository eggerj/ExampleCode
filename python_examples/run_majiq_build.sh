#!/bin/sh
#SBATCH --job-name=majiq_build
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -p exacloud 
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=500M
#SBATCH -A zheng_lab
#SBATCH --exclude=exanode-6-16,exanode-0-1,exanode-0-7
#SBATCH --error=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/slurm_scripts/log_files/majiq_build_v5.%J.err
#SBATCH --output=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/slurm_scripts/log_files/majiq_build_v5.%J.out

######################################################################################
#
#   Slurm script for constructing alternative splicing graphs
#    using MAJIQ build function.
#
#   Build needs to be performed on all samples simultaneously, 
#    however each bam file needs to be sorted first.
#   This script first sorts bam files in parallel and copies them to single
#    work node for build step. Due to failing work nodes, a second
#    step is performed to ensure all files are sorted and copies.
#
#######################################################################################


# Set variables
DIR=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material
BAM_FILES=$DIR/beatAML_files/rnaSeq_aml_plus_control_bamPaths.txt
SCRATCH_TMP=/mnt/scratch/majiq_build
MAJIQ_OUT=$DIR/splicing_analysis/beatAML_splicing/majiq_build_v5
BAM_TMP=$DIR/splicing_analysis/beatAML_splicing/bam_tmp_v5
GFF=$DIR/gff_files/Homo_sapiens.GRCh37.75.gff3
#CONFIG=$DIR/beatAML_files/settings.beatAML.ini
CONFIG=$DIR/beatAML_files/settings.aml_plus_control.ini

# Try getting node list
nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST) # Getting the node names
nodes_array=( $nodes )
node1=${nodes_array[0]}

# Create temporary directory in requested scratch space to store files
srun -N 1 -n 1 -c 1 -w $node1 mkdir $SCRATCH_TMP
srun -N 1 -n 1 -c 1 -w $node1 mkdir $SCRATCH_TMP/majiq_out

# Create directory in Lustre for storing subsets of files during parallel sorting
srun -N 1 -n 1 -c 1 mkdir $BAM_TMP

##################################################################################################################
#
# PART 1: Sort bam files (in parallel) and store in scratch space 
#         (wait until all are finished before part 2)
#
##################################################################################################################

# Create an array for storing bam files
readarray -t SAMPLES < $BAM_FILES

# Parse array of bam files and sort/index in parallel using sambamba 
# Reserve first node for moving files, use rest for sorting/indexing
n=0
m=$(($SLURM_JOB_NUM_NODES - 1))
for F in ${SAMPLES[@]}; do
    if [[ $n -lt m ]]
    then
        n=$(( $n + 1))
    else
        n=1
    fi
    nodeN=${nodes_array[n]}
    fn="$(rev <<< "$F" | cut -d'/' -f 1 | rev)" # Remove path for writing output to new directory
    echo $fn
    { srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK -w $nodeN sambamba-0.7.1-linux-static sort -t $SLURM_CPUS_PER_TASK -m 4000M --tmpdir $BAM_TMP -o $BAM_TMP/$fn $F ;
      srun -N 1 -n 1 -c 1 -w $node1 mv $BAM_TMP/$fn* $SCRATCH_TMP/ ; } &
done
wait 

###############################################################################################
#
# Fail Safe: Make sure all files actually made it and rerun sorting for any files that din't
#
###############################################################################################
filecount=$(srun -N 1 -n 1 -c 1 -w $node1 find $SCRATCH_TMP -type f -name "*.bam" -printf x | wc -c)
numfiles=${#SAMPLES[@]}
while [[ $filecount != $numfiles ]]; do
    for F in ${SAMPLES[@]}; do
        fn="$(rev <<< "$F" | cut -d'/' -f 1 | rev)"
        { srun -N 1 -n 1 -c 1 -w $node1 test -f $SCRATCH_TMP/$fn || { echo "$fn does not exist" ;
          srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK sambamba-0.7.1-linux-static sort -t $SLURM_CPUS_PER_TASK -m 4000M --tmpdir $BAM_TMP -o $BAM_TMP/$fn $F ;
          srun -N 1 -n 1 -c 1 -w $node1 mv $BAM_TMP/$fn* $SCRATCH_TMP/ ; } } &
    done
    wait
    filecount=$(srun -N 1 -n 1 -c 1 -w $node1 find $SCRATCH_TMP -type f -name "*.bam" -printf x | wc -c)
done


# Check files actually made it before running MAJIQ
srun -N 1 -n 1 -c 1 -w $node1 ls -lh $SCRATCH_TMP

# Remove temp folder on lustre
srun -N 1 -n 1 -c 1 rm -r $BAM_TMP


##################################################################################################################
#
# PART 2: Run MAJIQ build (single task) using all bam files (after all have been indexed) 
#
##################################################################################################################

# Activate majiq virtual environment
source /home/exacloud/lustre1/zheng_lab/users/eggerj/majiq3.6/bin/activate

# Run MAJIQ build using all bam files (.ini file indicates that bam files are in temp directory)
srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK -w $node1 majiq build $GFF -c $CONFIG -j $SLURM_CPUS_PER_TASK -o $SCRATCH_TMP/majiq_out \
                                                             --min-experiments .08 \
                                                             --min-intronic-cov 1.0 --minreads 12 --minpos 6

# Move majiq output files from tmp directory to output directory on Lustre and remove 
srun -N 1 -n 1 -c 1 -w $node1 cp $SCRATCH_TMP/majiq_out/* $MAJIQ_OUT/
srun -N 1 -n 1 -c 1 -w $node1 rm -r $SCRATCH_TMP
