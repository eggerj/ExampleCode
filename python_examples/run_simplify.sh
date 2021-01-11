#!/bin/sh
#SBATCH --job-name=simplify
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:30:00
#SBATCH --mem-per-cpu=6G
#SBATCH -A zheng_lab
#SBATCH -p exacloud
#SBATCH --exclude=exanode-1-44,exanode-1-34,exanode-4-3,exanode-6-44,exanode-6-47,exanode-4-8,exanode-4-9,exanode-4-10,exanode-0-2
#SBATCH --error=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/slurm_scripts/log_files/simplify_v5.%J.err
#SBATCH --output=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/slurm_scripts/log_files/simplify_v5.%J.out

########################################################################################################
#
#  Slurm script to run pair of Python scripts for simplifying splicing graphs produced by MAJIQ
#
########################################################################################################

# Set variables
DIR=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis

# File of majiq file names with full paths
MAJIQ_FILES=$DIR/beatAML_splicing/majiq_build_v5/majiq_files.v5.txt  # <-------- NEEDS TO BE CREATED MANUALLY FIRST BEFORE RUNNING

# Other variables
OUT_DIR=$DIR/beatAML_splicing/majiq_simplify_v5
OUT_FILE=$OUT_DIR/majiq_simplified_files_v5.txt
PY_SCRIPT1=$DIR/py_scripts/optimize_simplify_majiq_files.py
PY_SCRIPT2=$DIR/py_scripts/write_new_majiq_files.py
TMP=$OUT_DIR/tmp

mkdir $TMP

# Activate majiq virtual environment
source /home/exacloud/lustre1/zheng_lab/users/eggerj/majiq3.6/bin/activate

srun -N 1 -n 1 -c 1 python $PY_SCRIPT1 --majiq_files $MAJIQ_FILES --junction_min_reads 10 --junction_min_experiments .08 --lsv_min_reads 20 --lsv_min_experiments .16 \
                                       --single_sample_junction_min_reads 50 --out_dir $OUT_DIR

# Create an array for MAJIQ files
readarray -t SAMPLES < $MAJIQ_FILES

# Parse through majiq files and run script to create new majiq files
for F in ${SAMPLES[@]}; do
    echo $F
    srun -N 1 -n 1 -c 1 python $PY_SCRIPT2 --majiq_file $F --out_dir $OUT_DIR &
done
wait


# Change file extensions from .npz to .majiq
srun -N 1 -n 1 -c 1 rename npz majiq $OUT_DIR/*.npz

# Make file listing for majiq psi
srun -N 1 -n 1 -c 1 ls -1 $OUT_DIR/*.majiq > $OUT_FILE
