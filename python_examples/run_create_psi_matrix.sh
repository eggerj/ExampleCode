#!/bin/sh
#SBATCH --job-name=create_psi_matrix
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH -A zheng_lab
#SBATCH -p exacloud
#SBATCH --error=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/slurm_scripts/matrix_log_files/create_psi_matrix_v5.%J.err
#SBATCH --output=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/slurm_scripts/matrix_log_files/create_psi_matrix_v5.%J.out

######################################################################################
#
#  Slurm job to run final Python processing script for collecting and annotation
#   LSV quantities across samples. Needs paths to sample TSV files as well as
#   exon and junction coordinates (produced manually).
#
######################################################################################


# Set variables
DIR=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis

TSV_FILES=$DIR/beatAML_splicing/voila_tsv_v5/tsv_files.filtered.n441.v5.txt  # Paths to sample TSV files (post-filtering)
OUT_DIR=$DIR/beatAML_splicing/lsv_data_final
PY_SCRIPT=$DIR/py_scripts/create_psi_matrix.py

# *** Need to make junction coordinates before hand from TSV files (UNIX command line functions)
JUNC_COORDS=$DIR/beatAML_splicing/voila_tsv_v5/junction_coords.v5.txt

# *** Need to make exon coordinates before hand from non-simplified TSV files (UNIX command line functions)
EXON_COORDS=$DIR/beatAML_splicing/non_simplified_voila_tsv_v5/exon_coordinates.v5.txt

# GTF file used for splicing event annotations
GTF=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/gff_files/hs.grch37.75.exons.gtf

# Create identifier label for dataset
GROUP_NAME=x441.v5.all



# Activate majiq virtual environment
source /home/exacloud/lustre1/zheng_lab/users/eggerj/majiq3.6/bin/activate

# Now run using work node on Exacloud
srun -N 1 -n 1 -c 1 python $PY_SCRIPT --tsv_files $TSV_FILES --sample_ratio 0.0 --out_dir $OUT_DIR --group_name $GROUP_NAME --junc_coords $JUNC_COORDS --exon_coords $EXON_COORDS --gtf $GTF --redundant keep


