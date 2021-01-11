#!/bin/sh
#SBATCH --job-name=aml_splicing_lasso
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-122
#SBATCH -A zheng_lab
#SBATCH -p exacloud
#SBATCH --time=2:00:00
#SBATCH --exclude=exanode-6-24
#SBATCH --error=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/beatAML_Rscripts/beatAML_slurms_4_R/log_files/aml_splicing_lasso.%A_%a.err
#SBATCH --output=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/beatAML_Rscripts/beatAML_slurms_4_R/log_files/aml_splicing_lasso.%A_%a.out

##################################################################
#
#  Slurm array job for modeling drug response with LASSO
#
#   Call preprocess script for each drug of array to load
#    data and run LASSO modeling
#
##################################################################


# Load in array of drug names
PARAM_SETS=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/beatAML_Rscripts/beatAML_drug_sets/all_drugs.txt
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $PARAM_SETS)

# Run script for current drug of array
SCRIPT=/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/beatAML_Rscripts/preprocess_lasso_splicing.R
srun Rscript --vanilla $SCRIPT $LINE
