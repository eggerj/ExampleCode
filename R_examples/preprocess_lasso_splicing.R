###########################################################################################################
#
#   Pre-processing script to load sample data (splicing modules, mutation calls, drug response) 
#    for current drug and model using iterative LASSO regression (using lasso_drug_response.R script)
#
#   Script is ran for each element (drug) within a Slurm array
#    --> Initialized from run_aml_splicing_lasso.sh Slurm script
#
###########################################################################################################   

# Load LASSO iterative modeling functions
source('/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/beatAML_Rscripts/lasso_drug_response.R')


##########################################
#
#   LOAD ALL SAMPLE DATA
#
##########################################


# Get name of current drug from Slurm array
args = commandArgs(trailingOnly=TRUE)
drugName = args[1]

# Load clinical data for AML samples from Tyner et al. 2017 (already pre-processed)
load("/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/beatAML_files/BeatAML_github_data/preproc_data.RData")

# Load co-splicing module data for AML samples
load('/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/rData_files/beatAML_rData/aml_network_files/amlNetworkData.RData')


#############################################
#
#   FORMAT SAMPLE PARAMETERS FOR LASSO
#
#############################################

# RNA-seq sample IDs sometimes differ from clinical sample IDs 
#  --> Need to rename samples for including clinical data
modExpr.aml.originalID <- modExpr.aml
rownames(modExpr.aml.originalID) <- sample_summary.toPlot[rownames(modExpr.aml),]$Original_LabID
modExpr.aml.originalID <- as.matrix(modExpr.aml.originalID)


# Create named vector for retreiving full name when input is just the first name of drug (for use with Slurm array)
drugDict <- names(split.drugs)
names(drugDict) <- vapply(strsplit(names(split.drugs)," "), `[`, 1, FUN.VALUE=character(1))


# Specific folder to store LASSO model results of each drug
outPath <- "/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/rData_files/beatAML_rData/lasso_data/splicing_mutations.v3/"

# Output file name of current drug
outRez <- paste(outPath, drugName, ".RData", sep = '')

# For log file on Exacloud
print(drugName)

# Minimum positive calls for each mutation to be included during modeling
# based on samples having response data for current drug
mut_thres <- 0.03

####################################################
#
#   NOW RUN LASSO MODELING FOR CURRENT DRUG 
#
####################################################

run_combined_datatypes(split.drugs[drugDict[drugName]], rna.clin, wes.clin, 
                       modExpr.aml.originalID, var.dt, clin.dt, 
                       outFile = outRez, minMut = 0.03)





