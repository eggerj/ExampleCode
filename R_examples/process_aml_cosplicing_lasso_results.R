####################################################################################
#
#   Process LASSO Drug Response Modeling Results
#
#     Script to load in results of lasso regression models for each drug
#
####################################################################################

# Set file paths
mainPath <- '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/'

# Load project library
source(paste(mainPath, 'project_Rscripts/shared_Rscripts/project_helper_functions.R',sep = ''))

# Load LASSO library 
source(paste(mainPath, 'project_Rscripts/beatAML_Rscripts/lasso_drug_response.R', sep = ''))

# Load custom helper script for creating color palettes for sample annotations
source(paste(mainPath, 'project_Rscripts/beatAML_Rscripts/aml_color_scheme.R',sep = ''))

# Load library for analyzing results of LASSO drug response modeling
source(paste(mainPath, 'project_Rscripts/shared_Rscripts/buzzsaw_plotting_functions.R',sep = ''))


# Load module set
load(paste(mainPath, 'rData_files/beatAML_rData/aml_network_files/amlNetworkData.RData',sep = ''))
modules.aml <- unique(str_to_title(amlMods))

# Load remaining sample data
load(paste(mainPath, 'beatAML_files/BeatAML_github_data/preproc_data.RData',sep = ''))


# Load in names of all drugs used in modeling (full names)
drugNames <- as.vector(read.csv(file = paste(mainPath, 'project_Rscripts/beatAML_Rscripts/beatAML_drug_sets/all_drugs.txt', sep = ''),
                                header = FALSE, stringsAsFactors = FALSE)[,1])


# Create named vector for retreiving full name when input is just the first 
#  name of drug (short name required for Slurm array)
drugDict <- drugNames
names(drugDict) <- vapply(strsplit(drugNames," "), `[`, 1, FUN.VALUE=character(1))



#########################################################
#
#     PLOT NUMBER OF SAMPLES USED FOR EACH DRUG
#
#########################################################

#  Need to convert sample names again
modExpr.aml.originalID <- modExpr.aml
rownames(modExpr.aml.originalID) <- sample_summary.toPlot[rownames(modExpr.aml),]$Original_LabID
modExpr.aml.originalID <- as.matrix(modExpr.aml.originalID)

# Now plot 
plot_sample_sizes(split.drugs, clin.dt, modExpr.aml.originalID, wes.clin)






##########################################################
#
#  LOAD AND MERGE ALL LASSO DRUG RESULTS
#
##########################################################

# ~~~~ CURRENT SET OF FOCUS ~~~~~ #
# Splicing Modules & Mutations (version 3: min mutation = .03 * N)

# Lasso results folder for splicing modules and mutations
inPath <- paste(mainPath,'rData_files/beatAML_rData/lasso_data/splicing_mutations.v3/', sep = '')

# Now load all data and store
#pred_data <- load_data(inPath, drugDict)

# Save and use again later
#save(pred_data, file = paste(mainPath,'rData_files/beatAML_rData/lasso_data/splicing_mutations.v3/pred_data.splicing.mut.RData',sep = ''))
load(paste(mainPath,'rData_files/beatAML_rData/lasso_data/splicing_mutations.v3/pred_data.splicing.mut.RData', sep =''))

# Get list of all features (modules and mutations)
all_features <- get_features(coef.list.all)

# Store in variables
drugRsquaredVals <- pred_data$r2
coef.list.all <- pred_data$coefs




############################################
#
#   PREDICTION RESULTS ON TEST SET
#
############################################

# Plot distribution of R^2 values from 100 model iterations for each drug
plot_rSquared(drugRsquaredVals)




######################################################################################
#
#   See if any mutated genes used are found in co-splicing modules
#
#######################################################################################

mutation_to_module(node_stats.amlSpliceNet, all_features)


###########################################################################################
#
#   LOOK AT FEATURE CORRELATIONS ACROSS ALL LASSO MODELS
#
############################################################################################


# Get non-zero rate for each feature and filter drugs not having at least one coefficient ratio greater than 25%
drugCoefficients <- make_coefficient_matrix(coef.list.all, all_features, thres = 0.25)


# Filter drugCoefficients for modules with Rsquared >= 0.0
medianR2 <- aggregate(drugRsquaredVals$Rsquared, list(drugRsquaredVals$Inhibitor), median)
selectR2 <- as.character(medianR2[medianR2$x >= 0.0,1])
drugCoefficients <- drugCoefficients[selectR2,]
dim(drugCoefficients)



# Sort models by highest coefficients across splicing modules
#  *** USE THIS FOR SPLICING MODULE PLOT ****
subCoef <- data.frame(Coef = as.vector(as.matrix(drugCoefficients[,modules.aml])),
                      Inhibitor = rep(rownames(drugCoefficients),length(modules.aml)),
                      Module = rep(modules.aml, each = length(rownames(drugCoefficients))))
drug_order <- as.character(unique(subCoef[order(-subCoef$Coef),]$Inhibitor))
drugCoefficients <- drugCoefficients[drug_order,]


# Get model ranks based on high splicing module correlations (save for later)
drugRanksSplicing <- data.frame(drug_order, 1:length(drug_order))
drugRanksSplicing$drug_order <- as.character(drugRanksSplicing$drug_order)



# Subset out splicing modules
moduleDrugCoefficients <- as.matrix(drugCoefficients[,modules.aml])
# Subset out mutations
mutationDrugCoefficients <- as.matrix(drugCoefficients[,!colnames(drugCoefficients) %in% modules.aml])

# Count models with high module correlations
length(rownames(moduleDrugCoefficients)[rowSums(moduleDrugCoefficients >= 0.5, na.rm = TRUE) > 0])
# Count modules with high model correlations
length(colnames(moduleDrugCoefficients)[colSums(moduleDrugCoefficients >= 0.5, na.rm = TRUE) > 0])


## PLOT LASSO MODEL FEATURE CORRELATIONS 

select_rows <- 1:16
modPreds <- coefficient_ratio_heatmap(moduleDrugCoefficients, title = 'Splicing Modules', select_rows = select_rows, fontSize = 10)
mutPreds <- coefficient_ratio_heatmap(mutationDrugCoefficients[,-c(27:28)], title = 'Mutations', select_rows = select_rows, fontSize = 10)
ComplexHeatmap::draw(modPreds + mutPreds, heatmap_legend_side = "bottom")




###############################################################################
#
#   LOOK AT FREQUENCY OF CO-OCCURRENCES BETWEEN MODULES AND MUTATIONS
#
################################################################################


# Select drug to use for Buzzsaw co-occurence plot
drug <- 'Venetoclax'
drug <- 'Sorafenib'
drug <- 'Dasatinib'
drug <- 'Flavopiridol'
drug <- 'Foretinib (XL880)'
drug <- 'PD173955'
drug <- 'Nilotinib'
drug <- 'Ponatinib (AP24534)'


# MAKE BUZZSAW PLOT

# Use if saving to file  
#d <- names(drugDict[drugDict == drug])
#fn <- paste('/home/users/eggerj/tmp_transfer/buzzsaw.', d, '.splice.mut.png', sep = '')
#png(filename = fn, width = 3200, height = 2900, res = 240) # Need to adjust for each drug

# Create dataframe showing co-occuring ratios
df <- make_co_occurence_dataframe(coef.list.all[[drug]])

# Initialize circos and create barplots
make_circoBP(df)

# Add label for largest co-occurence
label_largest_co_occurence(df)

# Get module and mutation colors 
cols <- buzzsaw_colors(all_features)

# Add colored segment around circle for each feature
add_colored_segments(df, cols)

# Add links between barplots
add_links_seq(df, cols)

# Clear and reset
circos.clear()
circos.par(RESET = TRUE)

# For saving
#dev.off()
  




####################################################################################
#
#   Make dot plots showing LASSO coefficient ratios grouped by drug families
#
####################################################################################


# Load drug family groupings
drug_fams <- read.csv(file = paste(mainPath, 'beatAML_files/drug_family_assignments.csv', sep = ''),
                      header = TRUE, stringsAsFactors = FALSE)


# Group data into drug families and plot frequencies
plot_drug_family_frequencies(drugCoefficients, drug_fams)

