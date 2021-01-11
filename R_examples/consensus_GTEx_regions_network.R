#################################################################################
#
#   De novo network inference of co-splicing across human tissues
#
#     This script infers co-splicing networks for tissue using SVRS and
#      identifies a consensus module set for analysis.
#     
#   --> Run this script on cluster using multi-CPUs (currently set at 20)
#
#################################################################################


# Set file paths
mainPath <- '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/'
dataPath <- paste(mainPath, 'splicing_analysis/gtex_splicing/lsv_data_v3/', sep = '')

# Load project library
source(paste(mainPath, 'project_Rscripts/shared_Rscripts/project_helper_functions.R',sep = ''))

# Load color scheme
source(paste(mainPath,'project_Rscripts/gtex_Rscripts/gtex_color_scheme.R', sep = ''))



# Load SVRs 
load(paste(dataPath, 'gtexSVRS.all.v3.n1621.RData',sep = ''))

# All tissues LSV dictionary
gtex_lsv_dict <- data.frame(read.csv(paste(dataPath, 'lsv_data_dictionary.gtex.v3.n1621.all.csv',sep=''), 
                                     header = TRUE, row.names = 1), stringsAsFactors = FALSE)

# Now load gtex_sample_summary (version 3) --> Formatted in load_GTEx_data.v3.R script
load(paste(dataPath, 'gtex_sample_summary.v3.n1621.RData',sep = ''))



############################################################
#
#  Co-splicing network inference of human tissues
#
############################################################


# Parameters from soft-thresholding analysis: B = 6, signed, bicor
B = 6
corType = 'bicor'
netType = 'signed'

# Put tissue data into multiset expr format for WGCNA
multiExpr.gtex = list()
# Initialize an appropriate array to hold the adjacencies
tissue_adjacencies = list()
# Initialize an appropriate array to hold the TOMs
tissue_TOMs = list()

# Calculate adjacencies in each individual tissue
regions <- unique(gtex_sample_summary[order(gtex_sample_summary$Tissue_Type,gtex_sample_summary$Region),]$Region)

for (tissue in regions) {
  
  # Select sample SVRs from current tissue and infer network
  tissue_samples <- rownames(gtex_sample_summary[gtex_sample_summary$Region == tissue,])
  svr.varPlot.tissue <- gtexSVRs[rownames(gtexSVRs) %in% tissue_samples,]
  tissue_adjacencies[[tissue]] = adjacency(svr.varPlot.tissue, power = B, corFnc = corType, type = netType)
  
  # Store expression matrix in multiExpr.gtex format
  multiExpr.gtex[[tissue]] = list(data = svr.varPlot.tissue)
  # Calculate TOM for current tissue
  tissue_TOMs[[tissue]] = TOMsimilarity(tissue_adjacencies[[tissue]])
} 



##########################################################################################
#
#    Calculate consensus topological overlap matrix, cluster consensus network, 
#      and identify consensus splicing module set
#
##########################################################################################

consensusTOM.gtex = pmin(tissue_TOMs[[1]], 
                         tissue_TOMs[[2]], 
                         tissue_TOMs[[3]], 
                         tissue_TOMs[[4]],
                         tissue_TOMs[[5]], 
                         tissue_TOMs[[6]],
                         tissue_TOMs[[7]], 
                         tissue_TOMs[[8]], 
                         tissue_TOMs[[9]], 
                         tissue_TOMs[[10]]
)

# Set row and column names to new consensus TOM
rownames(consensusTOM.gtex) <- rownames(tissue_adjacencies[[1]])
colnames(consensusTOM.gtex) <- colnames(tissue_adjacencies[[1]])

####################################################################
#
#   Use hierarchical clustering with dynamic tree cutting to
#     determine initial module count for spectral clustering
#
####################################################################

# Clustering of SVR TOM
ds = 2
minSize = 60 
dissConsTOM.gtex <- 1-consensusTOM.gtex
conSVRTree.gtex = hclust(as.dist(dissConsTOM.gtex), method = "average")
hcColors = labels2colors(cutreeDynamic(dendro = conSVRTree.gtex, distM = dissConsTOM.gtex,
                                             deepSplit = ds, minClusterSize = minSize,
                                             pamRespectsDendro = FALSE))
table(hcColors)
length(table(hcColors))

# Merge close modules from hierarchical clustering and determine K
mergeThres <- 0.25
mergedHCColors <- merge_consensus_modules(hcColors, multiExpr.gtex, mergeThres, FALSE)
table(mergedHCColors)
length(table(mergedHCColors))


# Now cluster using spectral clustering
tau = 0
spherical = TRUE
k <- length(table(mergedHCColors))

# Run normalized spectral clustering on TOM matrix
registerDoParallel(cores=20)
spColors <- iterate_spectral_clustering(consensusTOM.gtex, k = k, tau = 0, nstart = 100, 
                                        iterations = 101, numcores = 20)



# Merge correlated modules
gtexMods <- merge_consensus_modules(spColors, multiExpr.gtex, 0.25, TRUE)
table(gtexMods)
length(table(gtexMods))
min(table(gtexMods))
max(table(gtexMods))



######################################################################
#
#   Calculate consensus splicing module eigenvectors and 
#    consensus splicing module clustering (dendrogram)
#
######################################################################

# Contains a list with length of number of tissue regions
# conMods.multiSet.gtex$data contains the sample X eigengene expression matrix for each tissue
conMods.multiSet.gtex = multiSetMEs(multiExpr.gtex, universalColors = gtexMods, excludeGrey = TRUE)
# Create list of sample X module expression matrices (one for each tissue)
#  Names are now colors (no "ME" substring)
conMods.gtex <- multi_to_colors(conMods.multiSet.gtex)

# Calculate consensus dissimilarity of consensus modules to cluster and find hierarchical structure
conDiss.gtex = consensusMEDissimilarity(conMods.multiSet.gtex)
# Cluster consensus modules and create consensus splicing module tree
conTree.gtex = hclust(as.dist(conDiss.gtex), method = "average")


#######################################################################################
#
#     Create set-specific and across tissue module expression dataframes
#
#######################################################################################

# Merge consensus module expression into single dataframe
#  Splicing module values are set-specific 
modExpr.setSpecific.gtex <- multi_to_single_expr(conMods.gtex)

# Create mondule expression dataframe with splicing module
#  values determined across all samples simultaneously
modExpr.gtex = moduleEigengenes(gtexSVRs, gtexMods, excludeGrey = TRUE)$eigengenes
#modExpr.gtex = moduleEigengenes(svr.varPlot.allTissues, gtexMods, excludeGrey = TRUE)$eigengenes
names(modExpr.gtex) <- str_to_title(substring(names(modExpr.gtex), 3))


# Make tree labels the same as new mod labels
conTree.gtex$labels <- colnames(modExpr.gtex)



####################################################################
#
#   Get node connectivities and module memberhsip for each tissue
#
####################################################################

node_stats.gtexSpliceNet <- get_node_stats(tissue_adjacencies, gtexMods, multiExpr.gtex, isList = TRUE, 
                                           convert_SVRs = TRUE, lsv_dict = gtex_lsv_dict, getEntrez = TRUE)


########################################################
#
#   Save all data to file for module analysis
#
#######################################################
save(gtexSVRs, gtexMods, gtex_sample_summary, tissue_adjacencies, multiExpr.gtex, tissue_TOMs, conMods.multiSet.gtex,
     consensusTOM.gtex, conTree.gtex, modExpr.gtex, modExpr.setSpecific.gtex, conMods.gtex, node_stats.gtexSpliceNet, gtex_lsv_dict,
     file = paste(mainPath, 'rData_files/gtex_rData/consensus_network_files/gtexConsData.RData', sep = ''))


# Save results of hierarhical clustering (and merging) as well
save(mergedHCColors, file = paste(mainPath, 'rData_files/gtex_rData/consensus_network_files/hClustModules.gtex.RData', sep = ''))

