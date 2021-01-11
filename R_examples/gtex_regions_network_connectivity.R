##########################################################################
#
#  Soft-thresholding Beta Value Selection
#
#     This script first formulates SVRs from LSVs of GTEx tissues
#     and then evaluates scale-free fit of each tissue network after 
#     varying soft-thresholding values. Soft-thresholding used to
#     minimize effect of spurious correlations during network inference.
#
############################################################################

# Set file paths
mainPath <- '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/'
dataPath <- paste(mainPath, 'splicing_analysis/gtex_splicing/lsv_data_v3/', sep = '')

# Load project library
source(paste(mainPath, 'project_Rscripts/shared_Rscripts/project_helper_functions.R',sep = ''))

# Load color scheme
source(paste(mainPath,'project_Rscripts/gtex_Rscripts/gtex_color_scheme.R', sep = ''))



# Load expression data and lsv dictionary
load(paste(dataPath, 'psi_matrix.gtex.v3.n1621.all.RData',sep = ''))

# Only select LSVs quantified in all samples
gtexExpr.all <- gtexExpr[,colSums(is.na(gtexExpr)) == 0]
dim(gtexExpr.all)

# All tissues LSV dictionary
gtex_lsv_dict <- data.frame(read.csv(paste(dataPath, 'lsv_data_dictionary.gtex.v3.n1621.all.csv',sep=''), 
                                     header = TRUE, row.names = 1), stringsAsFactors = FALSE)

# Now load gtex_sample_summary (version 3) --> Formatted in load_GTEx_data.v3.R script
load(paste(dataPath, 'gtex_sample_summary.v3.n1621.RData',sep = ''))



#########################################################################
#
# Create splice variant regions from PSI matrix and LSV dictionary
#
##########################################################################

#gtexSVRs <- make_SVRs(gtexExpr.all, gtex_lsv_dict)
# Save for later to avoid recomputing 
#save(paste(dataPath, 'gtexSVRS.all.v3.n1621.RData',sep = ''))
# Load SVRs already made to save time
load(paste(dataPath, 'gtexSVRS.all.v3.n1621.RData',sep = ''))


#######################################################################
#
#   Plot scale-free fit across beta values for each network
#
#######################################################################

# Parameter full names
corNames <- list('cor' = 'Pearson Correlations', 'bicor' = 'Biweight Midcorrelations')    

# Select network type and R-squared cutoff minimum
corType <- 'bicor'    
netType <- 'Signed'
rCut <- 0.905
group_name <- 'GTEx Regions'

# Parse each tissue, create tissue-specific SVR matrix, and calculate scale-free fit and connectivity metrics  
# Store results of connectivity metrics here   
powerTables = vector(mode = "list")
minSF = vector(mode = "list")
powers = c(1:10)
for (region in sort(unique(gtex_sample_summary$Region))) {
  
  # Select samples of current region
  print(region)
  region_samples <- rownames(gtex_sample_summary[gtex_sample_summary$Region == region,])
  svr.region <- gtexSVRs[rownames(gtexSVRs) %in% region_samples,]
  
  # Compute SF fit across beta values
  sft = pickSoftThreshold(svr.region, powerVector = powers, verbose = 5, networkType = tolower(netType), 
                          corFnc = corType, RsquaredCut = rCut)
  
  # Store result for each region type
  powerTables[[region]] = sft$fitIndices[,c(1,2,3,5)]
  # Store beta value for which threshold was reached
  # If threshold was not reached, select beta as value for highest scale-free fit
  if (is.na(sft$powerEstimate)) {
    minSF[[region]] <- powerTables[[region]][which.max(powerTables[[region]]$SFT.R.sq),]$Power
  } else {  
    minSF[[region]] = sft$powerEstimate
  }
}  

# Format dataframes for input into network connectivity plotting function
powerDF <- bind_rows(powerTables, .id = "Region")
minDF <- merge(data.frame("Region" = names(minSF), 
                          "Power" = unlist(minSF,use.names=F)),
               powerDF,by=c("Region","Power"),all.x = T)

# Create connectivity network for current correlation type
if (corType == "cor") {
  conn_plot <- multi_network_connectivity_plots(rCut, powerDF, minDF, corNames[[corType]], netType, 
                                                group_name, "bottom", dim(gtexSVRs)[2]) 
} else {
  conn_plot <- multi_network_connectivity_plots(rCut, powerDF, minDF, corNames[[corType]], netType, 
                                                group_name, "bottom", dim(gtexSVRs)[2]) 
}


conn_plot

# Save to file
#fn <- paste('sfPlot.gtex', corType, tolower(netType), sep ='.')
#ggsave(filename = paste('/home/users/eggerj/tmp_transfer/', fn, '.png', sep = ''), width = 18, height = 9, dpi = 300)
#dev.off()