#############################################################################################
#
#   Beat AML Junction Read and LSV Coverage EDA 
#
#    Description:
#         
#    Plotting script to evaluate coverage of AML samples and decide on threshold 
#     for removing low coverage samples prior to network analysis
#
#    - The sample summary includes 411 unique patient samples and 32 control samples
#    - A single patient sample will dropped due to low coverage prior to analysis
#
#############################################################################################

# Set file paths
mainPath <- '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/'
dataPath <- paste(mainPath, 'splicing_analysis/beatAML_splicing/lsv_data_v5/', sep = '')

# Load project library
source(paste(mainPath, 'project_Rscripts/shared_Rscripts/project_helper_functions.R',sep = ''))


####################################
#
#   LOAD SAMPLE DATA
#
####################################


# Load the initial expession matrix containing all LSVs (including redundant and missing in samples)
#  --> amlExpr.initial  (dataframe)
load(paste(dataPath, 'psi_matrix.x442.v5.all.RData', sep = ''))

# Load LSV read counts for intial set (can use for final set too)
# --> lsv_read_counts.aml (dataframe)
load(paste(dataPath, 'lsv_read_counts.aml.x442.RData', sep = ''))

# --> aml_sample_summary_all (dataframe)
load(paste(dataPath, 'aml_sample_summary_all.x442.v5.RData', sep = ''))



# Load AML annotation colors for plotting (pre-filtering version)
#   --> Must have aml_sample_summary_all dataframe loaded
source(paste(mainPath, 'project_Rscripts/beatAML_Rscripts/initial_dataset_scripts/aml_color_scheme_preFiltering.R', sep = ''))



#################################################################
#
#  Plot sample counts and quantifiable LSVs by disease group
#
#################################################################

# Sample counts for each disease group
disease_type_counts <- ggplot(aml_sample_summary_all, aes(x = fct_infreq(DiseaseType), fill = DiseaseType)) +
  geom_bar(stat = 'count') + scale_y_continuous(breaks=seq(0,225,25)) + 
  scale_fill_manual(values = disease_colors) +
  geom_text(stat='count', aes(label=..count..), position=position_stack(0.5)) + 
  labs(title=paste("Total Samples (N = ", as.character(dim(aml_sample_summary_all)[1]), ")", sep = ''),
       x ="Disease Type", y = "Count") 
disease_type_counts



# Box plots showing quantifiable LSVs in each sample by disease group
total_lsvs <- dim(amlExpr.initial)[2] 
disease_type_coverage <- ggplot(aml_sample_summary_all, aes(x=fct_infreq(DiseaseType),y=Quantified_LSVs, fill=DiseaseType)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) + scale_fill_manual(values = disease_colors) + 
  scale_y_continuous(breaks=seq(10000,24000,2000), limits = c(10000,24000), sec.axis = sec_axis(~./total_lsvs, breaks = seq(0.4, 1.0, 0.1))) + 
  labs(title=paste("Quantified LSVs Per Sample (", as.character(total_lsvs), " Total LSVs)", sep = ''),
       x ="Disease Type", y = "Count") +
  geom_hline(aes(yintercept = total_lsvs*.9), color = "red", size = .3) 
disease_type_coverage


ggarrange(disease_type_counts, disease_type_coverage, nrow = 1, common.legend = TRUE, legend = "bottom")




###########################################################################################
#
#    Look at junction read coverage and quantified LSV coverage across AML samples
#
###########################################################################################

# Order samples by quantifiable LSVs (highest to lowest)
summary.toPlot <- aml_sample_summary_all[order(-aml_sample_summary_all$coverage_ratio),]

# Match PSI matrix and LSV read count matrix
expr.toPlot <- amlExpr.initial[rownames(summary.toPlot),]
lsv_read_counts.toPlot <- lsv_read_counts.aml[colnames(expr.toPlot), rownames(expr.toPlot)]

# Create color annotations for disease type, specimen type, and sequencing site
top_annot <- HeatmapAnnotation(DiseaseType = summary.toPlot$DiseaseType,
                               TissueType = summary.toPlot$SpecimenType,
                               Site = summary.toPlot$SpecimenSite,
                               col = list(DiseaseType = disease_colors,
                                          TissueType = tissue_colors,
                                          Site = site_colors),
                               annotation_legend_param = list(DiseaseType = list(ncol = 2),
                                                              TissueType = list(ncol = 2),
                                                              Site = list(nrow = 4)))


# Create plot object
covPlot <- coverage_plot(expr.toPlot, summary.toPlot, top_annot, lsv_read_counts.toPlot, 
                         reads_axis = c(seq(0,250000000,50000000)), lsv_axis = c(seq(0,2000,500)), drop_breaks = 4000)

# Create tilte
ht_title <- paste('Beat AML LSV Coverage (N = ', as.character(dim(summary.toPlot)[1]), '); ', 
                  as.character(dim(expr.toPlot)[2]), ' Total LSVs', sep = '')

# Draw heatmap object
draw(covPlot$ht, annotation_legend_list = covPlot$legend_list, annotation_legend_side = "bottom", 
     heatmap_legend_side = "bottom", column_title = ht_title)

# Save to file (difficult to automate figure with appropriate dimensions; easier to view within Rstudio)
#png(file="/home/users/eggerj/tmp_transfer/lsvCoverage.aml.png", width = 28, height = 20, res=600, units = 'in')
#draw(covPlot$ht, annotation_legend_list = covPlot$legend_list, annotation_legend_side = "bottom", 
#     heatmap_legend_side = "bottom", column_title = ht_title)
#dev.off()


############################################################
#
#  Subset to only look at lowest coverage samples
#
#############################################################
covPlot <- coverage_plot(expr.toPlot, summary.toPlot, top_annot, lsv_read_counts.toPlot, c(seq(0,250000000,50000000)), c(seq(0,2000,500)), 
                         drop_breaks = 4000, drop_height = unit(2.5, "cm"), quant_height = unit(1.5, "cm"), read_coverage_height = unit(2.5, "cm"), 
                         lsv_coverage_height = unit(1.5, "cm"))
n <- 20
e <- dim(covPlot$ht)[2]
s <- e - n - 1

draw(covPlot$ht[,s:e], annotation_legend_list = covPlot$legend_list, annotation_legend_side = "bottom", heatmap_legend_side = "bottom", 
     column_title = paste('Lowest', as.character(n), 'Samples (Ratio of Quantified LSVs)', sep = ' '))



############################################################################################
#
#    Now remove low coverage samples and save new sample summary analysis set to file
#
############################################################################################

# Setting sample LSV coverage minimum to 90% (removes single sample)
thres <- 0.9
aml_sample_summary <- aml_sample_summary_all[aml_sample_summary_all$coverage_ratio >= thres, ] 

# Save new sample summary to file
#save(aml_sample_summary, file = '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/beatAML_splicing/lsv_data_final/aml_sample_summary.x441.v5.RData')
