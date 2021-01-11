#################################################################################
#
#   GTEx Splicing Network Module Analysis
#
#     Load consensus module set and characterize inter-modular variation
#       and differential splicing of modules across tissue types
# 
################################################################################

# Set file paths
mainPath <- '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/'
dataPath <- paste(mainPath, 'rData_files/gtex_rData/consensus_network_files/', sep = '')

# Load project library
source(paste(mainPath, 'project_Rscripts/shared_Rscripts/project_helper_functions.R',sep = ''))

# Load color scheme
source(paste(mainPath,'project_Rscripts/gtex_Rscripts/gtex_color_scheme.R', sep = ''))


# Load all module data (LARGE) for consensus GTEx splicing network
load(paste(dataPath, 'gtexConsData.RData', sep = ''))




###################################################################
#
#   FIRST EVALUATE QUALITY OF MODULE DETECTION USING ZSCORES
#
###################################################################

# Load from cluster (gtex_mod_quality_netrep.R)
load(paste(dataPath, 'gtex_netrep_modQuality.n10000.RData', sep = ''))

# Select appropriate module statistics
#modStatSelect <- which(colnames(gtex_zScoreMats.spectral[[tissue]]) %in% c('avg.weight', 'coherence', 'avg.cor', 'avg.contrib'))
modStatSelect <- which(colnames(gtex_zScoreMats.spectral[[1]]) %in% c('Zdensity'))


# Compute heatmaps showing Zscores from summary statistics
# Initiate vertical heatmap stack with dendrogram
zHTs <- module_dendrogram(modExpr.setSpecific.gtex, modTree = conTree.gtex, 
                          modColors = gtexMods, ht_height = unit(3, "cm")) 
maxZ <- 10 # Cap heatmap intensity at 10 to show modules below high quality
# Parse each network and compute heatmap and concatenate
for (tissue in names(gtex_zScoreMats.spectral)) { 
  # Color row annotation 
  numStats <- length(modStatSelect)
  groupAnno <- rowAnnotation(Region = rep(unique(droplevels(gtex_sample_summary[gtex_sample_summary$Region == tissue,]$Tissue_Subtype)), numStats),
                             col = list(Region = gtex_region_colors),
                             annotation_legend_param = list(Region = list(ncol = 1)), 
                             show_annotation_name = FALSE)
  # Create heatmap
  zHTs <- zHTs %v% module_quality_zScores_heatmap(gtex_zScoreMats.spectral[[tissue]], conTree.gtex, group = tissue, 
                                                  htColor = gtex_region_colors_short[tissue],  lgd_max = maxZ, 
                                                  rowAnno = groupAnno, column_select = modStatSelect, remove_row_names = TRUE)
}

# Save to file
#png(file=paste("/home/users/eggerj/tmp_transfer/modDensityZscores.spectral.gtex.png"), width = 6400, height = 4000, res=402)
draw(zHTs, heatmap_legend_side = "bottom", column_title = "Module Quality Scores")
#dev.off()


#######################################################################################
#
#    CLUSTER CONSENSUS MODULES USING CONSENSUS HIERARCHY
#
#######################################################################################

# Create dendrogram with module sizes object
mod_dend.gtex <- module_dendrogram(modExpr.setSpecific.gtex, modTree = conTree.gtex, 
                                   modColors = gtexMods, ht_height = unit(3, "cm"))

# Create second dendrogram with unique gene count per module
gene_count.gtex <- gene_counts_dendrogram(node_stats.gtexSpliceNet[[1]])

# Now create heatmap (cluster rows first using tree)
gtexHT <- module_heatmap_only(modExpr.setSpecific.gtex, cluster_rows = conTree.gtex)

# Save to file
png(file="/home/users/eggerj/tmp_transfer/ConsensusNet.gtex.png", width = 6800, height = 4800, res=420)
draw(mod_dend.gtex %v% gene_count.gtex %v% gtexHT, heatmap_legend_side = "bottom", column_title = 'Consensus Co-splicing Module Network')
dev.off()


#########################################################
#
#   PLOT GENE DISTRIBUTIONS ACROSS MODULE ASSIGNMENTS
#
#########################################################
p1 <- plot_gene_to_module_distributions(node_stats.gtexSpliceNet[[1]], countType = 'Genes')
p2 <- plot_gene_to_module_distributions(node_stats.gtexSpliceNet[[1]])
p3 <- plot_gene_to_module_distributions(node_stats.gtexSpliceNet[[1]], countType = 'Within')
ggarrange(p1, p2, p3, nrow = 1)

# Save to file
#ggsave(filename = "/home/users/eggerj/tmp_transfer/gene2ModuleCounts.gtex.png", width = 20, height = 6, dpi = 300)
#dev.off()



###########################################################################
#
#   LOOK AT PRESERVATION OF CONSENSUS MODULE NETWORKS ACROSS TISSUES
#
###########################################################################


# Start with subset of tissue types  --> heart tissues and two brain tissues
prezPlot <- plot_consensus_module_network_preservation(conMods.gtex[c("Hippocampus","Cortex","Atrial Appendage","Left Ventricle")], 
                                                       conTree.gtex, figSize = 3.85, txtSize = 11)
grid.arrange(prezPlot$full)

# Save to file
#png(file="/home/users/eggerj/tmp_transfer/prezSubset.gtex.png", width = 4600, height = 3800, res=410)
grid.arrange(prezPlot$full)
#dev.off()


# Now look at module preservation across all tissues using heatmap instead
prezPlot <- plot_consensus_module_network_preservation(conMods.gtex, conTree.gtex, figSize = 1.7, txtSize = 8.5)
#grid.arrange(prezPlot$full)  # <-- LARGE IMAGE

# Show only heatmap of densities
regions <- unique(names(conMods.gtex))

# Create colored heatmap annotations
col_anno <- HeatmapAnnotation(Tissue = gtex_sample_summary$Tissue_Type[match(regions, gtex_sample_summary$Region)],
                              Region = droplevels(gtex_sample_summary$Tissue_Subtype)[match(regions, gtex_sample_summary$Region)],
                              col = list(Tissue = gtex_tissue_colors, Region = gtex_region_colors), show_legend = TRUE, show_annotation_name = FALSE)
row_anno <- rowAnnotation(Tissue = gtex_sample_summary$Tissue_Type[match(regions, gtex_sample_summary$Region)],
                          Region = droplevels(gtex_sample_summary$Tissue_Subtype)[match(regions, gtex_sample_summary$Region)],
                          col = list(Tissue = gtex_tissue_colors, Region = gtex_region_colors), show_legend = FALSE, show_annotation_name = FALSE)
# Now create heatmap and clustering results
dPrezHT <- network_preservation_density_heatmap(prezPlot$Dvals, col_anno = col_anno, row_anno = row_anno)

# Save to file
#png(file="/home/users/eggerj/tmp_transfer/densityPrezHeatmap.gtex.png", width = 6800, height = 4800, res=420)
draw(dPrezHT, heatmap_legend_side = "bottom", column_title = 'Co-splicing Module Network Preservation Across Tissues')
#dev.off()



#######################################################################
#
#   PLOT MODULE SPLICING LEVELS ACROSS GTEX TISSUES AND TEST FOR 
#     DIFFERENTIAL SPLICING OF MODULES ACROSS TISSUE GROUPS
#
#######################################################################

# Sort samples by tissue group
df_gtex.sorted <- gtex_sample_summary[order(gtex_sample_summary$Tissue_Type),]
modExpr.gtex.sorted <- modExpr.gtex[rownames(df_gtex.sorted),]

# Create tissue group annotation colors for heatmap
anno_gtex <- rowAnnotation(Tissue = droplevels(df_gtex.sorted$Tissue_Type), 
                           col = list(Tissue = gtex_tissue_colors),
                           annotation_legend_param = list(Tissue = list(ncol = 1)), 
                           show_annotation_name = FALSE)


# Plot sample module expression as heatmap
exprHt <- module_expression_heatmap(modExpr.gtex.sorted, conTree.gtex, anno_gtex, row_clust = FALSE, 
                                    scale_mat = TRUE, mat_height = unit(12, "cm"))


# Create matching colored annotations for Wilcoxon rank sum test results
tissue_anno <- rowAnnotation(Tissue = sort(unique(df_gtex.sorted$Tissue_Type)),
                             col = list(Tissue = gtex_tissue_colors),
                             annotation_legend_param = list(Tissue = list(ncol = 1)), 
                             show_annotation_name = FALSE)

# Perform Wilcoxon rank sum test between tissue groups and show results with heatmap
tissueDS <- module_group_test_heatmap(modExpr.gtex.sorted, conTree.gtex, df_gtex.sorted, c(2), 
                                      annotation = tissue_anno, font_size = 10, mat_height = unit(5, "cm"))

# Save to file
#png(file=paste("/home/users/eggerj/tmp_transfer/diffSplicing.groupedTissues.gtex.png"), width = 6400, height = 4000, res=402)
draw(mod_dend.gtex %v% exprHt %v% tissueDS, heatmap_legend_side = "bottom", 
     column_title = "Differential Splicing of Network Modules Across Tissues")
#dev.off()







