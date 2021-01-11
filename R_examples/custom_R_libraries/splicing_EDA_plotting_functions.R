###################################################################
#
#    SAMPLE EDA PLOTTING FUNCTIONS
#
#     Group of functions for EDA plots of sample splicing data
#
###################################################################





# Function to plot PCs 1 through 4 and annotate given color code
plot_PCAs <- function(pc, summary.toPlot, batchPlot, pc_colors) {
  pcPlot <- ggarrange(
    autoplot(pc, x = 1, y = 4, data = summary.toPlot, colour = batchPlot) + 
      scale_color_manual(values = pc_colors),
    autoplot(pc, x = 2, y = 4, data = summary.toPlot, colour = batchPlot) + 
      scale_color_manual(values = pc_colors),
    autoplot(pc, x = 3, y = 4, data = summary.toPlot, colour = batchPlot) + 
      scale_color_manual(values = pc_colors),
    autoplot(pc, x = 1, y = 3, data = summary.toPlot, colour = batchPlot) + 
      scale_color_manual(values = pc_colors),
    autoplot(pc, x = 2, y = 3, data = summary.toPlot, colour = batchPlot) + 
      scale_color_manual(values = pc_colors),
    NA,
    autoplot(pc, x = 1, y = 2, data = summary.toPlot, colour = batchPlot) + 
      scale_color_manual(values = pc_colors),
    nrow = 3, ncol = 3,
    common.legend = TRUE, legend = "right")
  return(pcPlot)
}

#  Function to create alignment and LSV read coverage plot across sample set
#  - Requires a variety of coverage data in form of dataset sample summary (sample order must match)
#  - Drop coverage plot requires drop coverage data to be pre-computed given specific sample ordering
#  - If not using LSV data, parameters can be set to NULL
#
coverage_plot <- function(exprMat, summary_table, top_annot, lsv_read_counts, reads_axis, lsv_axis, bw = 0.6, 
                          drop_breaks = 2000, drop_height = unit(2.0, "cm"), quant_height = unit(1.0, "cm"), 
                          read_coverage_height = unit(3.5, "cm"), lsv_coverage_height = unit(2.5, "cm"), read_coverage_only = FALSE) {
  
  # Make into matrix
  lsv_read_mat <- as.matrix(lsv_read_counts[,rownames(summary_table)])
  
  # Calculate increase in LSVs in sample is removed (based on ordering)
  dropCov <- summary_table$drop_coverage*dim(exprMat)[2]
  dropCov95 <- (summary_table$drop_coverage_95*dim(exprMat)[2]) - (summary_table$drop_coverage*dim(exprMat)[2])
  allLSVs <- dim(exprMat)[2] - dropCov - dropCov95
  
  # Use if only including junctions reads prior to LSV quantification
  if (read_coverage_only == TRUE) {
    bottom_annot = HeatmapAnnotation(show_annotation_name = F,
                                     ReadCoverage = anno_barplot(cbind(summary_table$JunctionReads, 
                                                                       summary_table$AlignedReads - summary_table$JunctionReads),
                                                                 height = read_coverage_height, gp = gpar(fill = c(5,6), col = c(5,6)), 
                                                                 ylim = c(reads_axis[1], tail(reads_axis,1)),
                                                                 axis_param = list(side = "right",
                                                                                   at = reads_axis,
                                                                                   labels = reads_axis)),
                                     fill1 = anno_empty(border = FALSE, height = unit(0.2, "cm")),
                                     JunctionRatio = anno_barplot(cbind(summary_table$JunctionRatio, 1 - summary_table$JunctionRatio),
                                                                  height = read_coverage_height, gp = gpar(fill = c(5,6), col = c(5,6)), 
                                                                  axis_param = list(side = "right", 
                                                                                    at = c(0, 0.25, 0.5, 0.75, 1), 
                                                                                    labels = c(0, 0.25, 0.5, 0.75, 1.0))))
  } else {
    bottom_annot = HeatmapAnnotation(show_annotation_name = F,
                                     ReadCoverage = anno_barplot(cbind(summary_table$JunctionReads, 
                                                                       summary_table$AlignedReads - summary_table$JunctionReads),
                                                                 height = read_coverage_height, gp = gpar(fill = c(5,6), col = c(5,6)), 
                                                                 ylim = c(reads_axis[1], tail(reads_axis,1)),
                                                                 axis_param = list(side = "right",
                                                                                   at = reads_axis,
                                                                                   labels = reads_axis)),
                                     fill1 = anno_empty(border = FALSE, height = unit(0.2, "cm")),
                                     JunctionRatio = anno_barplot(cbind(summary_table$JunctionRatio, 1 - summary_table$JunctionRatio),
                                                                  height = read_coverage_height, gp = gpar(fill = c(5,6), col = c(5,6)), 
                                                                  axis_param = list(side = "right", 
                                                                                    at = c(0, 0.25, 0.5, 0.75, 1), 
                                                                                    labels = c(0, 0.25, 0.5, 0.75, 1.0))),
                                     fill2 = anno_empty(border = FALSE, height = unit(0.1, "cm")),
                                     LSVCoverage = anno_boxplot(lsv_read_mat, height = lsv_coverage_height, outline = FALSE, gp = gpar(fill = 'purple'), box_width = bw,
                                                                ylim = c(lsv_axis[1], tail(lsv_axis,1)),
                                                                axis_param = list(side = "right",
                                                                                  at = lsv_axis,
                                                                                  labels = lsv_axis)),
                                     fill3 = anno_empty(border = FALSE, height = unit(0.1, "cm")),
                                     QuantifiedLSVs = anno_barplot(summary_table$coverage_ratio,
                                                                   height = quant_height, gp = gpar(fill = 8, col = 1),
                                                                   axis_param = list(side = "right", 
                                                                                     at = c(0, 0.5, 1), 
                                                                                     labels = c(0, 0.5, 1))),
                                     fill4 = anno_empty(border = FALSE, height = unit(0.1, "cm")),
                                     DropSamples = anno_barplot(cbind(dropCov, dropCov95, allLSVs),
                                                                height = drop_height, gp = gpar(fill = c(2,3,4), col = c(2,3,4)), 
                                                                ylim = c(0, round_any(dim(exprMat)[2], 1000, f = ceiling)),
                                                                axis_param = list(side = "right",
                                                                                  at = c(seq(0,round_any(dim(exprMat)[2], 1000, f = ceiling),drop_breaks)), 
                                                                                  labels = c(seq(0,round_any(dim(exprMat)[2], 1000, f = ceiling),drop_breaks)))))
  }
  
  lgd1 = Legend(labels = c('Aligned Reads', 'Junction Spanning Reads'), title = "Sample Read Coverage", 
                legend_gp = gpar(fill = c(6,5)))
  lgd2 = Legend(labels = c('Non-junction Reads', 'Junction Spanning Reads'), title = "Junction Read Ratio",legend_gp = gpar(fill = c(6,5)))
  lgd3 = Legend(labels = c('Read Counts'), title = "Sample LSV Read Coverage", 
                legend_gp = gpar(fill = 'purple'))
  lgd4 = Legend(labels = c('Ratio'), title = "Ratio of LSVs Quantified", legend_gp = gpar(fill = 8))
  lgd5 = Legend(labels = c('100% of Samples','95% of Samples','Removed LSVs'), title = "LSV Analysis Set", 
                legend_gp = gpar(fill = c(2,3,4)), grid_width = unit(1, "cm"))
  
  if (read_coverage_only == TRUE) {
    legend_list <- list(lgd1, lgd2)
  } else {
    legend_list <- list(lgd1, lgd2, lgd3, lgd4, lgd5)
  }
  
  ht <- Heatmap(matrix(ncol = dim(summary_table)[1], nrow = 0), cluster_columns = FALSE, cluster_rows = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE,
                top_annotation = top_annot, bottom_annotation = bottom_annot,
                show_heatmap_legend = FALSE)
  
  return(list("ht" = ht, "legend_list" = legend_list))
}
