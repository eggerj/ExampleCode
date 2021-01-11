#####################################################################
#
#   Network Inference Plotting Functions
#
#     Group of plotting functions for de novo network inference
#
######################################################################





##########################################################
#
#   Scale-free connectivity across beta value plots
#
##########################################################


# Plot scale-free fit and mean connectivity for multiple networks
multi_network_connectivity_plots <- function(rCut, powerDF, minDF, corType, netType, 
                                             group_name, legend, P) {
  
  # Network size
  netSize = paste('P = ', as.character(P), sep = '')
  
  SF <- ggplot(powerDF, aes(x=Power, y=SFT.R.sq)) + 
    geom_point(aes(col=Region)) + 
    scale_color_manual(values = gtex_region_colors_short) +
    scale_y_continuous(name = "Scale Free Fit", breaks = seq(0,1,.05), 
                       limits = c(0,1.0), minor_breaks = seq(0,1,.05)) +
    scale_x_continuous(name = "Soft Threshold (Power)", breaks = seq(1,max(powerDF$Power),1), 
                       limits = c(1,max(powerDF$Power)), minor_breaks = seq(1,max(powerDF$Power),1)) +
    geom_hline(yintercept=rCut, color = "red", size = .3)  + 
    geom_segment(aes(x = Power, xend = Power, y = SFT.R.sq, yend = -Inf, col=Region),
                 data = minDF, linetype="dotted", size=.4) + 
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Scale Free Topology Model Fit")
  
  
  MC <- ggplot(powerDF, aes(x=Power, y=mean.k.)) + 
    geom_point(aes(col=Region)) + 
    scale_color_manual(values = gtex_region_colors_short) +
    scale_y_continuous(name = "Mean Connectivity", 
                       breaks = seq(0,round_any(max(powerDF$mean.k.), 100, f = ceiling),100), 
                       limits = c(0,round_any(max(powerDF$mean.k.), 100, f = ceiling)), 
                       minor_breaks = seq(0,round_any(max(powerDF$mean.k.), 100, f = ceiling),100)) +
    scale_x_continuous(name = "Soft Threshold (Power)", breaks = seq(1,max(powerDF$Power),1), 
                       limits = c(1,max(powerDF$Power)), minor_breaks = seq(1,max(powerDF$Power),1)) +
    geom_segment(aes(x = Power, xend = Power, y = mean.k., yend = -Inf, col=Region),
                 data = minDF, linetype="dotted", size=.4) + 
    geom_segment(aes(x = Power, xend = -Inf, y = mean.k., yend = mean.k., col=Region),
                 data = minDF, linetype="dotted", size=.4) + 
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Mean Connectivity")
  
  SFMC.plot <- ggarrange(SF, MC, nrow = 1, ncol = 2, common.legend = TRUE, legend = legend)
  SFMC.plot <- annotate_figure(SFMC.plot, top = text_grob(paste(corType, netSize, netType, group_name, sep = ' - ')))
  
  return(SFMC.plot) 
  
}

# Plot scale-free fit and mean connectivty for single network (can use multiple correlation and network types)
single_network_connectivity_plots <- function(rCut, powerDF, minDF, corType, 
                                              netType, group_name, legend, P) {
  
  # Network size
  netSize = paste('P = ', as.character(P), sep = '')
  
  SF <- ggplot(powerDF, aes(x=Power, y=SFT.R.sq)) + 
    geom_point(aes(col=Network_Type)) + scale_color_manual(breaks = c("Unsigned", "Signed Hybrid", "Signed"),
                                                           values=c("darkorange", "darkgreen", "blue")) +
    scale_y_continuous(name = "Scale Free Fit", breaks = seq(0,1,.05), 
                       limits = c(0,1.0), minor_breaks = seq(0,1,.05)) +
    scale_x_continuous(name = "Soft Threshold (Power)", breaks = seq(1,max(powerDF$Power),1), 
                       limits = c(1,max(powerDF$Power)), minor_breaks = seq(1,max(powerDF$Power),1)) +
    geom_hline(yintercept=rCut, color = "red", size = .3)  + 
    geom_segment(aes(x = Power, xend = Power, y = SFT.R.sq, yend = -Inf, col=Network_Type),
                 data = minDF, linetype="dotted", size=.4) + 
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Scale Free Topology Model Fit")
  
  MC <- ggplot(powerDF, aes(x=Power, y=mean.k.)) + 
    geom_point(aes(col=Network_Type)) + scale_color_manual(breaks = c("Unsigned", "Signed Hybrid", "Signed"),
                                                           values=c("darkorange", "darkgreen", "blue")) +
    scale_y_continuous(name = "Mean Connectivity", 
                       breaks = seq(0,round_any(max(powerDF$mean.k.), 100, f = ceiling),100), 
                       limits = c(0,round_any(max(powerDF$mean.k.), 100, f = ceiling)), 
                       minor_breaks = seq(0,round_any(max(powerDF$mean.k.), 100, f = ceiling),100)) +
    scale_x_continuous(name = "Soft Threshold (Power)", breaks = seq(1,max(powerDF$Power),1), 
                       limits = c(1,max(powerDF$Power)), minor_breaks = seq(1,max(powerDF$Power),1)) +
    geom_segment(aes(x = Power, xend = Power, y = mean.k., yend = -Inf, col=Network_Type),
                 data = minDF, linetype="dotted", size=.4) + 
    geom_segment(aes(x = Power, xend = -Inf, y = mean.k., yend = mean.k., col=Network_Type),
                 data = minDF, linetype="dotted", size=.4) + 
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Mean Connectivity")
  
  SFMC.plot <- ggarrange(SF, MC, nrow =1, ncol = 2, common.legend = TRUE, legend = "bottom")
  SFMC.plot <- annotate_figure(SFMC.plot, top = text_grob(paste(corType, netSize, netType, group_name, sep = ' - ')))
  return(SFMC.plot)
}


# Plot topological overlap matrix
plot_TOM_matrix <- function(dissTOM, adj, spliceTree = NULL, dColors = NULL) {
  colnames(dissTOM) <- colnames(adj)
  rownames(dissTOM) <- colnames(adj)
  
  if (!is.null(dColors)) {
    nodes_to_colors <- data.frame(cbind(rownames(adj), dColors), row.names = rownames(adj))
    mod_annot <- HeatmapAnnotation(Modules = nodes_to_colors$dColors, 
                                   col = list(Modules = setNames(unique(dColors), unique(dColors))),
                                   show_legend = FALSE)
    mod_left <- rowAnnotation(Modules = nodes_to_colors$dColors, 
                              col = list(Modules = setNames(unique(dColors), unique(dColors))),
                              show_legend = FALSE, show_annotation_name = FALSE)
  } else {
    mod_annot <- NULL
    mod_left <- NULL
  }
  
  ht_colors = colorRamp2(c(0,.5,1), c("red", "yellow", "white"))
  
  htMat <- scale(dissTOM)
  col_fun = colorRamp2(seq(quantile(htMat, 0.01), quantile(htMat, 0.99), length = 3), 
                       c("red", "yellow", "white"))
  
  if (!is.null(spliceTree)) {
    plotTree <- spliceTree
    plotTree$height = (plotTree$height - min(plotTree$height))/(1.15 * (max(plotTree$height) - min(plotTree$height)))
  } else {
    plotTree <- FALSE
  }
  
  ht <- Heatmap(htMat, cluster_columns = plotTree, cluster_rows = plotTree, top_annotation = mod_annot, 
                show_row_names = FALSE, show_column_names = FALSE, col = col_fun, use_raster = FALSE,
                left_annotation = mod_left,
                column_dend_height = unit(3, "cm"), row_dend_width = unit(3, "cm"),
                #column_title = 'Consensus Topological Overlap Dissimilarity and Splicing Module Set',
                show_heatmap_legend = FALSE) 
  return (ht)
}



