###################################################################################
#
#   Module Quality Evaluation Plotting Functions
#
#     Plotting functions for evaluating module formulation in 
#      co-splicing networks including bootstrapping and module
#      permutations
#
#####################################################################################

# Function for plotting heatmap of module preservation p-values from NetRep
module_quality_scores_heatmap <- function(zScoreMats, modTree, metrics = NULL, group = NULL, htColor = "purple") {
  
  # Transpose NetRep matrix and set module order to that of original module matrix
  zScoreMats <- t(zScoreMats[modTree$labels,])
  
  if (!is.null(metrics)) zScoreMats <- zScoreMats[metrics,]
  htMat <- -log10(zScoreMats)
  
  maxVal <- max(abs(min(htMat)), abs(max(htMat)))
  ht_colors = colorRamp2(c(min(htMat), 0, maxVal), c("white", "turquoise", "blue")) 
  
  ht_colors = colorRamp2(seq(quantile(htMat, 0.01),  quantile(htMat, 0.99), length = 3), 
                         c("white", "turquoise", "purple"))
  
  ht_colors = colorRamp2(c(min(htMat),max(htMat)), c("white", htColor))
  
  ht <- Heatmap(htMat, cluster_columns = modTree, cluster_rows = FALSE, 
                show_column_names = FALSE, show_row_names = TRUE, row_names_side = "left",
                name = "-log10 (Quality P-value)", col = ht_colors, row_title = group, row_title_side = "right",
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(8, "cm"),
                                            legend_height = unit(10, "cm"),  title_position = "topcenter"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.g", zScoreMats[i, j]), x, y, gp = gpar(fontsize = font_size))})
  return(ht)
}

# Function for plotting zScores of module preserveation 
module_quality_zScores_heatmap <- function(zScoreMats, modTree = NULL, group = NULL, htColor = "purple", font_size = 10, column_select = c(8:10),
                                           lgd_max = NULL, show_column_names = FALSE, column_title = NULL, rowAnno = NULL, remove_row_names = FALSE) {
  
  # Automate figure legends depending on select Zscore types
  if (length(column_select) == 1) {
    remove_row_names <- TRUE
    lgd_name <- colnames(zScoreMats)[column_select]
  } else {
    lgd_name <- 'Zscore'
  }
  
  # Transpose NetRep matrix and set module order to that of original module matrix
  if (!is.null(modTree)) {
    zScoreMats <- t(zScoreMats[modTree$labels,column_select, drop = remove_row_names])
  } else {
    zScoreMats <- t(zScoreMats[,column_select, drop = remove_row_names])
  }
  htMat <- zScoreMats
  
  # Set heatmap intensity   
  if (is.null(lgd_max)) {lgd_max <- round_any(max(htMat), 5, f = ceiling) }
  ht_colors = colorRamp2(c(0,lgd_max), c("white", htColor))
  
  if (is.null(modTree)) {modTree <- FALSE}
  ht <- Heatmap(htMat, cluster_columns = modTree, cluster_rows = FALSE, 
                show_column_names = show_column_names, show_row_names = TRUE, 
                row_names_side = "left", name = lgd_name, col = ht_colors, 
                right_annotation = rowAnno, column_title = column_title,
                heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(8, "cm"),
                                            legend_height = unit(10, "cm"),  title_position = "topcenter"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", zScoreMats[i, j]), x, y, gp = gpar(fontsize = font_size))})
  return(ht)
}

# Function for plotting boxplots of module quality zScores
module_quality_zScores_boxplots <- function(zScores, ylim = NULL, title = NULL, zSummary_only = TRUE) {
  
  if (zSummary_only == TRUE) {
    melt_scores <- melt(zScores[,10])
    melt_scores$zScore <- 'Zsummary'
    colnames(melt_scores) <- c('Score', 'zScore')
  } else{
    
    melt_scores <- melt(zScores[,c(8:10)])
    colnames(melt_scores) <- c('Module', 'zScore', 'Score')
  }
  
  if (is.null(ylim)) {ylim = round_any(max(melt_scores$Score), 5, ceiling)}
  
  p <- ggplot() + geom_boxplot(data = melt_scores, aes(x = zScore, y = Score, fill = zScore)) +
    #scale_x_discrete(limits = unique(melt_scores$Module)) +
    #scale_fill_manual(values = module_colors) +
    scale_y_continuous(limits = c(0,ylim), breaks = seq(0,ylim,5)) + 
    theme_bw() + theme(legend.position = "none") +
    ggtitle(title)
  return(p) 
  
}

# Function for plotting barplots of module quality zScores
module_quality_zScores_barplots <- function(zScores, ylim = NULL, title = NULL, zSummary_only = TRUE) {
  
  if (zSummary_only == TRUE) {
    melt_scores <- melt(zScores[,10])
    melt_scores$Module <- rownames(melt_scores)
    melt_scores$zScore <- 'Zsummary'
    colnames(melt_scores) <- c('Score', 'Module','zScore')
  } else{
    
    melt_scores <- melt(zScores[,c(8:10)])
    colnames(melt_scores) <- c('Module', 'zScore', 'Score')
  }
  
  if (is.null(ylim)) {ylim = round_any(max(melt_scores$Score), 5, ceiling)}
  
  p <- ggplot() + geom_bar(data = melt_scores, aes(x = reorder(Module,-Score), y = Score, fill = Module), stat = "identity") +
    #scale_x_discrete(limits = unique(melt_scores$Module)) +
    scale_fill_manual(values = module_colors) +
    scale_y_continuous(limits = c(0,ylim), breaks = seq(0,ylim,5)) + 
    theme_bw() + theme(legend.position = "none") +
    ggtitle(title)
  return(p) 
  
}

# Function to plot results from module stat permutations
module_permutation_stats_heatmap <- function(obs_stats, perm_stats, modTree, mStat, font_size = 10, annotation = NULL, mat_height = unit(4, "cm"), ht_name = NULL, maxVal = NULL) {
  
  nGroups <- dim(obs_stats[[mStat]])[1]
  nMods <- dim(obs_stats[[mStat]])[2]
  traitModPVals <- matrix(0, nGroups, nMods)
  traitModScores <- matrix(0, nGroups, nMods)
  for (i in 1:nGroups) {
    for (j in 1:nMods) {
      obs <- obs_stats[[mStat]][i,j]
      nullvals <- c()
      for (n in 1:length(perm_stats)) {
        nullvals <- c(nullvals, perm_stats[[n]][[mStat]][i,j])
      }  
      
      if (obs > mean(nullvals)) {
        pval <- (sum(nullvals > obs) + 1) / (length(nullvals) + 1)
        traitModScores[i, j] = abs(-log10(pval))
      } else {
        pval <- (sum(nullvals < obs) + 1) / (length(nullvals) + 1)
        traitModScores[i, j] = -(abs(-log10(pval)))
      }  
      traitModPVals[i,j] <- pval
      traitModScores[i, j] = -log10(pval)
    }
  }
  
  #traitModPVals[] <- p.adjust(traitModPVals, method = 'BH')
  
  rownames(traitModScores) <- rownames(obs_stats[[mStat]])
  htMat <- traitModScores
  
  if (is.null(maxVal)) {maxVal <- max(abs(min(htMat)), abs(max(htMat)))} 
  print(maxVal)
  ht_colors = colorRamp2(c(0, maxVal), c("white", "purple")) 
  
  if (is.null(annotation)) {name_rows = TRUE} else {name_rows = FALSE}
  
  ht <- Heatmap(htMat, cluster_columns = modTree, cluster_rows = FALSE, 
                show_column_names = FALSE, show_row_names = name_rows, row_names_side = "left", height = mat_height, row_title = ht_name,
                name = "P-Value (Group Median - Module Assignment Permutations)", col = ht_colors, left_annotation = annotation,
                heatmap_legend_param = list(legend_direction = "horizontal",
                                            legend_width = unit(8, "cm"),
                                            legend_height = unit(10, "cm"), 
                                            title_position = "topcenter"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.g", traitModPVals[i, j]), x, y, gp = gpar(fontsize = font_size))})
  return(ht)
}


# Create an example of how permutations are calculated
module_permutation_histogram <- function(mod_stats, mod_perms, group, module, mStat = 'Medians') {
  
  g <- which(rownames(mod_stats[[mStat]]) == group)
  m <- which(colnames(mod_stats[[mStat]]) == module)
  obs <- mod_stats[[mStat]][g,m]
  null_vals <- c()
  for (i in 1:length(mod_perms)) {
    null_vals <- c(null_vals, mod_perms[[i]][[mStat]][g,m])
  }
  
  h_title <- paste('Observed Median Module Expression vs. Randomly Permuted Module Assignments (P = ', as.character(length(mod_perms)), ');\nGroup = ', 
                   group, '; Module = ', module, sep = '')
  p <- ggplot() + 
    geom_histogram(aes(null_vals), bins = 500, color="black", fill="white") +
    geom_vline(aes(xintercept = obs), color = "red", size = .4) +
    scale_x_continuous(limits = c(min(c(null_vals, obs)),c(max(c(null_vals,obs))))) +
    theme(text = element_text(size=20)) +
    ggtitle(h_title)
  return(p)
}

# Plot results of module stability
plot_bootstrap_similarity_results <- function(boot_data) {
  
  melt_scores <- boot_data$bootstrap_sim_scores
  hIndex_df <- boot_data$module_Hindex
  
  p <- ggplot() + geom_boxplot(data = melt_scores, aes(x = Module, y = Score, fill = Module)) +
    scale_x_discrete(limits = unique(melt_scores$Module)) +
    scale_fill_manual(values = module_colors) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) + 
    theme_bw() + theme(legend.position = "none") +
    geom_point(data = hIndex_df, aes(x = Module, y = H_index), color = 'red') +
    ggtitle(boot_data$param_set)
  return(p) 
}

# Plot connectivity across bootstrap network clustering
plot_bootstrap_connectivity_results <- function(boot_data, netType) {
  
  net_colors <- c("unsigned" = "darkorange", "signed hybrid" = "darkgreen", "signed" = "blue")
  
  sf <- ggplot(boot_data$sft, aes(x = Power, y = SFT.R.sq)) + geom_boxplot(fill = net_colors[netType]) +
    ylab("Scale Free Fit") + scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05)) + 
    geom_hline(yintercept=0.85, color = "red", size = .3) +
    geom_hline(yintercept=0.9, color = "purple", size = .3) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Scale Free Topology Model Fit")
  
  meanK <- ggplot(boot_data$meanK, aes(x = Power, y = meanK)) + geom_boxplot(fill = net_colors[netType]) +
    ylab("Mean Connectivity") + #scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05)) + 
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Mean Connectivity")
  
  SFMC.plot <- ggarrange(sf, meanK, nrow = 1, ncol = 2)
  SFMC.plot <- annotate_figure(SFMC.plot, str_to_title(netType))
  
  return(SFMC.plot)
  
}

# Plot Zscore distributions from clustering of bootstrap networks
plot_bootstrap_zScore_results <- function(boot_data, ylimit = NULL, zScores = NULL, sortBy = NULL) {
  
  if (is.null(zScores)) {zScores = names(boot_data$zScores)}
  if (is.null(ylimit)) {ylimit = round_any(get_max_bootstrap_zScore(boot_data, zScores),10,f = ceiling)}
  else {ylimit = round_any(ylimit,10,f = ceiling)}
  
  if (is.null(sortBy)) {
    if (length(zScores) == 1) {
      mod_order <- names(sort(unlist(lapply(boot_data$zScores[[zScores]], median, na.rm = TRUE)), decreasing = TRUE))
    } else {
      mod_order <- names(sort(unlist(lapply(boot_data$zScores$Zsummary, median, na.rm = TRUE)), decreasing = TRUE))
    }
  } else {
    mod_order <- names(sort(unlist(lapply(boot_data$zScores[[sortBy]], median, na.rm = TRUE)), decreasing = TRUE))
  }  
  plt_list <- list()
  for (z in zScores) {
    df <- melt(boot_data$zScores[[z]])
    colnames(df) <- c('Module', 'Zscore')
    df <- transform(df, Module = factor(Module, levels = mod_order))
    zplot <- ggplot(df, aes(x = Module, y = Zscore, fill = Module)) + geom_boxplot() +
      scale_fill_manual(values = module_colors) + ggtitle(z) + ylab(z) +
      scale_y_continuous(limits = c(0,ylimit), breaks = seq(0,ylimit,10)) + 
      geom_hline(yintercept=10, color = "red", size = .3, linetype = "dashed") +
      geom_hline(yintercept=2, color = "purple", size = .3, linetype = "dashed") +
      theme_bw() + theme(legend.position = "none")
    plt_list[[z]] <- zplot 
  }
  zplots <- ggarrange(plotlist = plt_list, nrow = length(zScores))
  zplots <- annotate_figure(zplots, boot_data$param_set)
  return(zplots)
}


