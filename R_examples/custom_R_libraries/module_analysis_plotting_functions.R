###############################################################
#
#   Module Analysis Plotting Functions
#
#     This script contains all functions for characterizing
#      co-splicing modules in different datasets
#
################################################################


# This is edited source code from ComplexHeatamp in order to customize annotation functions
source('/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/shared_Rscripts/custom_anno_text_function.R')



# Function that computes correlations between traits of interest and module set
#  Can take in a dendrogram object so heatmap syncs with module-module correlation heatmap's dendrogram
module_trait_correlation_heatmap <- function(modExpr, modTree, traitsDF, trait_index = NULL, trait_names = NULL, 
                                             numericTraits = FALSE, font_size = 6, sortCors = FALSE, topN = NULL) {
  # Convert traits matrix to numeric (factor first)
  # First subset to specified traits
  if (!is.null(trait_index)) {
    traitsMat <- traitsDF[,trait_index]
  } else {
    traitsMat <- traitsDF
  }
  if (dim(as.matrix(traitsMat))[2] > 1) {
    if (numericTraits == FALSE) {
      traitsMat_factor <- apply(traitsMat, 2, function(x) as.numeric(as.factor(x)))
    } else {
      traitsMat_factor <- traitsMat
    }  
    trait_names <- colnames(traitsMat) 
  } else {
    if (numericTraits == FALSE) {
      traitsMat_factor <- as.numeric(as.factor(traitsMat))
    } else {
      traitsMat_factor <- traitsMat
    }  
    trait_names <- trait_names
  }
  
  moduleTraitCor = cor(modExpr, traitsMat_factor, use = "pairwise.complete.obs")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(modExpr))
  colnames(moduleTraitCor) <- trait_names
  moduleTraitCor <- abs(moduleTraitCor)

  if (sortCors == TRUE) {
    trait_order <- names(sort(apply(moduleTraitCor, 2, max), decreasing = TRUE))
    moduleTraitCor <- moduleTraitCor[,trait_order]
    moduleTraitPvalue <- moduleTraitPvalue[,trait_order]
    if (!(is.null(topN))) {
      moduleTraitCor <- moduleTraitCor[,1:topN]
      moduleTraitPvalue <- moduleTraitPvalue[,1:topN]
    }
  }
  
  ht_colors = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
  btht <- Heatmap(t(moduleTraitCor), cluster_columns = modTree, cluster_rows = FALSE, 
                  row_names_side = "left", show_column_names = FALSE,
                  name = "Correlation", col = ht_colors,
                  heatmap_legend_param = list(legend_direction = "horizontal",
                                              legend_width = unit(8, "cm"),
                                              legend_height = unit(10, "cm"), 
                                              title_position = "topcenter"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f\n(%.g)", t(moduleTraitCor)[i,j], t(moduleTraitPvalue)[i, j]), x, y, gp = gpar(fontsize = font_size))})
  
  return(btht)
}

# Function to plot module expression across sample set
module_expression_heatmap <- function(modExpr, modTree, annotation = NULL, row_clust = FALSE, scale_mat = TRUE, mat_height = NULL, 
                                      ht_lgd_name = "Splicing Module Zscore", adjust_ht = FALSE, htLims = NULL) {
  
  if (scale_mat == TRUE) {
    htMat <- scale(modExpr)
    if (is.null(htLims)) {htLims <- seq(quantile(htMat, 0.01), quantile(htMat, 0.99), length = 2)}
    #print(seq(quantile(htMat, 0.01), quantile(htMat, 0.99), length = 2))
    ht_colors = colorRamp2(htLims, c("midnightblue", "yellow"))
  } else {
    htMat <- modExpr
    ht_colors = colorRamp2(c(min(htMat), max(htMat)), c("midnightblue", "yellow"))
  }  
  
  if (adjust_ht == TRUE) {
    print(htMat[1:5,1:5])
    htMat <- htMat*2
    print(htMat[1:5,1:5])
    }
  
  lgd_params <- list(legend_direction = "horizontal",
                     legend_width = unit(8, "cm"),
                     legend_height = unit(10, "cm"), 
                     title_position = "topcenter")
  
  if (is.null(mat_height)) {
    ht <- Heatmap(htMat, cluster_columns = modTree, cluster_rows = row_clust, 
                  show_column_names = FALSE, show_row_names = FALSE,
                  left_annotation = annotation, 
                  name = ht_lgd_name, col = ht_colors,
                  heatmap_legend_param = lgd_params)
  } else {
    ht <- Heatmap(htMat, cluster_columns = modTree, cluster_rows = row_clust, 
                  show_column_names = FALSE, show_row_names = FALSE,
                  left_annotation = annotation, height = mat_height,
                  name = ht_lgd_name, col = ht_colors,
                  heatmap_legend_param = lgd_params) 
  }
  
  return(ht)
}


# Function to calculate 1 vs. all Wilcoxon rank sum tests based on module expression values
module_group_test_heatmap <- function(modExpr, modTree, traitsDF, trait_index, annotation = NULL, font_size = 6, mat_height = NULL) {
  
  # Break up trait into binary variables
  df <- data.frame(binarizeCategoricalVariable(traitsDF[,trait_index], 
                                               includePairwise = FALSE, includeLevelVsAll = TRUE, nameForAll = ""))
  rownames(df) <- rownames(traitsDF)
  
  # Create matrix that will store test results for each module and trait
  traitModScores = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2])
  traitModPVals = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2])
  
  for (t in 1:dim(df)[2]) {
    testSet <- rownames(df[df[,t] == 1,])
    restSet <- rownames(df[df[,t] == 0,])
    # Parse modules of current null expression matrix
    for (m in 1:dim(modExpr)[2]) {
      testVals <- modExpr[testSet,m]
      restVals <- modExpr[,m]
      testMean <- mean(testVals)
      restMean <- mean(restVals)
      # Calculate p-value
      p.val <- wilcox.test(testVals, restVals)$p.value
      traitModPVals[t,m] = p.val
      # Calculate -log10 of p-value and flip sign accordingly
      if (testMean < restMean) {
        traitModScores[t, m] = -(abs(-log10(p.val)))
      } else {
        traitModScores[t, m] = abs(-log10(p.val))
      }
    }
  }
  
  # Adjust for multiple testing
  traitModPVals[] <- p.adjust(traitModPVals, method = 'BH')
  for (t in 1:dim(traitModPVals)[1]) {
    for (m in 1:dim(traitModPVals)[2]) {
      if (traitModScores[t, m] < 0) {
        traitModScores[t, m] = -(abs(-log10(traitModPVals[t,m])))
      } else {
        traitModScores[t, m] = abs(-log10(traitModPVals[t,m]))
      }
    }
  }
  
  rownames(traitModScores) <- colnames(df)
  htMat <- traitModScores
  
  maxVal <- max(abs(min(htMat)), abs(max(htMat)))
  ht_colors = colorRamp2(c(-maxVal, 0, maxVal), c("blue", "white", "red")) 

  lgd_params <- list(legend_direction = "horizontal",
                     legend_width = unit(8, "cm"),
                     legend_height = unit(10, "cm"), 
                     title_position = "topcenter")
  
  if (is.null(annotation)) {name_rows = TRUE} else {name_rows = FALSE}
  
  if (is.null(mat_height)) {
  ht <- Heatmap(htMat, cluster_columns = modTree, cluster_rows = FALSE, 
                show_column_names = name_rows, show_row_names = name_rows, row_names_side = "left",
                name = "-log10(Adjusted P-value)", col = ht_colors, left_annotation = annotation,
                heatmap_legend_param = lgd_params,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.g", traitModPVals[i, j]), x, y, gp = gpar(fontsize = font_size))})
  } else {
    ht <- Heatmap(htMat, cluster_columns = modTree, cluster_rows = FALSE, 
                  show_column_names = name_rows, show_row_names = name_rows, row_names_side = "left",
                  name = "-log10(Adjusted P-value)", col = ht_colors, height = mat_height, left_annotation = annotation,
                  heatmap_legend_param = lgd_params,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.g", traitModPVals[i, j]), x, y, gp = gpar(fontsize = font_size))})
  }
  
  return(ht)
  
}



# Box plot version of module group Wilcoxon tests
module_group_test_boxplot <- function(modExpr, traitsDF, trait_index, bp_colors, numRows = 2) {

  # Break up trait into binary variables
  df <- data.frame(binarizeCategoricalVariable(traitsDF[,trait_index], 
                                               includePairwise = FALSE, includeLevelVsAll = TRUE, nameForAll = ""))
  rownames(df) <- rownames(traitsDF)
  
  # Make sure rows match
  modExpr <- modExpr[rownames(traitsDF),]
  
  # Make merged dataframe of modules
  mergedDF <- cbind(modExpr, traitsDF[,trait_index])
  colnames(mergedDF) <- c(colnames(modExpr), colnames(traitsDF)[trait_index])
  print(head(mergedDF))
  
  # Create matrix that will store test results for each module and trait
  traitModScores = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2])
  traitModPVals = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2])
  
  #compare_means(Black ~ Tissue_Type,  data = mergedDF, ref.group = ".all.", method = "wilcox.test")
  
  bp_list = list()
  for (m in colnames(modExpr)) {
    bp <- ggboxplot(mergedDF, x = colnames(traitsDF)[trait_index], y = m, color = colnames(traitsDF)[trait_index], 
                    add = "jitter", legend = "none", palette = bp_colors) +
      #rotate_x_text(angle = 45)+
      geom_hline(yintercept = mean(mergedDF[,m]), linetype = 2)+ # Add horizontal line at base mean
      #stat_compare_means(method = "kruskal", label.y = (max(mergedDF[,m]) + (max(mergedDF[,m])*0.1)))+        # Add global annova p-value
      stat_compare_means(method = "wilcox.test",
                         ref.group = ".all.", aes(label=..p.adj..)) 
    bp_list[[m]] <- bp
  }  
  return(ggarrange(plotlist = bp_list, common.legend = TRUE, align = "v"))
} 


# Function to create heatmap object of adjacency preservation between two networks
adjacency_preservation_heatmap <- function(cor1, cor2, square_unit = NULL, dSize = 8, anno_size = unit(0.5, "cm")) {
  
  # Calculate preservation between both networks of each module adjacency
  prez_mat <- 1 - abs(cor1 - cor2)
  
  modSize_annot <- HeatmapAnnotation(Module_Size = colnames(prez_mat), 
                                     col = list(Module_Size = setNames(colnames(prez_mat), colnames(prez_mat))),
                                     show_legend = FALSE, show_annotation_name = FALSE,  simple_anno_size = anno_size)
  modSize_left <- rowAnnotation(Module_Size = colnames(prez_mat), 
                                col = list(Module_Size = setNames(colnames(prez_mat), colnames(prez_mat))),
                                show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = anno_size)

  ht <- Heatmap(prez_mat, cluster_columns = FALSE, cluster_rows = FALSE, 
                bottom_annotation = modSize_annot, left_annotation = modSize_left,
                show_column_names = FALSE, show_row_names = FALSE,
                show_column_dend = FALSE, show_row_dend = FALSE,
                heatmap_width = square_unit, heatmap_height = square_unit,
                name = "Preservation", col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "blue")),
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(8, "cm"),
                                            legend_height = unit(10, "cm"),  title_position = "topcenter"))
  # BARPLOT
  # Barplot column means
  # prez_mat was "dp"
  sum_dp = mean(prez_mat[upper.tri(prez_mat)])  ### THIS IS THE D VALUE
  prez_mat2 <- prez_mat
  diag(prez_mat2) <- 0
  means = apply(prez_mat2, 2, sum)/(ncol(prez_mat2) - 1)
  
  dText <- paste('D', as.character(round(sum_dp, 3)), sep = ' = ') # For labeling D value above barplot
  
  # For barplot, must create empty heatmap to use (dummyht)
  bp <- HeatmapAnnotation(Preservation = anno_barplot(means, gp = gpar(fill = names(means)),
                                                      axis_param = list(at = c(0.0, 0.2, 0.4, 0.6, 0.8), 
                                                                        labels = c(0.0, 0.2, 0.4, 0.6, 0.8))),
                          height = square_unit, #width = square_unit,
                          show_annotation_name = FALSE, show_legend = FALSE)
  dummyht <- Heatmap(matrix(nrow = 0, ncol = length(means)),  cluster_rows = FALSE, cluster_columns = FALSE,
                     top_annotation = bp, column_title = dText, column_title_gp = gpar(fontsize = dSize),
                     #heatmap_width = square_unit) #, heatmap_height = square_unit)
                     width = square_unit)
  
  return_items = list()
  return_items[['pht']] <- ht 
  return_items[['bp']] <- dummyht
  return_items[['D']] <- sum_dp
  return(return_items)
}

# Function to plot just overall density preservation across networks in form of a heatmap 
network_preservation_density_heatmap <- function(dList, col_anno = NULL, row_anno = NULL) {
  
  dMat <- matrix(0, nrow = length(dList), ncol = length(dList), dimnames = list(names(dList),names(dList)))
  for (i in names(dList)) {
    for (j in names(dList)) {
      if (i == j) {dMat[i,j] <- 1.0} else {dMat[i,j] <- dList[[i]][[j]]}
    }
  }
  
  ht <- Heatmap(dMat, cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = FALSE, show_row_names = FALSE,
                top_annotation = col_anno, left_annotation = row_anno,
                column_dend_height = unit(3, "cm"), row_dend_width = unit(3, "cm"),
                name = "Preservation", col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "blue")),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.3f", dMat[i,j]), x, y, gp = gpar(fontsize = 11))},
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(8, "cm"),
                                            legend_height = unit(10, "cm"), title_position = "topcenter"))
  return(ht)

}

# Function to make comprehensive plot of consensus module networks across different systems
# NOTE: Make sure multiExpr module ordering is consistent with consensus tree (hclust object) ordering indices
#         --> Use consensusOrderMEs from WGCNA
plot_consensus_module_network_preservation <- function(multExpr, mTree, net_names = NULL, figSize = 5, txtSize = 12, 
                                                       bp_scale = 0.85, dend_scale = 0.8, anno_scale = 0.06) {
  
  # TEST PARAMETERS FOR DEBUGGING
  #multExpr <- list(multiConsMods[[1]], multiConsMods[[2]], multiConsMods[[3]])
  #names(multExpr) <- names(multiConsMods)[1:3]
  #net_names <- names(multExpr)
  #mTree <- dissTree.brain
  #figSize = 5.5
  
  # Size of each figure
  square_unit <- unit(figSize, "cm")
  bar_unit <- unit(figSize*bp_scale, "cm")
  dend_unit <- unit(figSize*dend_scale, "cm")
  anno_unit <- unit(figSize*anno_scale, "cm")
  
  if (is.null(net_names)) net_names <- names(multExpr)
  # Get module heatmaps and correlations for each network
  mod_heatmaps <- list()
  mod_cors <- list()
  mod_dends <- list()
  for (t in 1:length(net_names)) {
    nName <- net_names[t]
    nExpr <- multExpr[[t]]
    nExpr <- nExpr[,mTree$order]
    modht <- module_network_heatmap(nExpr, column_colors = "bottom",
                                    col_dend = FALSE, row_dend = FALSE, cluster_mods = FALSE,
                                    ht_width = square_unit, ht_height = square_unit, anno_size = anno_unit)
  mod_heatmaps[[nName]] <- modht$ht
  mod_cors[[nName]] <- modht$cor
  mod_dends[[nName]] <- module_dendrogram(modExpr = nExpr, modTree = NULL, show_sizes = FALSE, 
                                          ht_width = square_unit, ht_height = dend_unit)
  }
  
  # Create pair-wise sets for preservation plots
  prez_heatmaps = list()
  prez_barplots <- list()
  prez_densities <- list()
  for (t1 in 1:length(net_names)) {
    nName1 <- net_names[t1]
    prez_heatmaps[[nName1]] = list()
    prez_barplots[[nName1]] = list()
    prez_densities[[nName1]] = list()
    for (t2 in 1:length(net_names)) {
      nName2 <- net_names[t2]
      if (t1 != t2) {
        prez_heatmaps[[nName1]][[nName2]] <- adjacency_preservation_heatmap(mod_cors[[nName1]], mod_cors[[nName2]], 
                                                                            square_unit = square_unit, anno_size = anno_unit)$pht
        bpItems <- adjacency_preservation_heatmap(mod_cors[[nName1]], mod_cors[[nName2]], square_unit = bar_unit, dSize = txtSize)
        prez_barplots[[nName1]][[nName2]] <- bpItems$bp
        prez_densities[[nName1]][[nName2]] <- bpItems$D
      }
    }
  }
  
  # Store all plots in correct order
  plot_list <- list()
  
  # FIRST ROW WILL ALWAYS BE THE DENDROGRAMS IN ORDER
  for (n in net_names) {
    plot_list[[n]] <- grid.grabExpr(draw(mod_dends[[n]], show_heatmap_legend = FALSE, 
                                         column_title = n, column_title_gp = gpar(fontsize = txtSize)))
  }
  
  num_rows <- length(net_names)
  num_cols <- length(net_names)
  l <- length(plot_list)
  # SECOND ROW WILL BE FIRST NET CORRELATION PLOT, BARPLOT BETWEEN FIRST NET AND REMAINING NETS
  
  for (r in 1:num_rows) {
    row_name <- net_names[r]
    for(c in 1:num_cols) {
      col_name <- net_names[c]
      # If row == column: add correlation heatmap
      if (r == c) {
        n <- row_name
        plot_list[[l+1]] <- grid.grabExpr(draw(mod_heatmaps[[row_name]], show_heatmap_legend = FALSE, 
                                               column_title = col_name, column_title_gp = gpar(fontsize = txtSize)#,
                                               #row_title = row_name, row_title_gp = gpar(font_size = txtSize)
                                               ))
      } 
      # If row < column: add pairwise barplot
      else if (r < c) {
        plot_list[[l+1]] <- grid.grabExpr(draw(prez_barplots[[row_name]][[col_name]], show_heatmap_legend = FALSE))
      }
      # If row > column: add pairwise heatmap
      else {
        n <- paste(row_name, col_name, sep = ' / ')
        plot_list[[l+1]] <- grid.grabExpr(draw(prez_heatmaps[[row_name]][[col_name]], show_heatmap_legend = FALSE, 
                                               column_title = n, column_title_gp = gpar(fontsize = txtSize)#,
                                               #row_title = row_name, row_title_gp = gpar(font_size = txtSize)
                                               ))
      }
      l <- l + 1
    }
  }
  
  plt <- arrangeGrob(grobs = plot_list, ncols = length(net_names))
  
  # Create legends
  lgd_list <- list()
  col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
  lgd1 = Legend(col_fun = col_fun, title = "Edge Weight (Signed)", direction = "horizontal", title_position = "topcenter", 
                at = c(0,0.5,1), legend_width = square_unit)
  lgd_list[[1]] <- grid.grabExpr(draw(lgd1))
  
  col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "blue"))
  lgd2 = Legend(col_fun = col_fun, title = "Preservation", direction = "horizontal", title_position = "topcenter", 
                at = c(0,0.5,1), legend_width = square_unit)
  lgd_list[[2]] <- grid.grabExpr(draw(lgd2))
  
  lgd_plot <- arrangeGrob(grobs = lgd_list, nrows = 2)
  
  
  # Create arrangement matrix
  n <- length(multExpr)
  hlay <- c(rep(1,n), NA)
  for (i in 1:(n-1)) {
    hlay <- rbind(hlay, c(rep(1,n), NA))
  }
  hlay <- rbind(hlay, c(rep(1,n), 2))
  
  pltFinal <- arrangeGrob(plt, lgd_plot, layout_matrix=hlay)
  return_items <- list()
  return_items[['full']] <- pltFinal
  return_items[['Dvals']] <- prez_densities
  return(return_items)
}

# Function to plot module clustering dendrogram and 
# module color labels
# NOTE: For safety, module expression matrix must be included
#         - If tree (hclust object) is provided, tree order must match module expression matrix
#         - If tree not provided, module expression is used to create a tree specific to dataset  
module_dendrogram <- function(modExpr, modTree = NULL, modColors = NULL, show_sizes = TRUE, 
                              ht_width = NULL, ht_height = NULL) {

  # Make tree from expression matrix if no tree is provided
  if (is.null(modTree)) {
    # Compute module-module correlations
    moduleCor = cor(modExpr, modExpr, use = "p")
    modTree <- hclust(as.dist(1 - moduleCor), method = "average")
  }  
  
  if (show_sizes == TRUE) {
    # Calcuate module sizes for annotation (remove grey module first)
    modSizes <- data.frame(table(str_to_title(modColors[which(modColors != 'grey')])), row.names = 1)
    modSizes$Colors <- sapply(rownames(modSizes), FUN = isDark)
    modSizes <- modSizes[colnames(modExpr),]
    modSize_annot <- HeatmapAnnotation(ModuleSize = anno_text_custom(modSizes$Freq, rot = 0, location = 0.5, just = "center",
                                                               gp = gpar(fill = rownames(modSizes), col = modSizes$Colors,
                                                                         fontsize = 9),
                                                               height = unit(.75, "cm")),
                                       show_annotation_name = TRUE, annotation_name_side = "left")
  } else {
    modSize_annot <- HeatmapAnnotation(ModuleSize = colnames(modExpr), 
                                       col = list(ModuleSize = setNames(colnames(modExpr), colnames(modExpr))),
                                       show_legend = FALSE, show_annotation_name = FALSE)
  }
  
  if (is.null(ht_height)) { 
    d_height <- unit(10, "mm")
  } else {
    d_height <- ht_height
  }  
  return(Heatmap(matrix(nrow = 0, ncol = length(modTree$labels)), cluster_columns = modTree, cluster_rows = FALSE,
                 top_annotation = modSize_annot, show_column_names = FALSE, show_row_names = FALSE, 
                 column_dend_height = d_height, heatmap_width = ht_width))#, heatmap_height = ht_height))
}


# Function to create dendrogram but with unique genes in each module for annotation
#     Requires a node_stats dataframe
gene_counts_dendrogram <- function(node_stats, modTree = FALSE, ht_height = NULL, ht_width = NULL) {
  
  # Calcuate module sizes for annotation (remove grey module first)
  modSizes <- node_stats %>%  group_by(Module) %>% 
    summarise(Freq = n_distinct(GeneID)) %>% data.frame() 
  modSizes$Colors <- sapply(modSizes$Module, FUN = isDark)
  modSizes <- modSizes[order(as.character(modSizes$Module)),]
  modSize_annot <- HeatmapAnnotation(GeneCount = anno_text_custom(modSizes$Freq, rot = 0, location = 0.5, just = "center",
                                                                   gp = gpar(fill = modSizes$Module, col = modSizes$Colors,
                                                                             fontsize = 9),
                                                                   height = unit(.75, "cm")),
                                     show_annotation_name = TRUE, annotation_name_side = "left")
  
  if (is.null(ht_height)) { 
    d_height <- unit(10, "mm")
  } else {
    d_height <- ht_height
  }  
  return(Heatmap(matrix(nrow = 0, ncol = nrow(modSizes)), cluster_columns = modTree, cluster_rows = FALSE,
                 top_annotation = modSize_annot, show_column_names = FALSE, show_row_names = FALSE, 
                 column_dend_height = d_height, heatmap_width = ht_width))
  
}


# Function to reate heatmap using correlations between module set
module_network_heatmap <- function(modExpr, modTree = NULL, modColors = NULL, show_sizes = FALSE, column_colors = "top", 
                                   ht_width = NULL, ht_height = NULL, col_dend = TRUE, row_dend = TRUE, cluster_mods = TRUE,
                                   anno_size = unit(0.5, "cm"), d_height= unit(1, "cm"), d_width = unit(10, "mm")) {
  # Compute module-module correlations
  moduleCor = cor(modExpr, modExpr, use = "p")
  # If no tree provided, then make one using correlation distance
  # This would actually be better suited for the dendrogram only function
  if (is.null(modTree)) {
    if (cluster_mods == TRUE) {
      modTree <- hclust(as.dist(1 - moduleCor), method = "average")
    } else { 
      modTree = FALSE
    }
  }  
  # Make correlations a signed network
  moduleCor <- (moduleCor + 1) / 2
  
  
  if (show_sizes == TRUE) {
    # Calcuate module sizes for annotation (remove grey module first)
    modSizes <- data.frame(table(str_to_title(modColors[which(modColors != 'grey')])), row.names = 1)
    modSizes$Colors <- sapply(rownames(modSizes), FUN = isDark)
    modSizes <- modSizes[colnames(modExpr),]
    modSize_annot <- HeatmapAnnotation(ModuleSize = anno_text_custom(modSizes$Freq, rot = 0, location = 0.5, just = "center",
                                                               gp = gpar(fill = rownames(modSizes), col = modSizes$Colors,
                                                                         fontsize = 9),
                                                               height = unit(.75, "cm")),
                                       show_annotation_name = TRUE, annotation_name_side = "left")
  } else {
    modSize_annot <- HeatmapAnnotation(ModuleSize = colnames(modExpr), 
                                       col = list(ModuleSize = setNames(colnames(modExpr), colnames(modExpr))),
                                       show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = anno_size)
  }
  modSize_left <- rowAnnotation(ModuleSize = colnames(modExpr), 
                                col = list(ModuleSize = setNames(colnames(modExpr), colnames(modExpr))),
                                show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = anno_size)
  
  ht_colors = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
  if (column_colors == "top") {
    ht <- Heatmap(t(moduleCor), cluster_columns = modTree, cluster_rows = modTree, 
                  top_annotation = modSize_annot,
                  left_annotation = modSize_left,
                  column_dend_height = d_height, row_dend_width = d_width,
                  show_column_names = FALSE, show_row_names = FALSE,
                  show_column_dend = col_dend, show_row_dend = row_dend,
                  heatmap_width = ht_width, heatmap_height = ht_height,
                  name = "Edge Weight - Signed", col = ht_colors,
                  heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(8, "cm"),
                                              legend_height = unit(10, "cm"),  title_position = "topcenter"))
  } else {
    if (column_colors == FALSE) modSize_annot = NULL
    ht <- Heatmap(t(moduleCor), cluster_columns = modTree, cluster_rows = modTree, 
                  bottom_annotation = modSize_annot,
                  left_annotation = modSize_left,
                  column_dend_height = d_height, row_dend_width = d_width,
                  show_column_names = FALSE, show_row_names = FALSE,
                  show_column_dend = col_dend, show_row_dend = row_dend,
                  heatmap_width = ht_width, heatmap_height = ht_height,
                  name = "Edge Weight - Signed", col = ht_colors,
                  heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(8, "cm"),
                                              legend_height = unit(10, "cm"),  title_position = "topcenter"))
  }
  return_items = list()
  return_items[['ht']] <- ht
  return_items[['cor']] <- t(moduleCor)
  return(return_items)
}


# New alternative function that no longer automates module clustering so precomputed dendrogram
#  can be used with multiple heatmap objects. TODO: Merge with above
module_heatmap_only <- function(modExpr, cluster_columns = FALSE, cluster_rows = FALSE, modColors = NULL, annotate_left_colors = TRUE, 
                                         ht_width = NULL, ht_height = NULL, col_dend = FALSE, row_dend = TRUE,
                                         anno_size = unit(0.5, "cm"), d_height= unit(1, "cm"), d_width = unit(10, "mm")) {
  
  # Compute module-module correlations
  moduleCor = cor(modExpr, modExpr, use = "p")
  # Make correlations a signed network
  moduleCor <- (moduleCor + 1) / 2
  
  
  if (annotate_left_colors == TRUE) {
    modSize_left <- rowAnnotation(ModuleSize = colnames(modExpr), 
                                  col = list(ModuleSize = setNames(colnames(modExpr), colnames(modExpr))),
                                  show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = anno_size)
  }
  
  ht <- Heatmap(t(moduleCor), cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
                left_annotation = modSize_left,
                column_dend_height = d_height, row_dend_width = d_width,
                show_column_names = FALSE, show_row_names = FALSE,
                show_column_dend = col_dend, show_row_dend = row_dend,
                heatmap_width = ht_width, heatmap_height = ht_height,
                name = "Edge Weight - Signed", col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
                heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(8, "cm"),
                                            legend_height = unit(10, "cm"),  title_position = "topcenter"))
  
  return(ht)
}


# Function to produce heatmap showing correlations between kME and kWithin at once and produces a heatmap
plot_group_module_node_stat_correlations_heatmap <- function(node_stats, x, y, modules = NULL, grouping = 'Module', 
                                                             row_anno = row_anno, 
                                                             title = 'Correlation of Node Intra-modular Connectivty & Node Module Membership',
                                                             includeGrey = FALSE) {
  
  # TODO: Make absolute value, labels, and axis limits parameters
  print('Taking absolute value of x stat')
  
  if (is.null(modules)) {modules <- sort(unique(node_stats[[1]][,grouping]))}
  if (includeGrey == FALSE) {modules <- modules[modules != 'Grey']}
  
  groups <- names(node_stats)
  corMat <- matrix(0, nrow = length(groups), ncol = length(modules), dimnames = list(groups, modules))
  pMat <- matrix(0, nrow = length(groups), ncol = length(modules), dimnames = list(groups, modules))
  
  for (g in groups) {
    group_stats <- node_stats[[g]]  
    for (mod in modules) {
      mod_node_stats <- group_stats[group_stats[,grouping] == mod,]
      mod_node_stats[,x] <- abs(mod_node_stats[,x])
      modCor = cor(mod_node_stats[,x], mod_node_stats[,y], use = "pairwise.complete.obs")
      modpVal = corPvalueStudent(modCor, nrow(group_stats))
      corMat[g,mod] <- modCor
      pMat[g,mod] <- modpVal
    }
  }
  
  ht <- Heatmap(corMat, cluster_columns = FALSE, cluster_rows = FALSE, #column_dend_height = unit(3,'cm'), row_dend_width = unit(3,'cm'), 
                column_title = title, name = 'Correlation', 
                left_annotation = row_anno, show_row_names = FALSE,
                show_column_names = FALSE, #top_annotation = col_anno, 
                col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")), 
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.3f\n(%.g)", corMat[i,j], pMat[i, j]), x, y, gp = gpar(fontsize = 9))},
                heatmap_legend_param=list(legend_direction="horizontal",legend_width=unit(8,"cm"),
                                          title_position="topcenter"))
  return(ht)
  
}

# Function to plot distibution for number of modules in which a gene is represented in splicing network or
#  distribution for counts of each gene per module
#  --> countType can be 'Across' for modules per gene or 'Within' for gene counts within modules
plot_gene_to_module_distributions <- function(node_stats, countType = 'Across') {
  
  totalMods <- as.character(length(unique(node_stats$Module)))
  totalSVRs <- as.character(length(unique(node_stats$SVR)))
  totalGenes <- as.character(length(unique(node_stats$GeneID)))
  
  if (countType == 'Across') {
    module_counts <- node_stats %>% group_by(GeneID) %>%
      summarise(count = n_distinct(Module)) %>% data.frame()
    xName <- 'Modules Per Gene'
    plotName <- 'B) Distribution of Module Assignments Per Gene'
    
  } else if (countType == 'Genes') {
    
    module_counts <- node_stats %>% 
      count(GeneID, name = 'count') %>% data.frame()
    xName <- 'SVRs Per Gene'
    plotName <- 'A) Distribution of SVRs Per Gene'
  }
  
  else {
    module_counts <- node_stats %>% 
      count(Module, GeneID, name = 'count') %>% data.frame()
    xName <- 'Gene Occurrences Per Module'
    plotName <- 'C) Distribution of Gene Occurrences Within Modules'
  }

  p <- ggplot(module_counts, aes(x = count)) + geom_bar() +
    geom_text(stat='count', aes(label=..count..), vjust=-1) +
    scale_x_continuous(breaks = c(1:max(module_counts$count))) +
    scale_y_continuous(breaks = seq(0, round_any(max(table(module_counts$count)), 500, f = ceiling), 500),
                       limits = c(0, round_any(max(table(module_counts$count)), 500, f = ceiling))) +
    xlab(xName) + ylab('Count') + ggtitle(plotName) +
    theme_bw()
  return(p)
}

# Function for clustering multiple networks based on node values (e.g. connectivity)
#   Takes a dataframe with each row being a node and each column being a group
#     Default is that last column has module assignments and is removed before clustering
group_clustering_by_network_nodes_heatmap <- function(node_matrix, distance_method = "kendall", removeLastColumn = TRUE, column = NULL,
                                                      col_anno = NULL, row_anno = NULL, 
                                                      title = 'Concordance of Intra-modular Connectivity Across Tissue Networks',
                                                      ht_name = "Kendall's Rank Correlation") {
  
  if (removeLastColumn == FALSE) {
    if (!is.null(column)) {
      node_matrix <- node_matrix[,-column]
    } 
  } else {
    node_matrix <- node_matrix[,-dim(node_matrix)[2]]
  }
  
  connCor = cor(node_matrix, method = distance_method)
  connCor.pvalue = corPvalueStudent(connCor, nrow(node_matrix))
  connCorTree <- hclust(as.dist(1 - connCor), method = "average")
  
  ht <- Heatmap(connCor, cluster_columns = connCorTree, cluster_rows = connCorTree, 
                column_dend_height = unit(3,'cm'), row_dend_width = unit(3,'cm'), 
                column_title = title, name = ht_name, 
                top_annotation = col_anno, left_annotation = row_anno, show_column_names = FALSE, show_row_names = FALSE,
                col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")), 
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.3f\n(%.g)", connCor[i,j], connCor.pvalue[i, j]), x, y, gp = gpar(fontsize = 9))},
                heatmap_legend_param=list(legend_direction="horizontal",legend_width=unit(8,"cm"),
                                          title_position="topcenter"))
  draw(ht, heatmap_legend_side = "bottom")
  return(ht)
  
}


# Function to create correlation plots between node characteristics
#   Ex) module-membership values vs. intramodular connectivity
plot_module_node_stat_correlations <- function(node_stats, x, y, modules = NULL, grouping = 'Module', includeGrey = FALSE) {
  
  print('Taking absolute value of x stat')
  
  if (is.null(modules)) {modules <- sort(unique(node_stats[,grouping]))}
  if (includeGrey == FALSE) {modules <- modules[modules != 'Grey']}
  
  plot_list <- list()
  for (mod in modules) {
    mod_node_stats <- node_stats[node_stats[,grouping] == mod,]
    mod_node_stats[,x] <- abs(mod_node_stats[,x])
    corplot <- ggscatter(mod_node_stats, x = x, y = y, 
                         color = module_colors[mod],
                         add = "reg.line", conf.int = FALSE, fullrange = FALSE, ylim = c(0.0,1),
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coeff.args = list(method = "pearson", label.x.npc = "center", label.y.npc = "top", 
                                               label.sep = "\n", size = 6),
                         xlab = x, ylab = y)
    plot_list[[mod]] <- corplot
  }
  
  x <- length(plot_list)
  cols = round(sqrt(x),0)
  rows = ceiling(x/cols)
  
  modplot <- ggarrange(plotlist = plot_list, ncol = cols, nrow = rows)
  
  return(modplot)  
  
}

# Function for making box plots of gene length by modules
plot_module_gene_size <- function(module_geneLengths, title = "Distribution of Node Gene Length By Module") {
  
  xlabs <- table(module_geneLengths$Module)
  
  module_gene_length_boxplot <- ggplot(module_geneLengths, aes(x=Module,y=GeneLength, fill=Module)) +
    geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) + scale_fill_manual(values = module_colors) + 
    scale_y_continuous(breaks=seq(0,45000,5000), limits = c(0,45000)) + 
    scale_x_discrete(labels = xlabs) +
    xlab('Module Gene Count') + 
    ggtitle(title) +
    theme(text = element_text(size=14))
  return(module_gene_length_boxplot)
  
}


# Function to create a heatmap indicating ratio of shared hubs
group_specific_node_ratio_heatmap <- function(subset_nodes, background_nodes, module_node_assignments, group_name = c("Tissue"), net_groups = NULL, group_by_column = c("Module"), modules = NULL,
                                              title = 'Tissue-specific Module Hubs', ht_title = 'Ratio of Tissue-specific Hubs in Module', 
                                              modTree = FALSE, group_annotation = NULL, includeGrey = FALSE) {
  
  
  if (is.null(net_groups)) { net_groups <- unique(subset_nodes[,group_name])}
  
  if (is.null(modules)) {modules <- sort(unique(subset_nodes[,group_by_column]))}
  
  if (includeGrey == FALSE) {modules <- modules[modules != 'Grey'] }
  
  # Findinding Group-specific hubs
  groupSpecificHubMat <- matrix(0, nrow = length(net_groups), ncol <- length(modules),
                                dimnames = list(net_groups, str_to_title(modules)))
  groupSpecificHubRatio <- matrix(0, nrow = length(net_groups), ncol <- length(modules),
                                  dimnames = list(net_groups, str_to_title(modules)))
  
  groupSpecificNodeList <- list()
  for (t in 1:length(net_groups)) {
    nGroup <- net_groups[t]
    groupSpecificNodeList[[nGroup]] <- list() 
    for (m in 1:length(modules)) {
      mod <- str_to_title(modules)[m]
      test_hubs <- rownames(subset_nodes[subset_nodes[,group_by_column] == mod & subset_nodes[,group_name] == nGroup,])
      mod_background_hubs <- rownames(background_hubs[background_hubs[,group_by_column] == mod & background_hubs[,group_name] != nGroup,])
      uniqueHubs <- test_hubs[!(test_hubs %in% mod_background_hubs)]
      groupSpecificNodeList[[nGroup]][[mod]] <- uniqueHubs
      numUnique <- length(uniqueHubs)
      groupSpecificHubMat[t,m]  <- numUnique
      groupSpecificHubRatio[t,m] <- numUnique / table(module_node_assignments)[str_to_lower(mod)]
      
    }
  }
  
  ht <- Heatmap(groupSpecificHubRatio, cluster_columns = modTree, cluster_rows = FALSE, 
                column_dend_height = unit(2,'cm'), column_title = title,
                left_annotation = group_annotation, show_column_names = FALSE, show_row_names = FALSE,
                name = ht_title, col = colorRamp2(c(0, 0.06), c("white", "red")),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.1d\n(%.3f)", groupSpecificHubMat[i,j], groupSpecificHubRatio[i, j]), x, y, gp = gpar(fontsize = 9))},
                heatmap_legend_param=list(legend_direction="horizontal",legend_width=unit(8,"cm"),
                                          title_position="topcenter"))
  return_list <- list()
  return_list[['ht']] <- ht
  return_list[['nodeList']] <- groupSpecificNodeList
  return(return_list)
  
}

# Function to plot heatmap showing overlap of nodes between groups
node_overlap_heatmap <- function(node_groups, grouping = c('Tissue'), node_column = c('SVR'), row_anno = NULL, col_anno = NULL,
                                 cluster_groups = TRUE, title = 'Intra-modular Hub Node Similarity Across Tissues') {
  
  nGroups <- unique(node_groups[,grouping])
  nodeOverlaps <- matrix(0, nrow = length(nGroups), ncol = length(nGroups),dimnames = list(nGroups, nGroups))
  hubSim <- matrix(0, nrow = length(nGroups), ncol = length(nGroups),dimnames = list(nGroups, nGroups))
  
  for (i in nGroups) {
    for (j in nGroups) {
      iNodes <- node_groups[node_groups[,grouping] == i,]$SVR
      jNodes <- node_groups[node_groups[,grouping] == j,]$SVR
      sharedNodes <- intersect(iNodes, jNodes)
      nodeOverlaps[i,j] <- length(sharedNodes)
      hubSim[i,j] <- length(sharedNodes) / min(length(iNodes), length(jNodes))
    }
  }
  
  if (cluster_groups == TRUE) {
    cluster_columns = TRUE
    cluster_rows = TRUE
  } else {
    cluster_columns = FALSE
    cluster_rows = FALSE
  }
  
  # TODO: Make sure tissues match before adding annotations
  ht <- Heatmap(hubSim, cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
                column_dend_height = unit(2,'cm'), column_title = title,
                left_annotation = row_anno, top_annotation = col_anno, show_column_names = FALSE, show_row_names = FALSE,
                name = 'Similarity', col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.1d\n(%.3f)", nodeOverlaps[i,j], hubSim[i, j]), x, y, gp = gpar(fontsize = 9))},
                heatmap_legend_param=list(legend_direction="horizontal",legend_width=unit(8,"cm"),
                                          title_position="topcenter"))
  return(ht)
}


# Function create venn diagram showing gene list overlap of modules
module_gene_overlap <- function(node_stats, modules = NULL, grouping = c("Module"), geneID = c("GeneID"), includeGrey = FALSE,
                                alpha_val = 0.8, font_size = 2, plt_labels = NULL, title = "Module Gene Sets", color_assign = NULL) {
  if (is.null(modules)) {
    modules = str_to_title(sort(node_stats[,grouping]))
  } else {
    modules = str_to_title(modules)
  }
  if (includeGrey == FALSE) {modules = modules[modules != 'Grey']}
  
  module_gene_list <- list()
  mod_labels <- c()
  for (mod in modules) {
    module_nodes <- node_stats[node_stats[,grouping] == mod,]
    mod_size <- length(rownames(module_nodes))
    module_gene_list[[mod]] <- unique(module_nodes[,geneID])
    num_genes <- length(module_gene_list[[mod]])
    mod_labels <- c(mod_labels, paste(mod, '\n', as.character(num_genes), ' (', as.character(mod_size), ')', sep = ''))
  }
  # Venn diagram colors
  if (is.null(color_assign)) {
    venn_colors <- module_colors[modules]
  } else {
    venn_colors <- color_assign[modules]
  }
  if ('Black' %in% names(venn_colors)) {venn_colors['Black'] = 'white'}
  venn_colors <- unname(venn_colors)
  a <- rep(alpha_val, length(modules))
  
  if (!(is.null(plt_labels))) {mod_labels = plt_labels}
  
  venn.plot <- venn.diagram(module_gene_list, NULL, fill=venn_colors, alpha=a, cex = font_size, category.names=mod_labels, 
                            cat.cex = font_size, main=title, main.cex = font_size)
  
  return(venn.plot)
}