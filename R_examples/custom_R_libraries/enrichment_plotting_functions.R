####################################################################################
#
#   Functional Enrichment Plotting Functions
#
#     Functions for plotting results of gene ontology and pathway enrichment
#
###################################################################################


# Function to create dotplot of module enrichment terms (similar to dotplot of clusterProfiler)
enrichment_dotplot <- function(rez, module, topN = NULL, plt_title = module, p_thres = 0.05) {
  
  gene_hits <- as.numeric(str_split_fixed(rez$GeneRatio, '/', 2)[1,2])
  x <- data.frame(Description = rez$Description, p.adjust = rez$p.adjust, Count = rez$Count, 
                  GeneRatio = rez$Count/gene_hits)
  x <- x[order(x$p.adjust),]
  x <- x[x$p.adjust < p_thres,]
  if (!(is.null(topN))) {
    n <- min(dim(x)[1], topN)
  } else {
    n <- dim(x)[1]
  }  
  x <- x[1:n,]
  
  p <- ggplot(x, aes(x= reorder(Description, Count), y=Count)) +   
    geom_point(aes(colour=p.adjust, size=Count)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(y = paste("Associated Genes (", as.character(gene_hits), ")",sep = ""), 
         fill = "p.adjust", x = "") +
    theme_light() + 
    ggtitle(plt_title) +
    theme(plot.title = element_text(size = 16),
          legend.key.size = unit(2, "line"),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    scale_y_continuous(breaks=seq(min(x$Count), max(x$Count), 1), limits=c(min(x$Count), max(x$Count))) + 
    coord_fixed(ratio = 0.5)+
    coord_flip()    
  
  return(p)
  
} 

# Function to plot number of modules each enrichment term belongs to
plot_enrichment_to_module_counts <- function(rez_list, modules = NULL, p_thres = 0.05, startN = 1, stopN = NULL, y_lab = 'Enriched Terms', minGenes = 1,
                                             title = 'Enriched Terms By Module Count', grep_term = NULL, color_vector = NULL, xLab = c('Module')) {
  if (is.null(modules)) {
    modules = str_to_title(sort(names(rez_list)))
  } else {
    modules = str_to_title(modules)
  }
  
  # Filter terms by gene count
  rez_list <- filter_enriced_terms_by_gene_count(rez_list[modules], minGenes = minGenes)
  
  all_mod_rez <- data.frame(rez_list[[modules[1]]])
  all_mod_rez <- all_mod_rez[FALSE,]
  # Fill dataframe using (up to) top N enriched terms for each module
  for (mod in modules) {
    rez <- data.frame(rez_list[[mod]])
    rez <- rez[rez$p.adjust < p_thres,]
    n <- dim(rez)[1]
    if (n >= 1) {
      gene_hits <- as.numeric(str_split_fixed(rez$GeneRatio, '/', 2)[1,2])
      x <- data.frame(Description = rez$Description, p.adjust = rez$p.adjust, Count = rez$Count, 
                      GeneRatio = rez$Count/gene_hits, ID = rez$ID)
      x[,xLab] <- rep(mod, dim(x)[1])
      x <- x[order(x$p.adjust),]
      x <- x[1:n,]
      all_mod_rez <- rbind(all_mod_rez, x)
    }
  }
  
  if (!is.null(grep_term)) {
    select_terms <- grep(grep_term,all_mod_rez$Description,ignore.case=TRUE,value=FALSE)
    all_mod_rez <- all_mod_rez[select_terms,]
  }
  
  print(paste(as.character(length(unique(all_mod_rez$Description))), ' Enriched Terms', sep = ''))
  
  # Get module counts for each enriched term for subsetting
  module_counts <- sort(table(all_mod_rez[,c('Description')]),decreasing = TRUE)
  if (is.null(stopN)) {
    n <- dim(all_mod_rez)[1]
  } else {
    n <- min(stopN, dim(all_mod_rez)[1])
  }
  startN <- min(n, startN)
  top_terms <- names(module_counts)[startN:n]
  all_mod_rez <- all_mod_rez[all_mod_rez$Description %in% top_terms,]
  
  print(paste('Showing ', as.character(length(unique(all_mod_rez$Description))), ' Enriched Terms', sep = ''))
  
  if (is.null(color_vector)) {color_vector = module_colors}
  
  p <- ggplot(all_mod_rez) + geom_bar(aes(y = fct_rev(fct_infreq(Description)), 
                                          fill = all_mod_rez[,xLab])) +
    labs(fill=xLab[1])  +
    scale_fill_manual(values = color_vector) + theme_bw() +
    scale_x_continuous(breaks=seq(0,length(modules), 1), limits=c(0, length(modules))) +
    labs(y = y_lab, x = paste(xLab, 'Count', sep = ' ')) +
    ggtitle(title)
  return(p)
  
}

# Function for plotting list of genes associated with a given term
#   Term is based on a grep string from enrichment descriptions
#   Can be a list of terms
plot_term_specific_gene_counts <- function(modDFs, grep_term, p_thres = 0.05, split_column = 1, ylab = "Associated Genes", title = NULL, 
                                           color_vector = NULL, xLab = c('Module')) {
  
  if (length(grep_term) > 1) {grep_term <- paste(grep_term, collapse = "|")}
  module_term_genes <- list()
  for (mod in names(modDFs)) {
    modDF <- data.frame(modDFs[[mod]])
    modDF <- modDF[modDF$p.adjust < p_thres,]
    rowIndex <- grep(grep_term, modDF$Description,ignore.case=TRUE,value=FALSE)
    term_genes <- sort(unique(unlist(str_split(modDF[rowIndex,]$geneID, pattern = '/'))))
    module_term_genes[[mod]] <- term_genes
  }
  
  termGenesDF <- bind_rows(lapply(module_term_genes, as.data.frame, stringsAsFactors = FALSE), .id = "Module")
  colnames(termGenesDF)[2] <- 'Gene'
  
  # Sort by number of modules for facetting
  gene_counts <- sort(table(termGenesDF$Gene), decreasing = TRUE)
  gene_order <- names(gene_counts)
  gene_groups <- split(gene_order, sort(1:length(gene_order) %% split_column))
  num_mods <- max(gene_counts)
  
  if (is.null(color_vector)) {color_vector = module_colors}
  
  groupGenesDF <- termGenesDF[termGenesDF$Gene %in% gene_groups[[2]],]
  
  t <- ggplot(groupGenesDF) + 
    geom_bar(aes(y = fct_rev(fct_infreq(Gene)),
                 fill = groupGenesDF[,c('Module')])) + 
    labs(fill=xLab[1])  +
    scale_fill_manual(values = color_vector) +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,max(num_mods), 1), limits=c(0, max(num_mods))) +
    labs(y = ylab, x = paste(xLab, 'Count', sep = ' '))
  
  plt_list <- list()
  for (i in 1:length(gene_groups)) {
    
    groupGenesDF <- termGenesDF[termGenesDF$Gene %in% gene_groups[[i]],]
    
    plt_list[[i]] <- ggplot(groupGenesDF) + 
      geom_bar(aes(y = fct_rev(fct_infreq(Gene)),
                   fill = Module)) + 
      scale_fill_manual(values = color_vector,
                        name = xLab[1]) +
      theme_bw() +
      scale_x_continuous(breaks=seq(0,max(num_mods), 1), limits=c(0, max(num_mods))) +
      labs(y = ylab, x = paste(xLab, 'Count', sep = ' '))
    
  }
  
  plt <- ggarrange(plotlist = plt_list, nrow = 1, ncol = split_column, common.legend = TRUE, legend = "bottom")
  if (!is.null(title)) {plt <- annotate_figure(plt, top = title)}
  return(plt)
  
}

# Plot enrichment for specific term in each group
plot_term_by_group <- function(rez_list, grep_term, node_stats, modules = NULL, p_thres = 0.05,
                               title = 'Enrichment for RNA Splicing Related GO:BP Terms', includeAll = FALSE,
                               grouping = c("Module"), geneID = c("GeneID"), includeGrey = FALSE, plt_labels = NULL, color_vector) {
  
  if (is.null(color_vector)) {color_vector = module_colors}
  
  if (length(grep_term) > 1) {grep_term <- paste(grep_term, collapse = "|")}
  
  if (is.null(modules)) {
    modules = str_to_title(sort(names(rez_list)))
  } else {
    modules = str_to_title(modules)
  }
  if (includeGrey == FALSE) {modules = modules[modules != 'Grey']}
  
  # Inner function to get genes for given module
  get_genes <- function(node_stats, mod, grouping, geneID) {
    mod_nodes <- node_stats[node_stats[,grouping] == mod,]
    return(mod_nodes[,geneID])
  }
  
  all_mod_rez <- data.frame(rez_list[[modules[1]]])
  all_mod_rez <- all_mod_rez[FALSE,]
  # Fill dataframe using (up to) top N enriched terms for each module
  for (mod in modules) {
    rez <- data.frame(rez_list[[mod]])
    rez <- rez[rez$p.adjust < p_thres,]
    rowIndex <- grep(grep_term, rez$Description,ignore.case=TRUE,value=FALSE)
    rez <- rez[rowIndex,]
    if (dim(rez)[1] >= 1) {
      gene_hits <- as.numeric(str_split_fixed(rez$GeneRatio, '/', 2)[1,2])
      x <- data.frame(Description = rez$Description, p.adjust = rez$p.adjust, Count = rez$Count, 
                      GeneRatio = rez$Count/gene_hits, ID = rez$ID)
      x$Module <- rep(mod, dim(x)[1])
      x$logP <- -log10(x$p.adjust)
      x <- x[order(-x$p.adjust),]
      x$p.adjust <- format(x$p.adjust, digits = 2, scientific = TRUE)
      x$p.adjust <- paste(x$ID, '\n(', x$p.adjust, ')')
      x <- x[1,]
      x$Label <- paste(mod, '\n', as.character(length(unique(get_genes(node_stats, mod, grouping, geneID)))),
                       ' (', as.character(length(get_genes(node_stats, mod, grouping, geneID))), ')', sep = '')
      all_mod_rez <- rbind(all_mod_rez, x)
    } else {
      if (includeAll == TRUE) {
        gene_hits <- 0
        x <- data.frame(Description = NA, p.adjust = 1.0, Count = 0, 
                        GeneRatio = 0, ID = NA)
        x$Module <- rep(mod, dim(x)[1])
        x$logP <- -log10(1)
        x$p.adjust <- format(x$p.adjust, digits = 2, scientific = TRUE)
        x$p.adjust <- paste('(', x$p.adjust, ')')
        x$Label <- paste(mod, '\n', as.character(length(unique(get_genes(node_stats, mod, grouping, geneID)))),
                         ' (', as.character(length(get_genes(node_stats, mod, grouping, geneID))), ')', sep = '')
        all_mod_rez <- rbind(all_mod_rez, x)
        
      }
    }
    
  }
  
  # Fix ordering of enriched terms
  all_mod_rez$Description <- factor(all_mod_rez$Description, levels = unique(all_mod_rez$Description))
  all_mod_rez$Module <- factor(all_mod_rez$Module, levels = modules)
  if (is.null(plt_labels)) {
    all_mod_rez$Label <- factor(all_mod_rez$Label, levels = unique(all_mod_rez$Label))
    plt_labels = unique(all_mod_rez$Label)
  } 
  
  p <- ggplot(all_mod_rez, aes(x = logP, y=Module, fill = Module)) +
    geom_bar(stat='identity') +
    labs(y = paste(grouping[1], " Gene Count (Total Nodes)", sep = ''), fill = grouping[1], x = "-log10(p.adjust)") +
    xlim(c(0, round_any(max(all_mod_rez$logP) + 0.5, 0.5, ceiling))) +
    geom_text(aes(label = p.adjust), vjust = -0.5, position = position_dodge(0.9)) + 
    theme_light() + ggtitle(title) +
    theme(plot.title = element_text(size = 16), legend.key.size = unit(2, "line"),
          axis.text=element_text(size=12), axis.title=element_text(size=14)) +
    scale_fill_manual(values = color_vector,
                      name = grouping[1]) + 
    coord_flip() +
    scale_y_discrete(labels = plt_labels)
  
  return(p)
}


# Plot enrichment for specific term across module ratio threshold
plot_term_by_threshold <- function(rez_list, grep_term, node_stats, modules = NULL, p_thres = 0.05,
                                   title = 'Enrichment for RNA Splicing Related GO:BP Terms', 
                                   grouping = c("Module"), geneID = c("GeneID"), includeGrey = FALSE, plt_labels = NULL, color_vector) {
  
  if (is.null(color_vector)) {color_vector = module_colors}
  
  if (length(grep_term) > 1) {grep_term <- paste(grep_term, collapse = "|")}
  
  if (is.null(modules)) {
    modules = str_to_title(sort(names(rez_list)))
  } else {
    modules = str_to_title(modules)
  }
  if (includeGrey == FALSE) {modules = modules[modules != 'Grey']}
  
  # Inner function to get genes for given module
  get_genes <- function(node_stats, mod, grouping, geneID) {
    mod_nodes <- node_stats[node_stats[,grouping] == mod,]
    return(mod_nodes[,geneID])
  }
  
  all_mod_rez <- data.frame(rez_list[[modules[1]]])
  all_mod_rez <- all_mod_rez[FALSE,]
  # Fill dataframe using (up to) top N enriched terms for each module
  for (mod in modules) {
    rez <- data.frame(rez_list[[mod]])
    rez <- rez[rez$p.adjust < p_thres,]
    rowIndex <- grep(grep_term, rez$Description,ignore.case=TRUE,value=FALSE)
    rez <- rez[rowIndex,]
    print(dim(rez))
    if (dim(rez)[1] >= 1) {
      gene_hits <- as.numeric(str_split_fixed(rez$GeneRatio, '/', 2)[1,2])
      x <- data.frame(Description = rez$Description, p.adjust = rez$p.adjust, Count = rez$Count, 
                      GeneRatio = rez$Count/gene_hits, ID = rez$ID)
      x$Module <- rep(mod, dim(x)[1])
      x$logP <- -log10(x$p.adjust)
      x <- x[order(-x$p.adjust),]
      x$p.adjust <- format(x$p.adjust, digits = 2, scientific = TRUE)
      x$p.adjust <- paste(x$ID, '\n(', x$p.adjust, ')')
      x <- x[1,]
      x$Label <- paste(mod, '\n', as.character(length(unique(get_genes(node_stats, mod, grouping, geneID)))),
                       ' (', as.character(length(get_genes(node_stats, mod, grouping, geneID))), ')', sep = '')
      all_mod_rez <- rbind(all_mod_rez, x)
    } else {
      gene_hits <- 0
      x <- data.frame(Description = NA, p.adjust = 0, Count = 0, 
                      GeneRatio = 0, ID = NA)
      x$Module <- rep(mod, 1)
      
      x$logP <- 0
      x <- x[order(-x$p.adjust),]
      x$p.adjust <- format(x$p.adjust, digits = 2, scientific = TRUE)
      x$p.adjust <- ''
      x <- x[1,]
      x$Label <- paste(mod, '\n', as.character(length(unique(get_genes(node_stats, mod, grouping, geneID)))),
                       ' (', as.character(length(get_genes(node_stats, mod, grouping, geneID))), ')', sep = '')
      all_mod_rez <- rbind(all_mod_rez, x)
    }
  }
  
  
  # Fix ordering of enriched terms
  all_mod_rez$Description <- factor(all_mod_rez$Description, levels = unique(all_mod_rez$Description))
  all_mod_rez$Module <- factor(all_mod_rez$Module, levels = modules)
  if (is.null(plt_labels)) {
    all_mod_rez$Label <- factor(all_mod_rez$Label, levels = unique(all_mod_rez$Label))
    plt_labels = unique(all_mod_rez$Label)
  } 
  
  p <- ggplot(all_mod_rez, aes(x = Module, y=logP, group = 1)) +
    geom_line() + 
    geom_point() +
    labs(x = paste(grouping[1], " Gene Count (Total Nodes)", sep = ''),y = "-log10(p.adjust)") +
    ylim(c(0, round_any(max(all_mod_rez$logP) + 0.5, 0.5, ceiling))) +
    theme_light() + ggtitle(title) +
    theme(plot.title = element_text(size = 16), #legend.key.size = unit(2, "line"),
          axis.text=element_text(size=12), axis.title=element_text(size=14)) +
    scale_x_discrete(labels = plt_labels)
  
  return(p)
}



# Function to select top N enriched terms from each module specified and create a dotplot comparing modules
compare_module_enrichment <- function(rez_list, node_stats, modules = NULL, topN = NULL, p_thres = 0.05,
                                      title = 'Top Enriched Terms Across Modules', uniqueTerms = FALSE, maxShared = 1,
                                      grouping = c("Module"), geneID = c("GeneID"), minGenes = 1, includeGrey = FALSE, plt_labels = NULL) {
  if (is.null(modules)) {
    modules = str_to_title(sort(names(rez_list)))
  } else {
    modules = str_to_title(modules)
  }
  if (includeGrey == FALSE) {modules = modules[modules != 'Grey']}
  
  # Filter terms by gene count
  rez_list <- filter_enriced_terms_by_gene_count(rez_list[modules], minGenes = minGenes)
  
  if (uniqueTerms == TRUE) {rez_list <- unique_module_enrichment(rez_list[modules], maxShared = maxShared, p_thres = p_thres)$DFs}
  
  # Inner function to get genes for given module
  get_genes <- function(node_stats, mod, grouping, geneID) {
    mod_nodes <- node_stats[node_stats[,grouping] == mod,]
    return(mod_nodes[,geneID])
  }
  
  all_mod_rez <- data.frame(rez_list[[modules[1]]])
  all_mod_rez <- all_mod_rez[FALSE,]
  # Fill dataframe using (up to) top N enriched terms for each module
  for (mod in modules) {
    print(mod)
    rez <- data.frame(rez_list[[mod]])
    rez <- rez[rez$p.adjust < p_thres,]
    if (is.null(topN)) {
      n <- dim(rez)[1]
    } else {
      n <- min(dim(rez)[1], topN)
    }
    if (n >= 1) {
      gene_hits <- as.numeric(str_split_fixed(rez$GeneRatio, '/', 2)[1,2])
      x <- data.frame(Description = rez$Description, p.adjust = rez$p.adjust, Count = rez$Count, 
                      GeneRatio = rez$Count/gene_hits, ID = rez$ID)
      x$Module <- rep(mod, dim(x)[1])
      x <- x[order(x$p.adjust),]
      x <- x[1:n,]
      x <- x[order(-x$Count),]
      x$Label <- paste(mod, '\n', as.character(length(unique(get_genes(node_stats, mod, grouping, geneID)))),
                       ' (', as.character(length(get_genes(node_stats, mod, grouping, geneID))), ')', sep = '')
      all_mod_rez <- rbind(all_mod_rez, x)
    }
  }
  
  # Now add terms selected for each module where any module is enriched for such terms
  #all_mod_rez_init <- all_mod_rez
  #all_mod_rez <- all_mod_rez[FALSE,]
  for (mod in modules) {
    rez <- data.frame(rez_list[[mod]])
    rez <- rez[rez$p.adjust < p_thres,]
    if (is.null(topN)) {
      n <- dim(rez)[1]
    } else {
      n <- min(dim(rez)[1], topN)
    }
    if (n >= 1) {
      gene_hits <- as.numeric(str_split_fixed(rez$GeneRatio, '/', 2)[1,2])
      x <- data.frame(Description = rez$Description, p.adjust = rez$p.adjust, Count = rez$Count, 
                      GeneRatio = rez$Count/gene_hits, ID = rez$ID)
      x$Module <- rep(mod, dim(x)[1])
      x$Label <- paste(mod, '\n', as.character(length(unique(get_genes(node_stats, mod, grouping, geneID)))),
                       ' (', as.character(length(get_genes(node_stats, mod, grouping, geneID))), ')', sep = '')
      x <- x[x$ID %in% all_mod_rez$ID,]
      all_mod_rez <- rbind(all_mod_rez, x)
    }
  }
  
  # Fix ordering of enriched terms
  all_mod_rez$Description <- factor(all_mod_rez$Description, levels = unique(all_mod_rez$Description))
  all_mod_rez$Module <- factor(all_mod_rez$Module, levels = modules)
  if (is.null(plt_labels)) {
    all_mod_rez$Label <- factor(all_mod_rez$Label, levels = unique(all_mod_rez$Label))
    plt_labels = unique(all_mod_rez$Label)
  } 
  
  p <- ggplot(all_mod_rez, aes(x = Description, y=Module)) +
    geom_point(aes(colour=p.adjust, size=Count)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(y = paste(grouping[1], " Gene Count (Total Nodes)", sep = ''), fill = "p.adjust", x = "") +
    theme_light() + ggtitle(title) +
    theme(plot.title = element_text(size = 16), legend.key.size = unit(2, "line"),
          axis.text=element_text(size=12), axis.title=element_text(size=14)) +
    #coord_fixed(ratio = 0.5) + 
    coord_flip() +
    scale_x_discrete(limits = rev(levels(all_mod_rez$Description))) + 
    scale_y_discrete(labels = plt_labels)
  
  return(p)
}

# Functiont to plot number of enriched terms per module
plot_module_enrichment_counts <- function(modCounts, node_stats, y_lab, plt_title, grouping = c('Module')) {
  
  mod_labels <- c()
  for (mod in modCounts[,1]) {
    mod_nodes <- node_stats[node_stats[,grouping] == mod,]
    num_nodes <- length(rownames(mod_nodes))
    num_genes <- length(unique(mod_nodes$GeneID))
    mod_labels <- c(mod_labels, paste(mod, '\n', as.character(num_genes), ' (', as.character(num_nodes), ')', sep = ''))
  }
  
  modCounts <- cbind(modCounts, mod_labels)
  
  y_max <- round_any(max(modCounts[,2]), 25, f = ceiling) #max(modCounts[,2])
  
  mc_plot <- ggplot(modCounts, aes(x=reorder(mod_labels,-modCounts[,2]), y = modCounts[,2], fill = Module)) +
    geom_bar(stat="identity") + scale_y_continuous(breaks=seq(0,y_max,5)) + 
    scale_fill_manual(values = module_colors) +
    geom_text(stat='identity', aes(label=modCounts[,2]), position=position_dodge(0.5), vjust=-0.55) + 
    ggtitle(plt_title) +
    xlab("Module Gene Count (Total Nodes)") + ylab(y_lab) +
    theme_bw() +
    theme(text = element_text(size=16))#, axis.text.x=element_blank(), axis.ticks.x=element_blank())
  return(mc_plot)
}
