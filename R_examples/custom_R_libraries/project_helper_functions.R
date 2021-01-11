################################################################
#
#   Co-splicing Network Inference R Library
#
#     Collection of functions for processing and analyzing 
#     splicing data for de novo network inference
#
#     Also sources all plotting libraries used for project
#
#################################################################


set.seed(101)  

# Load required R libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(reshape2)
library(scales)
library(yarrr)
library(WGCNA)
library(stringr)  
library(RColorBrewer)
library(yarrr)
library(colorspace)
library(VennDiagram)
library(stringr)
library(ggpubr)
library(factoextra)
library(ggfortify)
library(gtools)
library(data.table)
library(doParallel)
library(MASS)
library(NetRep)
library(tibble)
library(sva)
library(limma)
library(impute)
library(igraph)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(biomaRt)
library(randnet)
library(GO.db)
library(parallel)
library(gridExtra)
library(ComplexHeatmap)

# Load project plotting functions
source("/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/shared_Rscripts/splicing_EDA_plotting_functions.R")
source("/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/shared_Rscripts/net_inference_plotting_functions.R")
source("/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/shared_Rscripts/module_quality_plotting_functions.R")
source("/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/shared_Rscripts/module_analysis_plotting_functions.R")
source("/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/shared_Rscripts/enrichment_plotting_functions.R")



# Function to perform logit transformation on PSI matrix during batch correction
to_logit <- function(dat){
  eps <- 0.0000000001 # Convert 0s with constant to ensure no -Inf values
  dat[dat == 0] <- eps
  return(logit(dat)) # use logit function from gtools
}

# For function below, search redundant LSVs to see if it is in PSI matrix
search_redundant_lsvs <- function (redundant_lsv, exprMat) {
  if (redundant_lsv != 'NaN') {
    # Check if there's multiple pairs
    # If not, see if it's in exprMat and return it
    if (!grepl('|', redundant_lsv, fixed = TRUE)) {
      if (redundant_lsv %in% colnames(exprMat)) {
        return(as.character(redundant_lsv))
      }
    } 
    # If so, we need to parse the string until we find one that is present
    else {
      split_lsvs <- strsplit(as.character(redundant_lsv), '|', fixed = TRUE)[[1]]
      for (i in 1:length(split_lsvs)) {
        lsv <- split_lsvs[i]
        if (lsv %in% colnames(exprMat)) {
          return(as.character(lsv))
        }
      }
    }
  }  
  return(NaN)
}


# Function to compute splice variant regions from non-redundant LSV sets
make_SVRs <- function(exprMat, lsvDict) {
  
  # Select LSVs that were initially labeled as non-redundant (python script)
  #   --> Initially selected based on PSI sample variance from MAJIQ
  nonRedundant <- rownames(lsvDict[lsvDict$IS_REDUNDANT == 'keep',])
  
  # If any non-redundant LSVs are absent due to coverage, check if there is a 
  # redundant mate that is present in PSI matrix
  lsv_set <- c()
  for (i in 1:length(nonRedundant)) {
    lsv <- nonRedundant[i]
    # If it's present, use it
    if (lsv %in% colnames(exprMat)) {
      lsv_set <- c(lsv_set, lsv) 
    }
    # If not, see if it exists in redundant pairing
    else {
      redundant_lsv <- lsvDict[lsv,]$REDUNDANT_PAIR
      search_lsv <- search_redundant_lsvs(redundant_lsv, exprMat)
      if (!is.na(search_lsv)) {
        lsv_set <- c(lsv_set, search_lsv)
      }
    }
  }

  # Create new PSI matrix of selected LSVs
  #  --> Use unique in case any duplicate LSVs were selected as non-redundant
  useExpr <- exprMat[,unique(lsv_set)]
  
  # Create subset of lsv dictionary with row ordering matching PSI matrix
  lsvDict_set <- lsvDict[colnames(useExpr),]
  
  # Use WGCNA's eigengene function for convenient calculation of PC1
  clusterExpr <- moduleEigengenes(useExpr, as.character(lsvDict_set$SVR))$eigengenes
  names(clusterExpr) <- substring(names(clusterExpr), 3) # Fix this ME crap
  
  return(clusterExpr)
}


# Function to compute module expression using original LSVs instead of the SVRs (doesn't do much)
# Need original lsv matrix, node_stats for module assignments, and lsv dictionary
lsv_based_mod_expression <- function(exprMat, node_stats, lsv_dict) {
  
  lsv_color_assign <- c()
  for(lsv in colnames(exprMat)) {
    svr <- as.character(lsv_dict[lsv,]$SVR)
    mod <- node_stats[svr,]$Module
    lsv_color_assign <- c(lsv_color_assign, mod)
  }
  
  # Calculate module expression for AML
  modExpr.lsvs = moduleEigengenes(exprMat, colors = lsv_color_assign, excludeGrey = TRUE)$eigengenes
  names(modExpr.lsvs) <- str_to_title(substring(names(modExpr.lsvs), 3))
  
  return(modExpr.lsvs)
  
}


# Simple function for taking a list of expression matrices
#  and joining them into single expression matrix
multi_to_single_expr <- function(exprList, isMulti = FALSE) {
  if (isMulti == TRUE) {
    singleExpr <- exprList[[1]]$data
    for (i in 2:length(exprList)) {
      singleExpr <- rbind(singleExpr, exprList[[i]]$data)  
    }
  } else {
    singleExpr <- exprList[[1]]
    for (i in 2:length(exprList)) {
      singleExpr <- rbind(singleExpr, exprList[[i]])  
    }
  }
  return(singleExpr)
}

## Function to remove "ME" substring on module colors
#   for multi-Expr lists placed on by WGCNA eigengene function
# Note this will just return single matrices ($data) containing 
# sample module valuesfor each list item
multi_to_colors <- function(exprList) {
  newExprList = list()
  for (i in 1:length(exprList)) {
    newExprList[[i]] <- exprList[[i]]$data
    names(newExprList[[i]]) <- str_to_title(substring(names(newExprList[[i]]), 3))
  }
  if (!is.null(names(exprList))) {
    names(newExprList) <- names(exprList)
  }
  return(newExprList)
}


# Use LSV dictionary to convert SVR or list of SVRs to their gene of origin
#  Can convert to EnsemblID or gene name
svr_to_gene <- function(svrs, lsv_dict, gene_name = FALSE) {
  if (gene_name == TRUE) {i <- 3} else {i <- 2}
  gene_list <- c()
  for (svr in svrs) {
    gene_list <- c(gene_list, as.character(unique(lsv_dict[lsv_dict[,7] == svr,i])))
  }
  return(gene_list)
}



# Function for computing node connectivity and module membership values
#   Nodes need to be in same order for adj and nodeExpr
#    If list, adj and nodeExpr must have names
get_node_stats <- function(adj, mod_assignments, nodeExpr, isList = FALSE, convert_SVRs = FALSE, 
                           lsv_dict = NULL, getEntrez = TRUE) {
  # Inner function for running computations
  compute_stats <- function(adj, mod_assignments, geneIDs, nodeExpr, modExpr, modMem.across = NULL) {
    # Compute node connectivity values
    node_connectivities <- intramodularConnectivity(adj, mod_assignments, scaleByMax = FALSE)
    node_connectivities$kWithinScaled <- intramodularConnectivity(adj, mod_assignments, scaleByMax = TRUE)$kWithin
    
    # Add module assignments and geneIDs to dataframe
    node_connectivities$Module <- str_to_title(mod_assignments)
    node_connectivities$GeneID <- geneIDs
    # Compute module membership values
    mod_membership <- signedKME(nodeExpr, modExpr)
    colnames(mod_membership) <- colnames(modExpr)
    # Select module membership value for module to which node was assigned
    node_connectivities$kME_Within <- 0.0
    for (node in rownames(node_connectivities)) {
      mod <- node_connectivities[node,c("Module")]
      kME_vals <- mod_membership[node,]
      node_connectivities[node,]$kME_Within <- kME_vals[,mod]
    }
    # Calculate module membership using module eigengenes calculated across dataset
    if (!(is.null(modMem.across))) {
      node_connectivities$kME_Across <- 0.0
      for (node in rownames(node_connectivities)) {
        mod <- node_connectivities[node,c("Module")]
        kME_vals <- modMem.across[node,]
        node_connectivities[node,]$kME_Across <- kME_vals[,mod]
      }
    } else {
      node_connectivities$kME_Across <- node_connectivities$kME_Within
    }
    return(node_connectivities[,c(7,6,8,9,1,2,5,3,4)])
  }
  if (isList == FALSE) {
    if (convert_SVRs == TRUE) {
      geneIDs <- svr_to_gene(rownames(adj), lsv_dict, gene_name = FALSE)
      geneNames <- svr_to_gene(rownames(adj), lsv_dict, gene_name = TRUE)
    } else {
      geneIDs <- rownames(adj)
      geneNames <- rownames(adj)
    }
    # If not list, use moduleEigengenes function
    modExpr <- moduleEigengenes(nodeExpr, mod_assignments, excludeGrey = FALSE)$eigengenes
    names(modExpr) <- str_to_title(substring(names(modExpr), 3))
    # Compute stats and return node stats dataframe
    node_connectivities <- compute_stats(adj, mod_assignments, geneIDs, nodeExpr, modExpr, modMem.across = NULL)
    node_connectivities$SVR <- rownames(node_connectivities)
    node_connectivities$GeneName <- geneNames
    if (getEntrez == TRUE) {
      entIDs <- ensembl_to_entrez(geneIDs)
      node_connectivities$ENTREZ <- entIDs
    }
    node_connectivities <- node_connectivities[,c("SVR", "GeneID", "GeneName", "ENTREZ","Module","kME_Within","kME_Across","kTotal",
                                                  "kWithin","kWithinScaled","kOut","kDiff")]
  } else {
    if (convert_SVRs == TRUE) {
      geneIDs <- svr_to_gene(rownames(adj[[1]]), lsv_dict, gene_name = FALSE)
      geneNames <- svr_to_gene(rownames(adj[[1]]), lsv_dict, gene_name = TRUE)
    } else {
      geneIDs <- rownames(adj[[1]])
      geneNames <- rownames(adj[[1]])
    }
    if (getEntrez == TRUE) {entIDs <- ensembl_to_entrez(geneIDs)}   
    # If list, use multiSetMEs function to compute module eigengenes (group specific)
    modExpr <- multi_to_colors(multiSetMEs(nodeExpr, universalColors = mod_assignments, excludeGrey = FALSE))
    # Compute module eigengenes across all samples modExpr.across (MEs across groups)
    modExpr.across <- moduleEigengenes(multi_to_single_expr(nodeExpr, isMulti = TRUE), mod_assignments, excludeGrey = FALSE)$eigengenes
    names(modExpr.across) <- str_to_title(substring(names(modExpr.across), 3))
    modMem.across <- signedKME(multi_to_single_expr(nodeExpr, isMulti = TRUE), modExpr.across)
    colnames(modMem.across) <- colnames(modExpr.across)
    node_connectivities <- list()
    for (n in names(adj)) {
      node_connectivities[[n]] <- compute_stats(adj[[n]], mod_assignments, geneIDs, nodeExpr[[n]]$data, modExpr[[n]], modMem.across = modMem.across)
      if (getEntrez == TRUE) {node_connectivities[[n]]$ENTREZ <- entIDs}
      node_connectivities[[n]]$SVR <- rownames(node_connectivities[[n]])
      node_connectivities[[n]]$GeneName <- geneNames
      node_connectivities[[n]] <- node_connectivities[[n]][,c("SVR", "GeneID", "GeneName", "ENTREZ","Module","kME_Within","kME_Across","kTotal",
                                                              "kWithin","kWithinScaled","kOut","kDiff")]
    }
  }
  return(node_connectivities)
}


# Helper function to take a list of node_stats dataframes and create a single dataframe containing
#  the stat of interest for nodes of each group (used for multi-network analysis)
group_node_stats <- function(node_stats_list, statType) {
  node_matrix <- matrix(nrow = dim(node_stats_list[[1]][1]),ncol = length(node_stats_list))
  for (i in 1:length(node_stats_list)) { node_matrix[,i] <- node_stats_list[[i]][[statType]] }
  node_matrix <- data.frame(node_matrix)
  rownames(node_matrix) <- rownames(node_stats_list[[1]])
  colnames(node_matrix) <- names(node_stats_list)
  # Add module assignments
  node_matrix$Module <- node_stats_list[[1]]$Module
  return(node_matrix)
}


# Function for taking node_stats dataframe (or list of dataframes) and subsetting nodes from each module 
#  based on a given paramemeter (e.g. Intra-modular node connectivity)
#   --> Default is to simply remove duplicate genes from each module for a given dataframe 
subset_module_nodes <- function(node_stats, num_nodes = NULL, modules = NULL, grouping = c("Module"), geneID = c("GeneID"),
                                uniqueModGenes = TRUE, stat_sort = NULL, useAbsolute = FALSE, reverseSort = FALSE, isList = FALSE, 
                                includeGrey = FALSE, keepShared = TRUE) {
  
  # Note that useAbsolute will order by absolute value, but values returned are not absolute value (contain negatives)
  if ((!is.null(stat_sort)) && (stat_sort == 'kME_Within' || stat_sort == 'kME_Across')) {useAbsolute = TRUE}
  
  # Inner function to select top nodes based on provided statistic and module set
  select_nodes <- function(node_stats, num_nodes, modules, grouping, stat_sort, useAbsolute, reverseSort) {
    node_list <- c() # Store nodes (rownames for each module)
    if (useAbsolute == TRUE) {
      if (reverseSort == TRUE) {
        sortDF <- node_stats[order(abs(node_stats[,stat_sort])),]
      } else { 
        sortDF <- node_stats[order(-abs(node_stats[,stat_sort])),]
      }
    } else {
      if (reverseSort == TRUE) {
        sortDF <- node_stats[order(node_stats[,stat_sort]),] # Order by given stat
      } else {  
        sortDF <- node_stats[order(-node_stats[,stat_sort]),] # Order by given stat
      }
    }
    for (mod in modules) {
      modNodes <- sortDF[sortDF[,grouping] == mod,]
      # Select top N nodes
      if (num_nodes < 1) {nNodes = ceiling((dim(modNodes)[1])*num_nodes)} else {nNodes = num_nodes}
      node_select <- modNodes[1:nNodes,]
      node_list <- c(node_list, rownames(node_select))
    }
    return(node_list)
  }
  
  # Inner function to remove genes appearing in more than one module
  removeShared <- function(newDF, grouping, gID) {
    modCounts <- aggregate(newDF[,grouping] ~ newDF[,gID], newDF, function(x) length(unique(x)))
    uniqueGenes <- modCounts[modCounts[,2] == 1,][,1]
    return(newDF[newDF[,gID] %in% uniqueGenes,])
  }
  if (isList == TRUE) {
    # If list, parse each dataframe and select nodes from each module
    # Store in list containing all nodes selected from each dataset 
    if (!is.null(modules)) {modules = str_to_title(modules)} else {modules = unique(str_to_title(node_stats[[1]][,grouping]))}
    if (includeGrey == FALSE) {modules <- modules[!(modules=="Grey")]}
    # Store in here
    select_nodes_list <- list()
    for (n in names(node_stats)) {
      sortDF <- node_stats[[n]][sort(rownames(node_stats[[n]])),] # Sort by rowname for consistency
      sortDF <- sortDF[sortDF[,grouping] %in% modules,] # Select nodes from specified modules
      # If no statistic (connectivity, kME), remove duplicate genes from each module
      if (is.null(stat_sort)) {
        select_nodes_list[[n]] <- rownames(sortDF)  # Store nodes to keep
      } else { # Else select top nodes for each module and return those nodes for current dataset
        select_nodes_list[[n]] <- select_nodes(sortDF, num_nodes, modules, grouping, stat_sort, useAbsolute, reverseSort)
      }
    }
    # Unpack list
    select_nodes_list <- unique(unlist(select_nodes_list, use.names = FALSE))
    # Now fill each node_stats dataframe using selected nodes
    new_node_stats <- list()  
    for (n in names(node_stats)) {
      newDF <- node_stats[[n]][select_nodes_list,]
      # Remove duplicate genes from each module
      if (uniqueModGenes == TRUE) {newDF <- newDF[!duplicated(newDF[,c(geneID, grouping)]),]}
      # Only select genes that are module-specific (not in more than one module)
      if (keepShared == FALSE) {newDF <- removeShared(newDF, grouping, geneID)}
      new_node_stats[[n]] <- newDF
    }
    
  } else {
    # If not list, just select nodes from each module based on stat
    if (!is.null(modules)) {modules = str_to_title(modules)} else {modules = unique(str_to_title(node_stats[,grouping]))}
    if (includeGrey == FALSE) {modules <- modules[!(modules=="Grey")]}
    sortDF <- node_stats[sort(rownames(node_stats)),] # Sort by rowname for consistency
    sortDF <- sortDF[sortDF[,grouping] %in% modules,] # Select nodes from specified modules
    # If no statistic (connectivity, kME), remove duplicate genes from each module
    if (is.null(stat_sort)) {
      select_nodes_list <- rownames(sortDF)
    } else { # Else select top nodes for each module and return those nodes for current dataset
      select_nodes_list <- select_nodes(sortDF, num_nodes, modules, grouping, stat_sort, useAbsolute, reverseSort)
    }
    newDF <- node_stats[select_nodes_list,]
    # Remove duplicate genes from each module
    if (uniqueModGenes == TRUE) {newDF <- newDF[!duplicated(newDF[,c(geneID, grouping)]),]}
    # Only select genes that are module-specific (not in more than one module)
    if (keepShared == FALSE) {newDF <- removeShared(newDF, grouping, geneID)}
    new_node_stats <- newDF
  }
  
  return(new_node_stats)
  
}



subset_nodes_by_group <- function(node_stats, groups_name, num_nodes = NULL, stat_sort = NULL, modules = NULL, net_groups = NULL, group_column = c("Module"), 
                                  geneID = c("GeneID"), uniqueModGenes = TRUE, useAbsolute = FALSE, reverseSort = FALSE, includeGrey = FALSE,
                                  keepShared = TRUE) {
 
  # Note that useAbsolute will order by absolute, but values returned are not absolute value (contain negatives)
  if ((!is.null(stat_sort)) && (stat_sort == 'kME_Within' || stat_sort == 'kME_Across')) {useAbsolute = TRUE}
  
  if (is.null(net_groups)) { net_groups <- names(node_stats) }
  
  # Get top nodes from each group
  node_set <- node_stats[[1]]
  node_set[groups_name] <- NA
  node_set <- node_set[FALSE,] # Initial empty node dataframe
  for (net in net_groups) {
    group_nodes <- node_stats[[net]]
    group_nodes[groups_name] <- rep(net, dim(group_nodes)[1])
    group_nodes <- subset_module_nodes(group_nodes, num_nodes = num_nodes, stat_sort = stat_sort, modules = modules, grouping = group_column, 
                                       geneID = c(geneID), uniqueModGenes = uniqueModGenes, useAbsolute = useAbsolute, reverseSort = reverseSort, 
                                       includeGrey = includeGrey)
    
    node_set <- rbind(node_set, group_nodes)
  }

  return(node_set)
    
}

# Function for computing rank product values for node stats given a list of node_stat dataframes
compute_node_ranks <- function(node_stats_list, select_stats = c('kME_Within','kME_Across','kTotal','kWithin','kWithinScaled','kOut','kDiff'), useAbsolute = TRUE,
                               makeAbsolute = c('kME_Within', 'kME_Across'), group_column = c('Module'), modules = NULL) {
  
  # Use if not selecting by any groups (e.g. module) group_column == None
  if (group_column == 'None') {
    tmp_stats_list <- list()
    group_column <- c('Hold')
    for (n in 1:length(node_stats_list)) {
      group_stats <- node_stats_list[[n]]
      group_stats[,group_column] <- 'Holder'
      tmp_stats_list[[n]] <- group_stats
    }
    node_stats_list <- tmp_stats_list
  }
  
  
  if (is.null(modules)) {modules = sort(unique(node_stats_list[[1]][,group_column]))}
  new_node_stats <- node_stats_list[[1]][FALSE,]
  
  for (mod in modules) {
    module_nodes <- node_stats_list[[1]][node_stats_list[[1]][,group_column] == mod,] # Get all nodes for module
    group_module_node_ranks <- list() # Store module node ranks for each group
    for (n in 1:length(node_stats_list)) {
      group_module_nodes <- node_stats_list[[n]][node_stats_list[[n]][,group_column] == mod,] # Nodes for current group of current module
      module_group_ranks <- group_module_nodes[,select_stats] %>%
        mutate_each(funs(dense_rank(-.))) %>% data.matrix() %>% apply(2, FUN=as.numeric)
      group_module_node_ranks[[n]] <- module_group_ranks
    }    
    
    # Now compute rank product across groups for current module
    module_node_rank_product <- Reduce("*", group_module_node_ranks) %>% data.frame() %>%
      mutate_each(funs(dense_rank(.))) %>% as.matrix() %>% apply(2, FUN=as.numeric)
    
    # Add module node ranks to new node_stats dataframe    
    new_module_node_stats <- module_nodes
    new_module_node_stats[,select_stats] <- module_node_rank_product[,select_stats]
    new_node_stats <- rbind(new_node_stats, new_module_node_stats)
  }
  if (group_column[1] == 'Hold') {new_node_stats <- new_node_stats[,-which(colnames(new_node_stats) == group_column)]}
  return(new_node_stats)
}


#################################################
#
#   Module Enrichment Functions
#
#################################################

# Function that takes a vector of module assisgnments and a vector of nodes (equal size) and returns
#  the Entrez IDs for all genes assigned to each module
module_to_entrez <- function(modules, genes, include_grey = FALSE, convert_SVRs = FALSE, lsv_dict = NULL) {
  if (include_grey == TRUE) {
    mod_names <- sort(unique(modules))
  } else {
    mod_names <- sort(unique(modules[modules != 'grey']))
  }
  
  if (convert_SVRs == TRUE) { genes <- svr_to_gene(genes, lsv_dict)  }
  
  gene_to_module <- data.frame(Gene = genes, Module = modules)
  
  #Get unique Entrez IDs for all genes in each module
  entrezIDs.modules <- c()
  mod_set <- c()
  for (mod in mod_names) {
    ensIDs <- gene_to_module[gene_to_module$Module == mod,]$Gene
    ensg_to_entrz <- bitr(ensIDs, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
    entIDs <- c()
    for(ensg in ensIDs) {
      entIDs <- c(entIDs, ensg_to_entrz[ensg_to_entrz$ENSEMBL == ensg,2][1]) # Select unique entrez ID incase of multiples for each Ensemble ID
    }
    
    entIDs <- unique(entIDs[!(is.na(entIDs))])
    entrezIDs.modules <- c(entrezIDs.modules, entIDs)
    mod_set <- c(mod_set, rep(mod, length(entIDs)))
  }
  mod2entrez <- data.frame(Module = mod_set, ENTREZ = entrezIDs.modules)
  return(mod2entrez)
}


# Function for getting ENTREZ ID from Ensembl Gene ID for enrichment analysis
ensembl_to_entrez <- function(geneIDs) {
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  rez <- getBM(attributes = c('ensembl_gene_id','entrezgene_id'),
               filters = 'ensembl_gene_id', 
               values = geneIDs,mart = mart)
  rez2 <- rez[!(is.na(rez$entrezgene_id)),]
  rez3 <- rez2[!duplicated(rez2[ , c("ensembl_gene_id")]),]
  entIDs <- c()
  for (gID in geneIDs) {
    if (gID %in% rez3$ensembl_gene_id) {
      eID <- rez3[rez3$ensembl_gene_id == gID,]$entrezgene_id
    } else {
      eID <- NA
    }
    entIDs <- c(entIDs, eID)
  }
  return(entIDs)
}
  
# Function for performing enrichment analysis of modules using clusterProfiler and ReactomePA
module_enrichment <- function(mod2entrez, enrich_function, modules = NULL, grouping = c("Module"), adj_method = "fdr", p_thres = 0.1, 
                              goType = "BP", includeGrey = FALSE, background = NULL) {
  if (!is.null(modules)) {modules = str_to_title(modules)} else {modules = unique(str_to_title(mod2entrez[,grouping]))}
  if (includeGrey == FALSE) {modules = modules[modules != 'Grey']}
  
  numEnriched <- c()
  dp_list <- list()
  df_list <- list()
  for (mod in modules) {
    print(mod)
    modEntrezIDs <- unique(mod2entrez[mod2entrez[,grouping] == mod,]$ENTREZ)
    if (enrich_function == 'enrichPathway') {
      modRez <- enrichPathway(gene=modEntrezIDs,pvalueCutoff=p_thres, readable=T, pAdjustMethod = adj_method, universe = background)
      plot_title <- paste(str_to_title(mod), ' Module: Top Enriched Reactome Pathways (FDR < ',p_thres, ')',sep = '')
    } else if (enrich_function == 'enrichGO') {
      modRez <- enrichGO(gene = modEntrezIDs, keyType = "ENTREZID", OrgDb = org.Hs.eg.db,
                         ont = goType, pAdjustMethod = adj_method, pvalueCutoff = p_thres, readable = TRUE, universe = background)
      #universe = net_entrezIDs
      if (goType == "BP") {
        goName <- 'Biological Process'
      } else if (goType == "MF") {
        goName <- 'Molecular Function'
      }
      plot_title <- paste(str_to_title(mod), ' Module: Top Enriched GO (',goName,') Terms (FDR < ',p_thres, ')',sep = '')
    }
    
    if (is.null(modRez)) {
      modDP <- NULL
      modCount <- 0
    } else {  
      modDP <- dotplot(modRez, showCategory=10, font.size = 16, title = plot_title)
      modDP <- modDP + theme(plot.title = element_text(size = 18, face = "bold"))
      modCount <- dim(modRez)[1]
    }
    dp_list[[mod]] <- modDP 
    df_list[[mod]] <- modRez
    numEnriched <- c(numEnriched, modCount)
  }
  return_list <- list()
  return_list[['modCounts']] <- data.frame(Module = modules, NumEnriched = numEnriched)
  return_list[['DPs']] <- dp_list
  return_list[['DFs']] <- df_list
  return_list[['p.adjust_thres']] <- p_thres
  return(return_list)
}


# Function that takes a list of result objects from clusterProfiler and creates a
#  new list of data frames containing unique enrichment terms for each module
#   Can specify a maximum shared across module set (e.g. 2 modules instead of 1)
unique_module_enrichment <- function(rez_list, modules = NULL, p_thres = 0.05, maxShared = 1) {
  if (is.null(modules)) {modules = names(rez_list)}
  
  all_rez <- data.frame(rez_list[[modules[1]]])
  for (mod in modules[2:length(modules)]) {
    all_rez <- rbind(all_rez, data.frame(rez_list[[mod]]))
  }
  all_rez <- all_rez[all_rez[,c('p.adjust')] < p_thres,]
  
  # Get terms found in only one module
  unique_terms <- names(table(all_rez$ID)[table(all_rez$ID) <= maxShared])
  
  unique_rez_list <- list()
  numEnriched <- c()
  for (mod in modules) {
    mod_rez <- data.frame(rez_list[[mod]])
    mod_rez <- mod_rez[mod_rez$ID %in% unique_terms,]
    unique_rez_list[[mod]] <- mod_rez    
    modCount <- dim(mod_rez)[1]
    numEnriched <- c(numEnriched, modCount)
  }
  return_list <- list()
  return_list[['modCounts']] <- data.frame(Module = modules, NumEnriched = numEnriched)
  return_list[['DFs']] <- unique_rez_list
  return_list[['p.adjust_thres']] <- p_thres
  return(return_list)
}

# Used to removed enriched terms having less than minimum gene count
filter_enriced_terms_by_gene_count <- function(rez_list, minGenes = 5, modules = NULL) {
  if (is.null(modules)) {modules = names(rez_list)}
  filtered_rez_list <- list()
  for (mod in modules) {
    rez <- data.frame(rez_list[[mod]])
    filtered_rez_list[[mod]] <- rez[rez[,c('Count')] >= minGenes,]
  }
  return(filtered_rez_list)
}

# Function that takes a list of result dataframes from clusterProfiler and filters them based on criteria
filter_and_sort_enrichment_results <- function(rez_list, modules = NULL, topN = 5, p_thres = 0.05, select_by = c('p.adjust'), 
                                               display_by = c('Count'), reverse_select = FALSE, reverse_display = TRUE) {
  
  if (is.null(modules)) {modules <- names(rez_list)}
  filtered_rez_list <- list()
  
  for (mod in modules) {
    rez <- data.frame(rez_list[[mod]])
    rez <- rez[rez$p.adjust <= p_thres,]
    if (reverse_select == TRUE) {
      rez <- rez[order(-rez[,select_by]),]
    } else {
      rez <- rez[order(rez[,select_by]),]
    } 
    n <- min(topN, dim(rez)[1])
    rez <- rez[1:n,]
    if (reverse_display == TRUE) {
      rez <- rez[order(-rez[,display_by]),]
    } else {
      rez <- rez[order(rez[,display_by]),]
    }
    
    filtered_rez_list[[mod]] <- rez
  }
  
  return(filtered_rez_list)
}


# Function to remove redundant GO terms based on similarity in GO tree
filter_redundant_GO_terms <- function(rez_list, cutoff = 0.7) {
  simplified_rez_list <- list()
  for (t in names(rez_list)) {
    rez <- rez_list[[t]]
    simplified_rez_list[[t]] <- clusterProfiler::simplify(rez, cutoff = cutoff, measure = "Wang")
  }
  return(simplified_rez_list)
}

########################################################################################
#
#     Functions for looking at distribution of edge weights in various forms
#
########################################################################################

# Return a matrix of unique edge weights given an adjacency matrix
get_edge_weights <- function(x) {
  x[lower.tri(x)] <- NA
  diag(x) <- NA
  pw_cors <- as.data.frame(as.table(x))
  pw_cors <- pw_cors[!is.na(pw_cors$Freq),]
  return(pw_cors)
}

get_edge_types <- function(pw_cors){
  pw_cors[,4] <- str_split_fixed(pw_cors$Var1, ":", 3)[,1]
  pw_cors[,5] <- str_split_fixed(pw_cors$Var2, ":", 3)[,1]
  same_genes <- pw_cors[pw_cors$Var1!=pw_cors$Var2 & pw_cors$V4==pw_cors$V5,]
  diff_genes <- pw_cors[pw_cors$Var1!=pw_cors$Var2 & pw_cors$V4!=pw_cors$V5,]
  same_genes$EdgeType <- 'Within Genes'
  diff_genes$EdgeType <- 'Between Genes'
  merged_cors <- rbind(same_genes[,c(1,2,3,6)], diff_genes[,c(1,2,3,6)])
  return(merged_cors)
}

get_edge_types_SVRs <- function(pw_cors){
  pw_cors[,4] <- str_split_fixed(pw_cors$Var1, "_", 3)[,1]
  pw_cors[,5] <- str_split_fixed(pw_cors$Var2, "_", 3)[,1]
  same_genes <- pw_cors[pw_cors$Var1!=pw_cors$Var2 & pw_cors$V4==pw_cors$V5,]
  diff_genes <- pw_cors[pw_cors$Var1!=pw_cors$Var2 & pw_cors$V4!=pw_cors$V5,]
  same_genes$EdgeType <- 'Within Genes'
  diff_genes$EdgeType <- 'Between Genes'
  merged_cors <- rbind(same_genes[,c(1,2,3,6)], diff_genes[,c(1,2,3,6)])
  return(merged_cors)
}

get_same_gene_edge_types <- function(same_genes, lsv_to_clust) {
  pw_cors <- same_genes[,c(1,2,3)]
  pw_cors[,4] <- lsv_to_clust[as.character(same_genes[,1]),7]
  pw_cors[,5] <- lsv_to_clust[as.character(same_genes[,2]),7]
  same_overlap <- pw_cors[pw_cors$V4 == pw_cors$V5,]
  diff_overlap <- pw_cors[pw_cors$V4 != pw_cors$V5,]
  same_overlap$EdgeType <- 'Within Overlaps'
  diff_overlap$EdgeType <- 'Between Overlaps'
  merged_cors <- rbind(same_overlap[,c(1,2,3,6)], diff_overlap[,c(1,2,3,6)])
  return(merged_cors)
}




#########################################################
#
#     Module quality functions
#
#########################################################

# Re-retrieve module colors from a single NetRep run
#   Not for bootstrap runs (function is in local library code)
reLabel_mods <- function(pb){
  colorOrder = c("grey", standardColors(50))
  new_labs <- colorOrder[as.numeric(c(rownames(pb[[1]])))+1]
  new_labs <- str_to_title(new_labs)
  for (i in 1:length(pb)){
    rownames(pb[[1]]) <- new_labs
    rownames(pb[[2]]) <- new_labs
    rownames(pb[[3]]) <- new_labs
  }
  return(pb)
}

# Run NetRep software
module_quality_NetRep <- function(exprMat, adjMat, modColors, corMat = NULL, numPerms = 10000, numThreads = 20) {
  
  # Construct list of necessary data (the network and itself)
  adj_list <- list(adjMat, adjMat)
  exp_list <- list(as.matrix(exprMat), as.matrix(exprMat))
  if (is.null(corMat)) {corMat <- cor(exprMat)}
  cor_list <- list(corMat, corMat)
  
  modColors <- str_to_title(modColors)
  names(modColors) <- colnames(exprMat)
  
  # Run NetRep on given network
  prez <- modulePreservation(network=adj_list, data=exp_list, 
                             correlation=cor_list, 
                             moduleAssignments=modColors, 
                             discovery=c(1), #test=seq(2,102), 
                             nPerm=numPerms, nThreads=numThreads, 
                             selfPreservation = TRUE)

  return(prez)
}  

# Function for use with NetRep package for module permutations
#  which uses numerical module identifiers
make_num_labels <- function(modColors, exprMat){
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50));
  num_labels = match(modColors, colorOrder)-1;
  # Convert module colors to numbers
  names(num_labels) <- colnames(exprMat)
  num_labels
}




# New version that creates named vector automatically (needs original expression matrix)
get_mod_groups <- function(mod_labs, exprMat){
  names(mod_labs) <- colnames(exprMat)
  df <- data.frame(mod_labs)
  df[,2] <- names(mod_labs)
  df <- df[order(df[,1]),]
  split(df[,2], df[,1])
}


# Function to compute similarity measure between two modules q and q'
sim_measure <- function(q1, q2){
  length(intersect(q1,q2)) / min(length(q1), length(q2))
}

# Function that assigns similarity score to every reference module given bootstrap modules
bootstrap_module_similarity_scores <- function(ref_mod_nodes, boot_mod_nodes) {
  
  # For each reference module, take score of most similar bootstrap module
  ref_mod_scores <- list() # Store similarity score for each assigned (highest) bootstrap module for each reference module
  for (ref_mod in names(ref_mod_nodes)) {
    ref_nodes <- ref_mod_nodes[[ref_mod]] # Nodes of current reference module
    best_score <- 0.0
    #best_match <- NULL
    for (boot_mod in names(boot_mod_nodes)) {
      boot_nodes <- boot_mod_nodes[[boot_mod]] # Nodes of current bootstrap module
      # Find which module it is most similar to
      sim <- sim_measure(ref_nodes, boot_nodes)
      if (sim >= best_score) {
        best_score <- sim
        #best_match <- names(boot_mods[k])
      }
    }  
    ref_mod_scores[[ref_mod]] <- best_score
  }
  return(ref_mod_scores)
}  

# Function to compute h-indices from series of boostrap runs against refernce modules
get_bootstrap_hIndex_scores <- function(ref_mod_nodes, boots) {
  # Merge bootstrap reference module scores into single dataframe while also computing h-scores
  ref_mod_bootstrap_scores <- list()
  h_scores <- list()
  for (mod in names(ref_mod_nodes)) {
    mod_scores <- c()
    for (b in 1:length(boots)) {
      mod_scores <- c(mod_scores, boots[[b]][['ref_sim_scores']][[mod]])
    }
    ref_mod_bootstrap_scores[[mod]] <- mod_scores
    scores <- sort(mod_scores, decreasing = TRUE)
    h <- 0.0
    for (j in 1:length(scores)){
      index <- j/length(scores)
      if (scores[j] > index) {
        h <- index
      }
    }
    h_scores[[mod]] <- h
  }
  
  
  # Create melted dataframes of reference module scores and h-indices
  sim_scores_df <- data.frame(ref_mod_bootstrap_scores)
  hIndex_df <- data.frame(cbind(t(data.frame(h_scores)), apply(sim_scores_df, 2, median)))
  names(hIndex_df) <- c('H_index', 'MedianSim')
  hIndex_df$Module <- str_to_title(rownames(hIndex_df))
  # Order both dataframes by H-index using median as tie-breaker
  hIndex_df <- hIndex_df[order(-hIndex_df$H_index, -hIndex_df$MedianSim),]
  sim_scores_df <- sim_scores_df[,rownames(hIndex_df)]
  melt_scores <- melt(sim_scores_df)
  # Set column names and capitalize module names
  names(melt_scores) <- c('Module', 'Score')
  melt_scores$Module <- str_to_title(melt_scores$Module)
  # Return data frames
  return_list <- list()
  return_list[['bootstrap_sim_scores']] <- melt_scores
  return_list[['module_Hindex']] <- hIndex_df
  return(return_list)
}


# Function to retreive scale free values and mean connectivity values from bootstrap runs
get_scale_free_bootstrap_results <- function(boots) {
  # Save scale free fit data
  sft <- matrix(nrow = length(boots), ncol = length(powers))
  for (i in 1:length(boots)) {
    sft[i,] <- boots[[i]]$sft$SFT.R.sq
  }
  sft <- data.frame(sft)
  colnames(sft) <- c(1:dim(sft)[2])
  sft <- stack(sft)
  names(sft) <- c('SFT.R.sq', 'Power')
  sft <- sft[,c(2,1)]
  # Save mean connectivity data
  meanK <- matrix(nrow = length(boots), ncol = length(powers))
  for (i in 1:length(boots)) {
    meanK[i,] <- boots[[i]]$sft$mean.k.
  }
  meanK <- data.frame(meanK)
  colnames(meanK) <- c(1:dim(meanK)[2])
  meanK <- stack(meanK)
  names(meanK) <- c('meanK', 'Power')
  meanK <- meanK[,c(2,1)]
  return_list <- list()
  return_list[['sft']] <- sft
  return_list[['meanK']] <- meanK
  return(return_list)
}


#################################################################################
#
#     Z-scores Module Quality Functions
#
##################################################################################

# Compute z-scores for stats from NetRep permutations for a given network
get_zScores <- function(bnet){
  nperms <- dim(bnet$nulls)[3]
  zmat <- matrix(nrow = dim(bnet$nulls)[1], ncol = dim(bnet$nulls)[2], 
                 dimnames = list(rownames(bnet$nulls), colnames(bnet$nulls)))
  for (i in 1:dim(bnet$nulls)[1]){
    for (j in 1:dim(bnet$nulls)[2]){
      obs <- bnet$observed[i,j]
      null_vals <- bnet$nulls[i,j,]
      null_vals[is.na(null_vals)] <- 0.0
      z <- (obs - mean(null_vals)) / sd(null_vals)
      zmat[i,j] <- z
    }
  }
  # Compute composite and summary z-scores
  newCols <- c(colnames(bnet$nulls), 'Zdensity', 'Zconnectivity', 'Zsummary')
  zmat2 <- matrix(nrow = dim(bnet$nulls)[1], ncol = 10, 
                  dimnames = list(rownames(bnet$nulls),newCols))
  for (i in 1:dim(bnet$nulls)[1]){
    zDens <- median(c(zmat[i,6], zmat[i,1], zmat[i,2], zmat[i,7]))
    zConn <- median(c(zmat[i,4], zmat[i,5], zmat[i,3]))
    zSumm <- (zDens + zConn)/2
    zmat2[i,1:7] <- zmat[i,1:7]
    zmat2[i,8] <- zDens
    zmat2[i,9] <- zConn
    zmat2[i,10] <- zSumm
  }  
  zmat2
}

# Functiont to merge list of zScores from bootstrap runs
merge_bootstrap_zScores <- function(boots) {
  # Create a list of matrices for each zScore type
  zScore_list <- list()
  #for (z in colnames(boots[[1]])) {
  # Just do summary zscores for now
  for (z in c('Zdensity', 'Zconnectivity', 'Zsummary')) {
    zMat <- matrix(nrow = length(boots), ncol = dim(boots[[1]][['zScores']])[1], 
                   dimnames = list(seq(1:length(boots)), 
                                   rownames(boots[[1]][['zScores']])))
    # Store zScores from each bootstrap run for current zScore type
    for (i in 1:length(boots)) {
      zMat[i,] <- boots[[i]][['zScores']][,z]
    }
    zScore_list[[z]] <- data.frame(zMat)
  }
  return(zScore_list)
}

# Get max summary zScore from list of module zScores (for plotting)
get_max_zScore <- function(zScore_list) {
  maxZ <- 5
  for (i in 1:length(zScore_list)) {
    # Get summary zScores
    zMat <- zScore_list[[i]][,8:10]
    if (max(zMat) > maxZ) {
      maxZ <- round_any(max(zMat), 5, f = ceiling)
    }
  }
  return(maxZ)
}

# Get max zScore from bootstrap zScore results (for plotting)
get_max_bootstrap_zScore <- function(boot_data, zScores = NULL) {
  if (is.null(zScores)) {zScores = names(boot_data$zScores)}
  all_scores <- unlist(boot_data$zScores[zScores])
  return(max(all_scores, na.rm = TRUE))
}


###########################################################################
#
#   Stats and permutation functions for comparing module expression
#    across multiple groupes (e.g. tissues)
#
###########################################################################

# Unparallized version (only for debugging prior to setting up parallelized verison)
group_module_permutations <- function(nodeExpr, modColors, traitsDF, trait_index, nperms = 1000) {
  
  # Break up trait into binary variables
  df <- data.frame(binarizeCategoricalVariable(traitsDF[,trait_index], 
                                               includePairwise = FALSE, includeLevelVsAll = TRUE, nameForAll = ""))
  rownames(df) <- rownames(traitsDF)
  
  # Make sure nodeExpr matches traits dataframe
  nodeExpr <- nodeExpr[rownames(df),]
  
  # Number of modules to store
  nMods <- length(table(modColors[which(modColors != 'grey')]))
  
  # For each null model, store various data for each module/group  
  null_mats <- list()
  for (n in 1:nperms) {
    if ((n %% 100) == 0) print(n)
    # Create matrices that will store  results for each module and trait
    traitModScores = matrix(0, nrow=dim(df)[2], ncol=nMods)
    traitModPVals = matrix(0, nrow=dim(df)[2], ncol=nMods)
    traitModMeans = matrix(0, nrow=dim(df)[2], ncol=nMods)
    traitModMedians = matrix(0, nrow=dim(df)[2], ncol=nMods)
    # Create null module expression matrix and store
    nullColors <- modColors[sample(length(modColors), replace = FALSE)]
    nullExpr = moduleEigengenes(nodeExpr, nullColors, excludeGrey = TRUE)$eigengenes
    names(nullExpr) <- str_to_title(substring(names(nullExpr), 3))
    null_mats[[n]] <- list()
    null_mats[[n]][['Expr']] <- nullExpr
    # Parse samples
    for (t in 1:dim(df)[2]) {
      testSet <- rownames(df[df[,t] == 1,])
      restSet <- rownames(df[df[,t] == 0,])
      # Parse modules of current null expression matrix
      for (m in 1:nMods) {
        testVals <- nullExpr[testSet,m]
        restVals <- nullExpr[,m]
        testMean <- mean(testVals)
        restMean <- mean(restVals)
        testMedian <- median(testVals)
        restMedian <- median(restVals)
        # Calculate p-value
        p.val <- wilcox.test(testVals, restVals)$p.value
        traitModPVals[t,m] = p.val
        # Calculate -log10 of p-value and flip sign accordingly
        if (testMean < restMean) {
          traitModScores[t, m] = -(abs(-log10(p.val)))
        } else {
          traitModScores[t, m] = abs(-log10(p.val))
        }
        # Store group/module mean and group/module median
        traitModMeans[t, m] <- testMean
        traitModMedians[t, m] <- testMedian
      }
    }
    # Store matrices for current null model set
    traitModPVals[] <- p.adjust(traitModPVals, method = 'BH')
    null_mats[[n]][['Scores']] <- traitModScores
    null_mats[[n]][['p.val']] <- traitModPVals
    null_mats[[n]][['Means']] <- traitModMeans
    null_mats[[n]][['Medians']] <- traitModMedians
  }
  return(null_mats)
} 

# Parallelized version to run random module permutations on computing cluster using multiple CPUs
group_module_permutations_parallel <- function(nodeExpr, modColors, traitsDF = NULL, trait_index = NULL, 
                                               nperms = 1000, numCores = 32) {
  
  if (!(is.null(traitsDF))) {
    # Break up trait into binary variables
    df <- data.frame(binarizeCategoricalVariable(traitsDF[,trait_index], 
                                                 includePairwise = FALSE, includeLevelVsAll = TRUE, nameForAll = ""))
    rownames(df) <- rownames(traitsDF)
  } else {
    df <- data.frame(rep(1, dim(nodeExpr)[1]))
    rownames(df) <- rownames(nodeExpr)
  }
  
  # Make sure nodeExpr matches traits dataframe
  nodeExpr <- nodeExpr[rownames(df),]
  
  # Number of modules to store
  nMods <- length(table(modColors[which(modColors != 'grey')]))
  
  # For each null model, store various data for each module/group  
  #null_mats <- list()
  print(detectCores())
  registerDoParallel(cores=numCores)
  n <- 0
  null_mats <- foreach(icount(nperms)) %dopar% {
  #for (i in 1:nperms) {
    n <- n + 1
    if ((n %% 100) == 0) print(n)
    # Create matrices that will store  results for each module and trait
    traitModScores = matrix(0, nrow=dim(df)[2], ncol=nMods)
    traitModPVals = matrix(0, nrow=dim(df)[2], ncol=nMods)
    traitModMeans = matrix(0, nrow=dim(df)[2], ncol=nMods)
    traitModMedians = matrix(0, nrow=dim(df)[2], ncol=nMods)
    # Create null module expression matrix and store
    nullColors <- modColors[sample(length(modColors), replace = FALSE)]
    nullExpr = moduleEigengenes(nodeExpr, nullColors, excludeGrey = TRUE)$eigengenes
    names(nullExpr) <- str_to_title(substring(names(nullExpr), 3))
    nullMat <- list()
    nullMat[['Expr']] <- nullExpr
    # Parse samples
    for (t in 1:dim(df)[2]) {
      if (dim(df)[2] > 1) {
        testSet <- rownames(df[df[,t] == 1,])
        restSet <- rownames(df[df[,t] == 0,])
      } else {
        testSet <- rownames(df)
        restSet <- rownames(df)
      }
      # Parse modules of current null expression matrix
      for (m in 1:nMods) {
        testVals <- nullExpr[testSet,m]
        restVals <- nullExpr[,m]
        testMean <- mean(testVals)
        restMean <- mean(restVals)
        testMedian <- median(testVals)
        restMedian <- median(restVals)
        # Calculate p-value
        p.val <- wilcox.test(testVals, restVals)$p.value
        traitModPVals[t,m] = p.val
        # Calculate -log10 of p-value and flip sign accordingly
        if (testMean < restMean) {
          traitModScores[t, m] = -(abs(-log10(p.val)))
        } else {
          traitModScores[t, m] = abs(-log10(p.val))
        }
        # Store group/module mean and group/module median
        traitModMeans[t, m] <- testMean
        traitModMedians[t, m] <- testMedian
      }
    }
    # Store matrices for current null model set
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
    nullMat[['Scores']] <- traitModScores
    nullMat[['p.val']] <- traitModPVals
    nullMat[['Means']] <- traitModMeans
    nullMat[['Medians']] <- traitModMedians
    nullMat
  }
  return(null_mats)
} 


# Compute various module stats from observed module set and sample expression
#   --> Stats include median and mean sample expression along with two-sided test results
group_module_stats <- function(modExpr, traitsDF = NULL, trait_index = NULL) {
  
  if (!(is.null(traitsDF))) {
    # Break up trait into binary variables
    df <- data.frame(binarizeCategoricalVariable(traitsDF[,trait_index], 
                                                 includePairwise = FALSE, includeLevelVsAll = TRUE, nameForAll = ""))
    rownames(df) <- rownames(traitsDF)
  } else {
    df <- data.frame(rep(1, dim(modExpr)[1]))
    rownames(df) <- rownames(modExpr)
  }

  # Make sure trait dataframe and module expression dataframe match
  modExpr <- modExpr[rownames(df),]
    
  # Create matrix that will store test results for each module and trait
  traitModScores = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2], 
                          dimnames = list(colnames(df), colnames(modExpr)))
  traitModPVals = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2], 
                         dimnames = list(colnames(df), colnames(modExpr)))
  traitModMeans = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2], 
                         dimnames = list(colnames(df), colnames(modExpr)))
  traitModMedians = matrix(0, nrow=dim(df)[2], ncol=dim(modExpr)[2], 
                           dimnames = list(colnames(df), colnames(modExpr)))

  for (t in 1:dim(df)[2]) {
    if (dim(df)[2] > 1) {
      testSet <- rownames(df[df[,t] == 1,])
      restSet <- rownames(df[df[,t] == 0,])
    } else {
      testSet <- rownames(df)
      restSet <- rownames(df)
    }
    # Parse modules of current null expression matrix
    for (m in 1:dim(modExpr)[2]) {
      testVals <- modExpr[testSet,m]
      restVals <- modExpr[,m]
      testMean <- mean(testVals)
      restMean <- mean(restVals)
      testMedian <- median(testVals)
      restMedian <- median(restVals)
      # Calculate p-value
      p.val <- wilcox.test(testVals, restVals)$p.value
      traitModPVals[t,m] = p.val
      # Calculate -log10 of p-value and flip sign accordingly
      if (testMean < restMean) {
        traitModScores[t, m] = -(abs(-log10(p.val)))
      } else {
        traitModScores[t, m] = abs(-log10(p.val))
      }
      # Store group/module mean and group/module median
      traitModMeans[t, m] <- testMean
      traitModMedians[t, m] <- testMedian
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
  mat_list <- list()
  mat_list[['Medians']] <- traitModMedians
  mat_list[['Means']] <- traitModMeans
  mat_list[['p.val']] <- traitModPVals
  mat_list[['Scores']] <- traitModScores
  
  return(mat_list)
  
}


####################################################################
#
#     Network clustring functions
#
#####################################################################



#### Merge consensus modules and get new color labels ##
merge_consensus_modules <- function(dColors, multiExpr, dist_thres, make_plot) {
  # Calculate module eigengenes
  unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = dColors)
  # Calculate consensus dissimilarity of consensus module eigengenes
  consMEDiss = consensusMEDissimilarity(unmergedMEs);
  if (make_plot == TRUE) {
    # Cluster consensus modules
    consMETree = hclust(as.dist(consMEDiss), method = "average");
    # Plot the result
    sizeGrWindow(7,6)
    par(mfrow = c(1,1))
    plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
         xlab = "", sub = "")
    abline(h=dist_thres, col = "red")
  }
  merge = mergeCloseModules(multiExpr, dColors, cutHeight = dist_thres, verbose = 0)
  merge$colors
}

#### Perform hierarchical clustering from TOM #####
hierclust_TOM <- function(adjMat, ds, minModSize, show_tree){
  # Turn adjacency into topological overlap
  dissTOM = 1-TOMsimilarity(adjMat)
  spliceTree = hclust(as.dist(dissTOM), method = "average")
  dynamicMods = cutreeDynamic(dendro = spliceTree, distM = dissTOM, deepSplit = ds, 
                              pamRespectsDendro = FALSE, minClusterSize = minModSize);
  dColors = labels2colors(dynamicMods)
  if (show_tree == TRUE){
    # Plot the dendrogram and colors underneath
    sizeGrWindow(8,6)
    plotDendroAndColors(spliceTree, dColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Splice event dendrogram and module colors")
  }
  dColors
}

#### Get merge close modules and get module color labels ####
merge_and_get_colors <- function(dColors, exprMat, dist_thresh, make_plot){
  # Calculate eigengenes
  MEList = moduleEigengenes(exprMat, colors = dColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  if (make_plot == TRUE){
    sizeGrWindow(7, 6)
    plot(METree, main = "Clustering of module eigengenes",
         xlab = "", sub = "")
    # Plot the cut line into the dendrogram
    abline(h=dist_thresh, col = "red")
  }
  ### MERGE ####
  # Call an automatic merging function
  merge = mergeCloseModules(exprMat, dColors, cutHeight = dist_thresh, verbose = 0)
  # The merged module colors
  merge$colors
}



#### Perform spectral clustering #####
cluster_sp <- function(adjMat, k, t, l, spherical = FALSE){
  if (spherical == TRUE) {
    rez_mem <- reg.SSP(adjMat, K=k, tau=t, lap = l)$cluster
  } else {
    rez_mem <- reg.SP(adjMat, K=k, tau=t, lap = l)$cluster
  }
  
  # Get to format workable in WGCNA and such
  member_table <- data.frame(table(rez_mem))
  clustered <- member_table[member_table$Freq > 1,]
  clustered <- clustered[order(clustered$Freq, decreasing = TRUE),]
  clustered$rank <- seq(nrow(clustered))
  members <- c()
  for (i in 1:length(rez_mem)){
    mod <- rez_mem[i]
    if (member_table[member_table$rez_mem == mod,]$Freq == 1) {
      members <- c(members, 0)
    } else {
      members <- c(members, clustered[clustered$rez_mem == rez_mem[i],]$rank)
    }
  }
  memberColors = labels2colors(members)
  dColors <- memberColors
  dColors
}  

# Function to run multiple iterations of spectral clustering and assign nodes to modules
#  based on majority voting. Use on Exacloud with multiple cpus.
iterate_spectral_clustering <- function(tom, k, tau = 0, nstart = 50, iterations = 100, numcores = 25) {
  
  cluster_matrix <- function(spColors) {
    cMat <- matrix(0, nrow = length(spColors), ncol = length(spColors), dimnames = list(names(spColors), names(spColors)))
    modules <- unique(unname(spColors))
    for (mod in modules) {
      mod_nodes <- which(unname(spColors) == mod)
      cMat[mod_nodes,mod_nodes] <- 1
    }
    return(cMat)
  }
  
  # Register CPUs for DoParallel
  print(detectCores())
  registerDoParallel(cores=numCores)
  # Begin parallized iterations
  rez_list <- foreach(icount(iterations)) %dopar% {
    spColors <- network_spectral_clustering(tom, K = k, tau = tau, nstart = nstart)
    names(spColors) <- colnames(tom)
    return_list <- list()
    cMat <- cluster_matrix(spColors)
    return_list[['cMat']] <- cMat
    return_list[['spColors']] <- spColors
    return_list
  }  
  
  cluster_mats <- list()
  minSize <- dim(tom)[1]
  for (i in 1:length(rez_list)) {
    cluster_mats[[i]] <- rez_list[[i]][['cMat']]
    spCols <- min(table(rez_list[[i]][['spColors']]))
    if (spCols < minSize) {minSize <- spCols}
  }
  
  # Now sum results of each clustering assignments and select most frequent
  mat_sum <- Reduce('+', cluster_mats)
  diag(mat_sum) <- 0.0
  mat_final <- matrix(0, nrow = dim(mat_sum)[1], ncol = dim(mat_sum)[2], dimnames = list(rownames(mat_sum), colnames(mat_sum)))
  # Highest neighbors must occur so many times
  mat_final[which(apply(mat_sum, 2, function(x) x == as.numeric(max(names(table(x)[table(x) >= minSize*.8])))))] <- 1.0
  
  # Find clusters from most frequent
  rez_mem <- unname(components(graph_from_adjacency_matrix(mat_final))$membership)

  # Format node to module assignments and label using arbitrary colors
  member_table <- data.frame(table(rez_mem))
  clustered <- member_table[member_table$Freq > 1,]
  clustered <- clustered[order(clustered$Freq, decreasing = TRUE),]
  clustered$rank <- seq(nrow(clustered))
  members <- c()
  for (i in 1:length(rez_mem)){
    mod <- rez_mem[i]
    if (member_table[member_table$rez_mem == mod,]$Freq == 1) {
      members <- c(members, 0)
    } else {
      members <- c(members, clustered[clustered$rez_mem == rez_mem[i],]$rank)
    }
  }
  return(labels2colors(members))
}

# Function to perform (normalized) spectral clustering for network given topological
# overlap matrix and value of k. This is a non-iterative version (not preferred).
network_spectral_clustering <- function (A, K, tau = 0, lap = TRUE, nstart = 50) {
  avg.d <- mean(colSums(A))
  A.tau <- A + tau * avg.d/nrow(A)
  if (!lap) {
    SVD <- irlba::irlba(A.tau, nu = K, nv = K)
    V <- SVD$v[, 1:K]
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm) %*% V
  }
  else {
    d.tau <- colSums(A.tau)
    L.tau <- diag(1/sqrt(d.tau)) %*% A.tau %*% diag(1/sqrt(d.tau))
    SVD <- irlba::irlba(L.tau, nu = K, nv = K)
    V <- SVD$v[, 1:K]
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm) %*% V
  }
  
  # Run k-means on normalized eigenvectors
  set.seed(101)
  km <- kmeans(V.normalized, centers = K, nstart = nstart, iter.max = 50)
  rez_mem <- km$cluster
  
  # Get to format workable in WGCNA and such
  member_table <- data.frame(table(rez_mem))
  clustered <- member_table[member_table$Freq > 1,]
  clustered <- clustered[order(clustered$Freq, decreasing = TRUE),]
  clustered$rank <- seq(nrow(clustered))
  members <- c()
  for (i in 1:length(rez_mem)){
    mod <- rez_mem[i]
    if (member_table[member_table$rez_mem == mod,]$Freq == 1) {
      members <- c(members, 0)
    } else {
      members <- c(members, clustered[clustered$rez_mem == rez_mem[i],]$rank)
    }
  }
  return(labels2colors(members))
}


######################################################
#
#    Helper functions for setting color schemes
#
######################################################

# Function to determine if text should be white or black
#  depending on shade of background color
isDark <- function(colr) { 
  if (sum( col2rgb(colr) * c(299, 587,114))/1000 < 123) {
    return("White")
  } else {
    return("Black")
  }
}

# Functions to set various color paletes for sample annotations
# Use piratepal package and specify color palette
set_colors <- function(x, col_name){
  my_colors <- setNames(piratepal(palette = col_name, 
                                  length.out = length(unique(x))),
                        unique(x))
  return(my_colors)
}

# Select random colors using brewer.pal
rand_colors <- function(x){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  n = length(unique(x))
  my_colors <- setNames(sample(col_vector, n), unique(x))
  return(my_colors)
}

# Select colors from rainbow color palette
rb_colors <- function(x){
  my_colors <- setNames(rainbow(length(unique(x))), unique(x))
  return(my_colors)
}


# Create named vector for color coding modules in plots
module_colors <- labels2colors((c(0:200)))
names(module_colors) <- str_to_title(module_colors)

