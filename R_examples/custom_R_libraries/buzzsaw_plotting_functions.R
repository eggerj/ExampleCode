#########################################################################################
#
#       PLOTTING FUNCTIONS SPECIFICALY FOR CREATING BUZZSAW PLOTS
#
#         Also contains other helper and plotting functions for processing 
#           results from LASSO drug prediction models
#
#########################################################################################


set.seed(101)



# Function to generate random coefficient ratio data given sample size
#   Returns dataframe and other info for plotting
#  Used for debugging buzzsaw plotting functions
generate_data <- function(n = 6, group_names = NULL) {
  
  if (is.null(group_names)) {
    groups <- letters[1:n]
  } else {
    groups <- group_names[1:n]
  }
  
  combs <- t(combn(groups, 2, simplify = TRUE))
  
  newMat <- matrix(0, nrow = n*(n-1), ncol = 4)
  
  coPos = runif(dim(combs)[1], min = 0, max = 0.5)
  coNeg = runif(dim(combs)[1], min = 0, max = 0.5)
  
  for (i in 1:dim(combs)[1]) {
    newMat[i,] <- c(combs[i,], coPos[i], coNeg[i])  
    j <- i + dim(combs)[1]
    x <- runif(1, min = 0.01, max = 0.99)
    p <- coPos[i] + coNeg[i]
    q <- p*x
    p <- p-q
    
    newMat[j,] <- c(combs[i,c(2,1)], p,q)
  }
  df <- data.frame(group1 = newMat[,1], group2 = newMat[,2],
                   coPos = as.numeric(newMat[,3]), coNeg = as.numeric(newMat[,4]),
                   coLink = (as.numeric(newMat[,3]) + as.numeric(newMat[,4])))
  df <- df[order(-df$coLink),]
  
  # Get ranks for link positions
  df <- df %>% group_by(group1) %>%
    mutate(coRank = (order(order(coLink, decreasing=TRUE)))-1) %>% 
    data.frame()
  
  # Use first for real data, use second more exaggerating random data
  #totalRatio <- data.frame(aggregate(coPos+coNeg~group1, df, sum))
  totalRatio <- data.frame(factors = groups, coefRatio = runif(n, min = 25.00001, max = 100.00))#sample(25:100, n))

  # Factors must be set so highest (num non-zero coefficients) group is first
  totalRatio <- totalRatio[order(-totalRatio[,2]),]
  df$group1 <- factor(df$group1, levels = totalRatio$factors[order(-totalRatio[,2])])
  
  # Add totalRatio to corresponding rows (duplicates for this column)
  df$coefRatio <- totalRatio$coefRatio[match(df$group1, totalRatio$factors)]
  
  sector.width <- totalRatio[,2] / sum(totalRatio[,2])
  
  return_list <- list()
  return_list[['df']] <- df
  return_list[['sector.width']] <- sector.width # We could maybe just do this part again in barplot function
  return(return_list)
}



# Create dataframe describing co-occurrences from coefficient matrix of given drug
make_co_occurence_dataframe <- function(coef_mat) {
  
  # Function to make sure modules are uppercase to match color assignments
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  rownames(coef_mat) <- firstup(rownames(coef_mat))
  
  # Number of bootstrap iterations
  boots <- dim(coef_mat)[2]
  
  # Select features that have at least 25% of non-zeros
  coef_ratios <- rowSums(coef_mat != 0) / boots    # keep for later so we can get total ratios
  to_keep <- names(coef_ratios[coef_ratios >= 0.25])
  coef_mat <- coef_mat[to_keep,]
  groups <- rownames(coef_mat)[rownames(coef_mat) != '(Intercept)']
  
  # Create new matrix for storing
  n <- length(groups)
  newMat <- matrix(0, nrow = n*(n-1), ncol = 6)
  
  # Generate all pairwise combinations (doesn't include reverse duplicate yet)
  combs <- t(combn(groups, 2, simplify = TRUE))
  # Now compute co-occuring frequences and store in new matrix
  for (i in 1:dim(combs)[1]) {
    
    # Select current pair
    g1 <- combs[i,1]
    g2 <- combs[i,2]
    sharedCo <- coef_mat[c(g1,g2),]
    sharedCo <- sharedCo[,colSums(sharedCo != 0.0) == 2]
    
    # Have to deal with R changing single columns to vectors
    if (!is.null(dim(sharedCo))) {
      coPos <- sum(sharedCo[g1,] > 0.0) / boots
      coNeg <- sum(sharedCo[g1,] < 0.0) / boots
    } else {
      coPos <- sum(sharedCo[g1] > 0.0) / boots
      coNeg <- sum(sharedCo[g1] < 0.0) / boots
    }
    
    # Add to new matrix
    newMat[i,] <- c(combs[i,], coPos, coNeg, (coPos + coNeg), coef_ratios[g1])
    
    # Do the same thing but flip the groups
    j <- i + dim(combs)[1]
    # Have to deal with R changing single columns to vectors
    if (!is.null(dim(sharedCo))) {
      coPos <- sum(sharedCo[g2,] > 0.0) / boots
      coNeg <- sum(sharedCo[g2,] < 0.0) / boots
    } else {
      coPos <- sum(sharedCo[g2] > 0.0) / boots
      coNeg <- sum(sharedCo[g2] < 0.0) / boots
    }
    newMat[j,] <- c(combs[i,c(2,1)], coPos, coNeg, (coPos + coNeg), coef_ratios[g2])
  }  
  
  # Store data in dataframe and compute ranks for plotting
  df <- data.frame(group1 = newMat[,1], group2 = newMat[,2],
                   coPos = as.numeric(newMat[,3]), coNeg = as.numeric(newMat[,4]),
                   coLink = as.numeric(newMat[,5]), coefRatio = as.numeric(newMat[,6]))
  df <- df[order(-df$coLink),]
  
  # Get ranks for link positions
  df <- df %>% group_by(group1) %>%
    mutate(coRank = (order(order(coLink, decreasing=TRUE)))-1) %>% data.frame()
  
  return(df)
  
}


# Function to initialize circos object and create barplots with group labels
make_circoBP <- function(df) {
  
  # First sort groups highest to lowest based on coefRatio
  df$group1 <- factor(df$group1, levels = unique(df[order(-df$coefRatio),]$group1))

  # Now create sector widths for plot
  sector.width <- df[!duplicated(df[,c('group1','coefRatio')]),]
  sector.width <- sector.width[order(sector.width$group1),]$coefRatio
  sector.width <- sector.width / sum(sector.width)
  
  nObs <- max(table(df$group1))
  
  # gap.degree currently 0.08
  circos.par(start.degree = 90, gap.degree = 0.00000, track.height = 0.15, # Default gap.degree = 1, track.height = 0.2
             track.margin = c(0.0025, 0.0025), cell.padding	= c(0.01, 1.0, 0.01, 1.0)) # Default track.margin = c(0.1, 0.1), cell.padding = c(0.02, 1.00, 0.02, 1.00)
  circos.initialize(df$group1, xlim = c(0,nObs), sector.width = sector.width)
  
  circos.track(df$group1, x = df$coPos, y = df$coNeg, ylim = c(0,1), # Sum of two boxplot categories should always be 1 or less
               panel.fun = function(x, y) {
                 
                 # Create matrix for stacked barplot
                 value = as.matrix(cbind(x,y))
                 # Make sure its sorted highest to lowest
                 value = value[order(rowSums(value),decreasing = TRUE),]
                 bar_cols <- c(2,4)
                 bar_cols <- c("#EE0011FF","#0C5BB0FF")
                 circos.barplot(value, 0:(nObs-1)+0.5, col = bar_cols, bar_width = 1)#col = ifelse(x > 0.25, 2, 4))
                 
                 # Make labels
                 name = get.cell.meta.data("sector.index")
                 xlim = get.cell.meta.data("xlim")
                 theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                 dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                 aa = c(1, 0.5)
                 if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                 
                 #plot sector labels
                 circos.text(x=mean(xlim), y=1.1, labels=name, facing = dd, cex=0.9,  adj = aa)

               }, bg.border = NA)
}

# Function to create a label indicating highest co-occurence partner for each feature
label_largest_co_occurence <- function(df) {
  # Now try adding lines, segments and text for annotation of values
  circos.track(df$group1, x = df$coLink, ylim = c(0, 0.8), track.height = 0.045, 
               panel.fun = function(x, y) {
                 
                 nObs <- max(table(df$group1))
                 
                 # Add label showing ratio of highest co-occurring co-efficient
                 maxCo <- max(x)
                 circos.text(x=0, y=0.4, labels=sprintf("%1.0f%%", 100*maxCo), facing = "inside", niceFacing = TRUE, cex=0.8,  adj = c(0,0.5))
                 
               }, bg.border = NA)
}


# Function to colored segments relative to frequency for each feature
add_colored_segments <- function(df, col_mat = NULL) {
  
  if (is.null(col_mat)) {
    col_mat = rand_color(length(unique(df$group1)))
    names(col_mat) <- unique(df$group1)
  }
  
  # Try adding colored segments
  circos.track(df$group1, x = df$coefRatio, ylim = c(0, 1), track.height = 0.08, 
               panel.fun = function(x, y) {
                 
                 nObs <- max(table(df$group1))
                 
                 # Get colors
                 i = get.cell.meta.data("sector.numeric.index")
                 name <- get.cell.meta.data("sector.index")
                 rec_col = col_mat[name]
                 txt_col = isDark(rec_col)
                 
                 #plot main sector
                 circos.rect(xleft=0, ybottom=0, xright=nObs, ytop=1,
                             col = rec_col, border=rec_col)
                 
                 coefRatio <- unique(x)
                 circos.text(x=nObs/2, y=0.5, labels=sprintf("%1.0f%%",coefRatio*100), col = txt_col, facing = "inside", niceFacing = TRUE, cex=0.8,  adj = c(0.5,0.5))
                 
                 
               }, bg.border = NA)
}

# Add links indicating co-occurrences between features
add_links_seq <- function(df, col_mat = NULL) {
  
  if (is.null(col_mat)) {
    col_mat = rand_color(length(unique(df$group1)))
    names(col_mat) <- unique(df$group1)
  }
  
  #  Create color matrix: first sort df by first two columns and assign color spectrum for each group
  #  highest co-occurence gets color as is, lowest gets lighten to highest level (scale between 0 and 1 for all values)
  link_order <- df[order(df$coefRatio, df$coLink),]

  assigned <- matrix(0, nrow = length(unique(link_order$group1)), ncol = length(unique(link_order$group1)),
                     dimnames = list(unique(link_order$group1), unique(link_order$group1)))
  diag(assigned) <- NA

  # Fill matrix from highest lowest co-occurence pairs; skip pair if color already assigned  
  for (g1 in unique(link_order$group1)) {
    cos <- link_order[link_order$group1 == g1,] # Current sector
    n <- sum(!is.na(assigned[g1,])) # Count of links to assign colors to
    inc <- 1 / dim(cos)[1] # Light intensity increment value of links for group
    if (n != 0) {
      lte <- (n-1)*inc
      for (g2 in cos$group2) {
        if (!is.na(assigned[g1,g2])) {
          assigned[g1,g2] <- adjustcolor(lighten(col_mat[g1], amount = lte), alpha.f = 0.9)
          assigned[g2,g1] <- NA   # Make sure we don't duplicate
          lte <- lte - inc
        }
      }
    }
  }
  
  # Parse all possible group pairings and add link using color assignment matrix; skip pair if link already exists (reverse pair)
  combs <- t(combn(unique(df$group1), 2, simplify = TRUE))
  for (i in 1:dim(combs)[1]) {
    g1 <- combs[i,1]
    g2 <- combs[i,2]
    
    if (!is.na(assigned[as.character(g1),as.character(g2)])) {
      
      g1I <- c(df[df$group1 == g1 & df$group2 == g2,]$coRank, df[df$group1 == g1 & df$group2 == g2,]$coRank + 1)
      g2I <- c(df[df$group2 == g1 & df$group1 == g2,]$coRank, df[df$group2 == g1 & df$group1 == g2,]$coRank + 1)
      circos.link(g1, g1I, g2, g2I, col = assigned[as.character(g1),as.character(g2)])
      
    } else {  
      
      g1 <- combs[i,2]
      g2 <- combs[i,1]
      
      if (!is.na(assigned[as.character(g1),as.character(g2)])) {
        
        g1I <- c(df[df$group1 == g1 & df$group2 == g2,]$coRank, df[df$group1 == g1 & df$group2 == g2,]$coRank + 1)
        g2I <- c(df[df$group2 == g1 & df$group1 == g2,]$coRank, df[df$group2 == g1 & df$group1 == g2,]$coRank + 1)
        circos.link(g1, g1I, g2, g2I, col = assigned[as.character(g1),as.character(g2)])
        
      }
    }
  }
}

##########################################################
#
#   Make buzzsaw color schemes
#
##########################################################

# This is from manually creating useful color annotations for mutations to go with module colors
#save(original_buzz_colors, file = '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/rData_files/beatAML_rData/lasso_data/original_buzz_colors.RData')
load(file = '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/rData_files/beatAML_rData/lasso_data/original_buzz_colors.RData')

# Function to assign colors for buzzsaw plot
buzzsaw_colors <- function(all_features, start = 35) {
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  all_features <- firstup(sort(as.character(all_features)))
  
  # Create named vector for color coding modules in plots
  module_colors <- labels2colors((c(0:300)))
  
  brewer.pal(n = 8, name = "Dark2")
  module_colors <- c(brewer.pal(n = 8, name = "Dark2"), unname(piratepal(palette = "basel")), unname(piratepal(palette = "eternal")), 
                     module_colors[start:300], module_colors[0:(start-1)])
  
  # Manually remove bad colors (done durring debugging stage)
  lose <- c(6,7,12,16,18,24,25,27,30,31,32,33,34,35,36,37,38,39,41,42,43,44,46, which(module_colors == 'honeydew1'))
  module_colors <- module_colors[-lose]
  names(module_colors) <- str_to_title(module_colors)
  
  # Assign colors: if it is a module use module color label
  i <- 1
  col_assign <- c()
  for (f in all_features) {
    if (f %in% names(module_colors)) {
      if (f == 'Lightcyan' || f == 'Lightyellow') {
        col_assign[f] = darken(module_colors[f], amount = 0.1)
      } else {
        col_assign[f] = module_colors[f]
      }
    } else {
      col_assign[f] = module_colors[i]
      i <- i + 1
    }
  }
  
  module_colors <- module_colors[!module_colors %in% original_buzz_colors]
  
  # Reassign
  i <- 1
  for (f in all_features) {
    if (f %in% names(original_buzz_colors)) {
      col_assign[f] = original_buzz_colors[f]
    } else {
      col_assign[f] = module_colors[i]
      i <- i + 1
    }
  }
  return(col_assign)
  
}


####################################################################
#
#       Other functions for processing LASSO results
#
####################################################################

# Look for missing files in directory
get_missing <- function(inPath, drugDict) {
  # First get missing files
  missing_files <- c()
  for (drug in names(drugDict)) {
    fn <- paste(inPath, drug, '.RData', sep = '')
    if (!file.exists(fn)) {
      #print(unname(drugDict[drug]))
      missing_files <- c(missing_files, unname(drugDict[drug]))
    }
  }
  print(length(missing_files))
  return(missing_files)
}


# Compute R-squared 
rSquared <- function(test, pred) {
  rss <- sum((test - pred) ^ 2)  ## residual sum of squares
  tss <- sum((test - mean(test)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}



# Load and merge data for each drug model
load_data <- function(inPath, drugDict) {
  
  coef.list.all <- list()
  drugRsquaredVals <- data.table(Inhibitor = character(), Iteration = numeric(), Rsquared = numeric())
  
  for (drug in names(drugDict)) {
    fn <- paste(inPath, drug, '.RData', sep = '')
    if (file.exists(fn)) {
      load(fn)
      if (!is.null(pred.sum)) {
        coef.list.all[[drugDict[drug]]] <- coef.list[[1]]
        iterRs <- sapply(split(pred.sum, pred.sum$index), function(X) rSquared(X$test_auc, X$pred_auc))
        drugVals <- data.frame(Inhibitor = rep(drug, 100), Iteration = names(iterRs), Rsquared = iterRs)
        drugRsquaredVals <- rbind(drugRsquaredVals, drugVals)
      }
    }
  }
  drugRsquaredVals <- data.frame(drugRsquaredVals)
  
  return_list <- list()
  return_list[['coefs']] <- coef.list.all
  return_list[['r2']] <- drugRsquaredVals
  return(return_list)
}


# Plot R-squared values
plot_rSquared <- function(drugRsquaredVals) {
  
  yMax <- round_any(max(drugRsquaredVals$Rsquared), .1, f = ceiling)
  yMin <- round_any(min(drugRsquaredVals$Rsquared), .1, f = floor)
  
  # Plot distribution of R^2 values from 100 model iterations for each drug
  p <- ggplot(drugRsquaredVals) + geom_hline(yintercept = 0, color = "black") + 
    geom_boxplot(aes(x = reorder(Inhibitor,-Rsquared,FUN = median), y = Rsquared)) + 
    ggtitle(expression(paste(R^2, ' Values of LASSO Drug Prediction Models (100 Iterations Per Drug Model)'))) +
    scale_y_continuous(limits = c(yMin,yMax), breaks = seq(yMin,yMax,0.1)) + labs(x = 'Inhibitor', y = expression(R^2)) +
    theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  return(p)
  
}


# Get all features 
get_features <- function(coef.list.all) {
  all_features <- c()
  for (drug in names(coef.list.all)) {
    all_features <- unique(c(all_features, rownames(coef.list.all[[drug]])))
  }
  # Remove intercept
  all_features <- all_features[all_features != '(Intercept)']
  
  return(all_features)
}


# Function to create a matrix of non-zero coefficient frequencies for each LASSO model
make_coefficient_matrix <- function(coef.list.all, all_features, thres = 0) {
  
  # Get non-zero rate for each feature
  drugCoefficients <- matrix(0, nrow = length(names(coef.list.all)), ncol = length(all_features),
                             dimnames = list(names(coef.list.all), all_features))
  for (drug in names(coef.list.all)) {
    coeff_ratios <- rowSums(coef.list.all[[drug]]!=0) / dim(coef.list.all[[drug]])[2]
    for (cr in names(coeff_ratios)) {
      if (cr != '(Intercept)') {
        drugCoefficients[drug, cr] <- coeff_ratios[cr]
      }
    }
    for (f in all_features) {
      if (!(f %in% names(coeff_ratios))) {drugCoefficients[drug, f] <- NA}
    }
  }
  drugCoefficients <- data.frame(drugCoefficients)
  
  # Filter drugs not having a coefficient ratio at least threshold defined
  drugCoefficients <- drugCoefficients[rowSums(drugCoefficients >= thres, na.rm = TRUE) >=1,]
  print(dim(drugCoefficients))
  
  return(drugCoefficients)
  
}


# Function to plot heatmap of LASSO model results
coefficient_ratio_heatmap <- function(mat, select_rows = NULL, cluster_columns = FALSE, title = '', fontSize = 9) {

  if (!(is.null(select_rows))) {mat <- mat[select_rows,]}
  
  ht <- Heatmap(mat, cluster_rows = FALSE, cluster_columns = cluster_columns, col = colorRamp2(c(0,1), c("white", "red")),
                name = 'Frequency of Non-zero Coefficients (100,000 Models)', column_title = title, row_title = 'Inhibitor LASSO Models',
                column_names_rot = 45,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = fontSize))},
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(10, "cm"),
                                            legend_height = unit(10, "cm"), title_position = "topcenter"))
  
  return(ht)

}


# Function to plot sample sizes used for each drug
plot_sample_sizes <- function(split.drugs, clin.dt, modExpr.aml.originalID, wes.clin, select_drugs = NULL) {
  
  if (!(is.null(select_drugs))) {split.drugs <- split.drugs[select_drugs]}
  
  sampleSizes <- sapply(split.drugs, function(x){
    cur.drug <- copy(x)
    use.samps <- get_earliest_samples(clin.dt, rownames(modExpr.aml.originalID) , wes.clin$SampleID, cur.drug$lab_id)
    use.samps[,.N]
  })
  
  print(length(sampleSizes))
  sampleSizes <- sampleSizes[sampleSizes >= 150]
  print(length(sampleSizes))
  print(max(sampleSizes))
  sampleSizes <- data.frame(Inhibitor = names(sampleSizes), Size = unname(sampleSizes))
  
  p <- ggplot(sampleSizes) + geom_bar(aes(x = Inhibitor, y = Size), stat = 'identity') + 
    geom_hline(yintercept = 150, color = "black") + ggtitle('Inhibitor Model Sample Sizes') +
    scale_y_continuous(limits = c(0,275), breaks = seq(0,275,25)) + labs(x = 'Inhibitor', y = 'Size') +
    theme_linedraw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(p)
}


# Function to group LASSO coefficient frequencies by drug family and plot ratios
plot_drug_family_frequencies <- function(drugCoefficients, drug_fams) {

    # Store Drug Family LASSO Frequencies
  DrugFamilyLassoCoefs <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("drugFamily", "Coef", "Freq"))))
  # Parse drug familes
  for (f in unique(drug_fams$family)) {
    
    # Select drugs in current family
    famDF <- drug_fams[drug_fams$family == f,]
    
    # Select drugs in current family that are also in LASSO models
    famDF <- famDF[famDF$inhibitor %in% rownames(drugCoefficients),]
    
    # Only select drug families with more than one drug
    if (dim(famDF)[1] > 1) {
      
      famDrugCoefs <- drugCoefficients[rownames(drugCoefficients) %in% famDF$inhibitor,]
      famDrugCoefs[is.na(famDrugCoefs)] <- 0
      
      # Calculate ratio in which each coefficient freq is greater than 0.5 for current family
      famCoefRatio <- colSums(famDrugCoefs > 0.5) / dim(famDrugCoefs)[1]
      newDF <- data.frame(drugFamily = rep(f, dim(drugCoefficients)[2]), Coef = colnames(drugCoefficients), Freq = famCoefRatio)
      
      if (max(newDF$Freq > 0)) {
        DrugFamilyLassoCoefs <- rbind(DrugFamilyLassoCoefs, newDF)
      }
    }
  }
  
  selectDrugFamilyCoefs <- DrugFamilyLassoCoefs[DrugFamilyLassoCoefs$Freq > 0,] 
  
  
  # Order columns highest to lowest non-zeros
  selectDrugFamilyCoefs$Coef <- factor(selectDrugFamilyCoefs$Coef, levels = names(sort(table(selectDrugFamilyCoefs$Coef), decreasing = TRUE)))
  # Order rows from lowest to highest non-zeros
  selectDrugFamilyCoefs$drugFamily <- factor(selectDrugFamilyCoefs$drugFamily, levels = names(sort(table(selectDrugFamilyCoefs$drugFamily),decreasing = TRUE)))
  
  # Make ggplot object and return
  p <- ggplot(selectDrugFamilyCoefs, aes(x = Coef, y=drugFamily)) +
    geom_point(aes(size=Freq)) +
    labs(y = "Drug Family", x = "LASSO Coefficient", size = 'Frequency') +
    theme_light() + ggtitle('Frequency of coefficients seen in > 50% of Lasso runs in drug families') +
    theme(plot.title = element_text(size = 16), legend.key.size = unit(2, "line"),
          axis.text=element_text(size=12), axis.title=element_text(size=14),
          axis.text.x = element_text(angle = 45, hjust=1))
  p
}


# Function to create color mapping indicating what modules any genes with mutations are found
mutation_to_module <- function(node_stats, all_features) {
  
  modules <- unique(node_stats$Module)
  
  mapMat <- matrix(0, nrow = length(modules), ncol = length(all_features[!all_features %in% modules]),
                   dimnames = list(modules, all_features[!all_features %in% modules]))
  
  for (g in colnames(mapMat)) {
    for (mod in rownames(mapMat)) {
      if (g %in% node_stats[node_stats$Module == mod,]$GeneName) {mapMat[mod,g] <- mod}
    }
  }
  
  mDF <- melt(mapMat)
  names(mDF) <- c('SplicingModule', 'Mutation', 'Module')
  p <- ggplot(mDF, aes(SplicingModule, Mutation)) + geom_tile(aes(fill = Module), colour = "black") + 
    scale_fill_manual(values=c(module_colors, c('0' = 'white')))+
    theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(p)
}
