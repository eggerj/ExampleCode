############################################################################################
#
#   LASSO Modeling of Drug Response
#
#    This script is based on methodology used in Tyner et al. 2017 
#     for modeling network modules and variant calls against drug response
#
#    Main function (run_combined_datatypes) takes all sample data and initiates 
#     remainder of functions
#
############################################################################################


source ("/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/project_Rscripts/shared_Rscripts/project_helper_functions.R")
library(glmnet)
set.seed(101)


# Function load and pre-process all sample data (can use ahead of time and save for all drug types)
preprocess <- function(){
  
  # Sample clinical info
  clin.dt <- fread("All_sample_data.csv", sep=",", header=T, skip=3)
  
  # Process drugs (from paper supplement/sage)
  final.drugs <- data.table(openxlsx::read.xlsx("TableS10_drug_response.xlsx"))
  split.drugs <- split(final.drugs, by="inhibitor")
  
  # Process variants from samples
  var.dt <- fread("All_variants.csv", sep=",", header=T, skip=2)
  
  # Collect only samples having whole-exome seq data
  wes.clin <- clin.dt[Included_2018_DNAseqAnalysis == "y"]
  
  # Add in the fusion data as well
  wes.clin[,new_value:=make.names(WHO_Fusion)]
  
  # Make fusion categories as symbol
  cat.dt <- var.dt[NA][seq_len(wes.clin[,.N])]
  cat.dt$sample_id <- wes.clin$SampleID
  cat.dt$gene <- wes.clin$new_value
  cat.dt <- cat.dt[gene %in% c("None", "Unknown")==F]
  
  # Merge
  var.dt <- rbind(var.dt, cat.dt)
  
  # Get gene co-expression modules from paper for comparison later
  use.mes <- as.matrix(read.delim("wgcna_modules.txt", sep="\t", header=T))
  rna.clin <- clin.dt[Included_2018_RNAseqAnalysis == "y"]
  
  
  # Save and load later
  save(split.drugs, rna.clin, wes.clin, use.mes, var.dt, clin.dt, 
       file="/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/beatAML_files/BeatAML_github_data/preproc_data.RData")
  
}


# Function to select training and test sets and then create bootstrap sets
get_boot_list <- function(use.mat, auc.vec, drug, stand, numIterations = 100, numBoots = 1000){
  
  # Select train and test (75/25)
  lapply(1:numIterations, function(y){
    train.samp <- sample.int(length(auc.vec), size=round(length(auc.vec)*.75), replace=F)
    test.samp <- setdiff(seq_along(auc.vec), train.samp)
    
    # Generate bootstrap samples
    boot.vals <- lapply(1:numBoots, function(y){
      
      boot.samps <- sample.int(n=length(train.samp), size=length(train.samp), replace=T)
      
      # Divide the bootstrapped samples into 10-folds 
      #   --> https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet_beta.pdf
      foldid <- sample(1:10,size=length(boot.samps),replace=TRUE)
      
      list(iter=y, boot=boot.samps,foldid=foldid)
    })
    
    list(rep=y, mat=use.mat, auc=auc.vec, train=train.samp, test=test.samp, boots=boot.vals, drug=drug, stand=stand)
  })
}

# Function to remove any duplicate patient samples
get_earliest_samples <- function(pat.info, ...){
  
  samp.lists <- list(...)
  
  use.samps <- copy(pat.info)
  
  for(i in seq_along(samp.lists)){
    
    # Allow for dataframes to be null for select type
    if (!(is.null(samp.lists[[i]]))) {
      use.samps <- use.samps[SampleID %in% samp.lists[[i]]]
    }
  }
  
  # Select data
  if (nrow(use.samps) > 0){
    
    use.samps[,collection_diff:=suppressWarnings(as.numeric(TimeOfCollectionRelativeToInclusion))]
    use.samps <- use.samps[is.na(collection_diff)==F]
    
    min.samps <- use.samps[,.(collection_diff=min(collection_diff), all_dates=paste(collection_diff, collapse=",")),by=PatientID]
    samps.m <- merge(use.samps, min.samps[,.(PatientID, collection_diff)], by=c("PatientID", "collection_diff"))
    
    # To get rid of duplicate patients (which are usual PB and BM), sort by 
    # tissue type such that BM is first and then !duplicated
    fin.samps <- samps.m[order(SpecimenType)][!duplicated(PatientID)]
    stopifnot(nrow(min.samps) == nrow(fin.samps))
    
  }else{
    
    fin.samps <- copy(use.samps)  
    
  }
  fin.samps
}

# Function to run LASSO for current set
run_lasso <- function(auc.vec, tmp.var.mat, train.samp, test.samp, boot.vals, stand, drug){
  
  stopifnot(length(auc.vec) == nrow(tmp.var.mat))
  
  train.mat <- tmp.var.mat[train.samp,]
  train.auc <- auc.vec[rownames(tmp.var.mat)][train.samp]
  
  #train models on the bootrap samples using 10-fold cv (The 10 folds provided)
  boot.res <- lapply(boot.vals, function(y){
    
    # Use try incase regression does not converge
    test.cv <- try(cv.glmnet(x=train.mat[y$boot,],y=as.numeric(train.auc[y$boot]),alpha=1, foldid=y$foldid, family="gaussian", standardize=stand, intercept=T))
    while (class(test.cv) == 'try-error') {
      newBoot <- sample.int(n=length(y$boot), size=length(y$boot), replace=T)
      newfoldid <- sample(1:10,size=length(newBoot),replace=TRUE)
      test.cv <- try(cv.glmnet(x=train.mat[newBoot,],y=as.numeric(train.auc[newBoot]),alpha=1, foldid=newfoldid, family="gaussian", standardize=stand, intercept=T))
    }
    
    list(boot=y$iter, model=test.cv)
    
  })
  
  # Predict using the test data
  test.res <- do.call(rbind, lapply(boot.res, function(y){
    
    preds.test <- predict(y$model, newx=tmp.var.mat[test.samp,], s="lambda.1se")
    data.table(boot=y$boot, test_samp=test.samp, pred_auc=preds.test[,1])
    
  }))
  
  # Average the predictions
  test.sum <- test.res[,.(pred_auc=mean(pred_auc)),by=test_samp]
  
  tmp.dt <- merge(data.table(test_samp=test.samp, test_auc=auc.vec[rownames(tmp.var.mat)][test.samp]), test.sum, by="test_samp")
  tmp.dt$inhibitor <- drug
  
  # Capture the resulting coefficients
  coef.mat <- do.call(cbind, lapply(boot.res, function(y){
    as.matrix(coef(y$model, s="lambda.1se"))
  } ))
  
  list(preds=tmp.dt, coefs=coef.mat, models=boot.res, drug=drug)
}


# Main function to setup training and test sets and perform LASSO
run_combined_datatypes <- function(split.drugs, rna.clin, wes.clin, use.mes, var.dt, clin.dt, outFile, 
                                   numIterations = 100, numBoots = 1000, numCores = 20, minMut = .05){
  
  # Merge data for lasso 
  run.comb <- sapply(split.drugs, function(x){
    
    # Get current drug
    cur.drug <- copy(x)
    message(cur.drug$inhibitor[1])
    
    # Select samples (remove duplicates)
    use.samps <- get_earliest_samples(clin.dt, rownames(use.mes) , wes.clin$SampleID, cur.drug$lab_id)
    
    # For log
    print(use.samps[,.N])
    
    # Run in at least 150 samples have all necessary data
    if (use.samps[,.N] >= 150){
      
      fin.drugs <- cur.drug[lab_id %in% use.samps$SampleID]
      
      # Allow for mutations dataframe to be NULL so only modules can be used for comparison of results
      if (!(is.null(var.dt))) {
        
        use.vars <- var.dt[sample_id %in% use.samps$SampleID]
        use.vars <- use.vars[,samp_count:=length(unique(sample_id)), by=gene][samp_count >= ceiling(minMut*use.samps[,.N])]
        use.vars[,samp_fac:=factor(sample_id, levels=use.samps$SampleID)]
        use.vars[,value:=1]
        tmp.var.mat <- reshape2::acast(samp_fac~gene,data=use.vars, value.var="value", fun.aggregate = function(x) as.integer(length(x) > 0), drop=F)
        
        # Allow for modules dataframe to be NULL so only mutations (var.dt) can be used
        if (!(is.null(use.mes))) {
          
          use.exprs <- use.mes[use.samps$SampleID,]
          stopifnot(all(rownames(tmp.var.mat) %in% rownames(tmp.var.mat)) && nrow(tmp.var.mat) == nrow(use.exprs))
          tmp.var.mat <- cbind(tmp.var.mat, use.exprs[rownames(tmp.var.mat),])
        } 
        
      } else {
        use.exprs <- use.mes[use.samps$SampleID,]
        tmp.var.mat <- use.exprs
      }
      
      auc.vec <- setNames(fin.drugs$auc, fin.drugs$lab_id)[rownames(tmp.var.mat)]
      
      # *Note the stand=T here means that we are having glmnet standardize the Xs here
      drug.run <- get_boot_list(tmp.var.mat, auc.vec, cur.drug$inhibitor[1], T, numIterations = numIterations, numBoots = numBoots)
      drug.run
      
    }else{
      
      NULL
      
    }
    
  }, simplify = FALSE)
  
  # Run for current drug
  drug.res.list <- lapply(run.comb, function(run.list){
    
    # Run in parallel
    res.list <- mclapply(run.list, function(x){  
      
      if (is.null(x)){
        # For debugging
        print('NULL')
        NULL
      }else{
        # Run LASSO and gather data
        tmp.list <- run_lasso(x$auc, x$mat, x$train, x$test, x$boots, x$stand, x$drug)
        tmp.list$index <- x$rep
        tmp.list$preds$index <- x$rep
        tmp.list
        
      }
    }, mc.cores = numCores)  

  })  
  
  # For debugging
  print(names(drug.res.list))
  
  # Now summarize results  
  res.list <- lapply(drug.res.list, function(tmp.list){
    
    pred.dt <- do.call(rbind, lapply(tmp.list, function(y){
      tmp.dt <- y$preds
      tmp.dt$inhibitor <- y$drug
      tmp.dt
    }))
    
    coef.mat <- do.call(cbind, lapply(tmp.list, function(y){
      tmp.mat <- y$coef
      colnames(tmp.mat) <- paste(y$preds$index[1], seq_len(ncol(tmp.mat)), sep="_")
      tmp.mat
    }))
    
    list(preds=pred.dt, coefs=coef.mat, drug=pred.dt$inhibitor[1], datatype="combined")
    
  })
  
  # Coerce into prediction results
  pred.sum <- do.call(rbind, lapply(res.list, "[[", "preds"))
  
  # Get list of coefficient matrices
  coef.list <- lapply(res.list, "[[" ,"coefs")
  
  # Save final results to file
  save(pred.sum, coef.list, file=outFile)
}

