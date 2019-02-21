
ProPCA <- function(data, max.iter=1000, err=1e-10){
	data <- as.matrix(data)

	Y <- Y.0 <- Y.1 <- Y.center <- data

	for(jj in 1:ncol(data)) Y[is.na(data)[,jj],jj] <- 0

	ct <- 0
	eps <- Inf
	eps.track <- c()
	obj.0 <- sum(as.vector(Y),na.rm=T)
	
	while(ct < max.iter & eps > err){
		Y.mean <- colMeans(Y)
		for(jj in 1:ncol(Y))
			Y.center[,jj] <- (Y[,jj] - Y.mean[jj])

		v <- eigen(Y.center %*% t(Y.center))$vectors[,1]
		Y.0 <- v%*%t(v)%*%Y.center
		for(jj in 1:ncol(Y.0))
			Y.0[,jj] <- Y.0[,jj] + Y.mean[jj]
	
		eps <- sum(abs(Y.0 - Y.1),na.rm=T)
		Y.1 <- Y.0
		Y[is.na(data)] <- Y.0[is.na(data)]
		ct <- ct + 1
	}

	val <- predict(lm(data[,1]~v))
	val
}

do_ProPCA = function(msnset, logSC, group_var = protein, keep_fData_cols = NULL){
  ## Very general because feature(names) and sample(names) variable will be found in every msnset
  ## Can also be used for multiple rounds of normalization, e.g. first from features to peptides, then from peptides to proteins
  
  group_var <- enquo(group_var) #group_var = quo(protein)
  
  fd <- msnset %>% fData %>% dplyr::select(c(rlang::quo_text(group_var), keep_fData_cols)) %>% distinct(!!group_var, .keep_all = TRUE)
  
  if((fd %>% pull(!!group_var) %>% unique %>% length) != nrow(fd %>% distinct)){
    stop("Values in the \"group_var\" column can only correspond to a single value in the \"keep_fData_cols\" column.")
  }
  
  system.time({
  
  exprs <- map_dfr(fd$protein,~{

    frame <- exprs(msnset)[fData(msnset) %>% pull(!!group_var) == .x,] %>% t %>% cbind(logSC = logSC[logSC$Proteins == .x, -1] %>% unlist,.)
    ProPCA(frame) %>% t %>% data.frame
    
  }) %>% as.matrix
  
  rownames(fd) <- fd %>% pull(!!group_var)
  rownames(exprs) <- rownames(fd)

  
  MSnSet(exprs(msnset), fData = AnnotatedDataFrame(fData(msnset)) , pData = pData(msnset))
  
  out <- MSnSet(exprs, fData = AnnotatedDataFrame(fd) , pData = pData(msnset))

  return(out)
  
  })
}

#
# load log(sc) and log(PPA) data for P00167 and run ProPCA (first 2 columns
# of P00167 contain sample ID and relative abundance -- discard these before 
# running ProPCA)
#

# P00167 <- read.table(paste0(wd,"/R/P00167.tsv"),sep="\t",header=TRUE)
# ProPCA.P00167 <- ProPCA(P00167[,-c(1:2)])
# 
# #
# # plot of log(rel. abund.) v. ProPCA
# #
# 
# ProPCA.P00167
# plot(ProPCA.P00167, log(P00167$rel.abund))
# cor(ProPCA.P00167, log(P00167$rel.abund))


do_lm = function(formula, msnset, group_var = protein,
                 contrasts = NULL, lfc = 0, p.adjust.method = "BH",
                 ## choose parallel_plan =sequential if you don't want parallelisation 
                 #, parallel_plan = multiprocess
                 parallel = TRUE){
 
  system.time({## can take a while
    if (parallel){
      cl <- makeClusterPSOCK(availableCores())
      plan(cluster, workers = cl)   
    }
    ## future::plan(parallel_plan,gc = TRUE)
    formula <- update(formula, expression ~ . )
    group_var <- enquo(group_var) # group_var = quo(protein)
    # target <- pData(msnset)
    df <- MSnSet2df(msnset)
    
    cat("Fitting linear models\n")
    ## select only columns needed for fitting
    df_prot <- select(df, !!group_var, one_of(all.vars(formula))) %>%
      group_by(!!group_var) %>% nest %>%
      mutate(model = furrr::future_map(data,~{gc();try(lm(formula, data = .x))}))
    
    ## Return also failed ones afterward
    df_prot_failed <- filter(df_prot, map_lgl(model,~{class(.x) != "lm"}))
    df_prot <- filter(df_prot, map_lgl(model, ~{class(.x)=="lm"}))
    
    if(nrow(df_prot) == 0) {print("No models could be fitted"); return()}
    
    df_prot <- mutate(df_prot,
                      ## get trace hat df for squeezeVar
                      df_protein = map_dbl(model, ~getDf(.x)),
                      var = map_dbl(model,~{MSqRob::getSigma(.x)^2}))
    
    
    ## Calculate fold changes and p values for contrast
    cat("Estimating p-values contrasts\n")
    ## If only within shrinkage, use residual squeezed df, else use conservative df
    df_prot <- df_prot %>% 
      transmute(!!group_var,
                contrasts = furrr::future_pmap(list(model = model, contrasts = list(contrasts), var = var,
                                                    df = df_protein, lfc = lfc), contEst))  %>%
      ## Calculate qvalues BH
      unnest %>%
      group_by(contrast) %>%
      mutate(qvalue = p.adjust(pvalue, method = p.adjust.method)) %>%
      group_by(!!group_var) %>% nest(.key = contrasts) %>%
      left_join(df_prot,.)
  }
  ) %>% print
  model = bind_rows(df_prot,df_prot_failed)
  result = model %>% dplyr::select(!!group_var, contrasts) %>% filter(map_lgl(contrasts,~{!is.null(.x)})) %>% unnest()
  if (parallel) stopCluster(cl)
  list(model = model, result = result)
}



