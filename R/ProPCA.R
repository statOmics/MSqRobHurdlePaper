
# Function obtained from https://www.mcponline.org/content/suppl/2010/09/07/M110.002774.DC1
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

do_ProPCA <- function(peptideset, logSCset, group_var = protein, keep_fData_cols = NULL){
  ## Very general because feature(names) and sample(names) variable will be found in every peptideset
  ## Can also be used for multiple rounds of normalization, e.g. first from features to peptides, then from peptides to proteins
  
  group_var <- enquo(group_var) #group_var = quo(protein)
  
  fd <- peptideset %>% fData %>% dplyr::select(c(rlang::quo_text(group_var), keep_fData_cols)) %>% distinct(!!group_var, .keep_all = TRUE)
  
  if((fd %>% pull(!!group_var) %>% unique %>% length) != nrow(fd %>% distinct)){
    stop("Values in the \"group_var\" column can only correspond to a single value in the \"keep_fData_cols\" column.")
  }
  
  system.time({
    
    exprs <- map_dfr(fd$protein,~{
      
      frame <- exprs(peptideset)[fData(peptideset) %>% pull(!!group_var) == .x,] %>% t %>% cbind(logSCset = exprs(logSCset)[fData(logSCset)$protein == .x, ] %>% unlist,.)
      ProPCA(frame) %>% t %>% data.frame
      
    }) %>% as.matrix
    
    rownames(fd) <- fd %>% pull(!!group_var)
    rownames(exprs) <- rownames(fd)
    
    
    MSnSet(exprs(peptideset), fData = AnnotatedDataFrame(fData(peptideset)) , pData = pData(peptideset))
    
    out <- MSnSet(exprs, fData = AnnotatedDataFrame(fd) , pData = pData(peptideset))
    
    return(out)
    
  })
}

### OLD: ###
# do_ProPCA = function(msnset, logSC, group_var = protein, keep_fData_cols = NULL){
#   ## Very general because feature(names) and sample(names) variable will be found in every msnset
#   ## Can also be used for multiple rounds of normalization, e.g. first from features to peptides, then from peptides to proteins
#   
#   group_var <- enquo(group_var) #group_var = quo(protein)
#   
#   fd <- msnset %>% fData %>% dplyr::select(c(rlang::quo_text(group_var), keep_fData_cols)) %>% distinct(!!group_var, .keep_all = TRUE)
#   
#   if((fd %>% pull(!!group_var) %>% unique %>% length) != nrow(fd %>% distinct)){
#     stop("Values in the \"group_var\" column can only correspond to a single value in the \"keep_fData_cols\" column.")
#   }
#   
#   system.time({
#   
#   exprs <- map_dfr(fd$protein,~{
# 
#     frame <- exprs(msnset)[fData(msnset) %>% pull(!!group_var) == .x,] %>% t %>% cbind(logSC = logSC[logSC$Proteins == .x, -1] %>% unlist,.)
#     ProPCA(frame) %>% t %>% data.frame
#     
#   }) %>% as.matrix
#   
#   rownames(fd) <- fd %>% pull(!!group_var)
#   rownames(exprs) <- rownames(fd)
# 
#   
#   MSnSet(exprs(msnset), fData = AnnotatedDataFrame(fData(msnset)) , pData = pData(msnset))
#   
#   out <- MSnSet(exprs, fData = AnnotatedDataFrame(fd) , pData = pData(msnset))
# 
#   return(out)
#   
#   })
# }