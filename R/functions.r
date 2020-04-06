################################################################
## load package or install automatically if not installed yet ##
################################################################
load_lib = function(lib){
  err <- try(library(lib, character.only = TRUE), silent = TRUE)
  if (class(err) == "try-error") {
    install.packages(lib, repos = "https://cloud.r-project.org")
    library(lib, character.only = TRUE)
  }
}

# Install BioConductor if not installed yet
if (!requireNamespace("BiocManager")){
  install.packages("BiocManager")
  BiocManager::install()
}

load_bioConductor <- function(lib, path){
  err <- try(library(lib, character.only = TRUE), silent = TRUE)
  if (class(err) == "try-error") {
    BiocManager::install(path)
    library(lib, character.only = TRUE)
  }
}

load_lib_GitHub <- function(lib, path){
  err <- try(library(lib, character.only = TRUE), silent = TRUE)
  if (class(err) == "try-error") {
    devtools::install_github(path)
    library(lib, character.only = TRUE)
  }
}

load_lib("devtools")
load_bioConductor("preprocessCore", "preprocessCore")
load_lib("remotes")
load_bioConductor("MSstats", "MSstats")
load_lib("tidyverse")
load_lib("limma")
load_lib("Matrix")
load_lib("lme4")
load_lib("survival")
load_lib("corpcor")
load_lib("zoo")
load_lib("colorspace")
load_lib("readxl")
load_lib("dplyr")
load_lib("here")
load_bioConductor("xcms", "xcms")
load_bioConductor("MSnbase", "lgatto/MSnbase")

# Install and load furrr
load_lib_GitHub("furrr", "DavisVaughan/furrr")

# Install and load stageR version 1.3.29
load_lib_GitHub("stageR", "statOmics/stageR@6595e4412d040890c187ae6fd962ae20b6c72942")

rename <- dplyr::rename
count <- dplyr::count
select <- dplyr::select
slice <- dplyr::slice
desc <- dplyr::desc

# Install and load MSqRob 0.7.6
load_lib_GitHub("MSqRob", "statOmics/MSqRob@MSqRob0.7.6")

# Install and load msqrobsum
load_lib_GitHub("msqrobsum", "statOmics/MSqRobSum")

load_lib_GitHub("proDA", "const-ae/proDA")

#######################
## Data manipulation ##
#######################

MSnSet2df = function(msnset, na.rm = TRUE){
  ## Converts Msnset to a tidy dataframe
  ## Always creates feature and vector column so these shouldn't be defined by user.
  ## convenient for downstream analysis steps.
  if(any(c("sample", "feature", "expression") %in% c(colnames(fData(msnset)),colnames(pData(msnset))))){
    stop("Column names in the \"fData\" or \"pData\" slot of the \"msnset\" object cannot be named
         \"sample\", \"feature\" or \"expression\". Please rename these columns before running the analysis.")
  }

  dt <- as.data.frame(Biobase::exprs(msnset)) %>% mutate(feature = rownames(.)) %>%
    gather(sample, expression, - feature, na.rm = na.rm)
  dt <- fData(msnset) %>% mutate(feature = rownames(.)) %>% left_join(dt,. , by = "feature")
  dt <- pData(msnset) %>% mutate(sample = rownames(.)) %>% left_join(dt,. , by = "sample")
  as_tibble(dt)
}

############################################################
## Empirical Bayes variance and overdispersion estimation ##
############################################################

squeezeVarRob <- function (var, df, covariate = NULL, robust = FALSE, winsor.tail.p = c(0.05, 0.1), min_df = 1, allow_underdispersion = TRUE) 
{
  n <- length(var)
  if (n == 0) 
    stop("var is empty")
  if (sum(!is.na(var)) == 1) {
    return(list(var.post = var, var.prior = var, df.prior = df))
  }
  if (length(df) == 1) {
    df <- rep.int(df, n)
  }  else {
    if (length(df) != n) 
      stop("lengths differ")
  }
  if (robust) {
    fit <- MSqRob:::fitFDistRobustly_LG(var, df1 = df, covariate = covariate, 
                                        winsor.tail.p = winsor.tail.p, min_df = min_df)
    df.prior <- fit$df2.shrunk
  }  else {
    fit <- MSqRob:::fitFDist_LG(var, df1 = df, covariate = covariate, 
                                min_df = min_df)
    df.prior <- fit$df2
  }
  var.prior <- fit$scale
  
  var.post <- MSqRob:::.squeezeVarRob(var = var, df = df, var.prior = var.prior, 
                                      df.prior = df.prior)
  
  if(!allow_underdispersion){
    var.post[var.post<1] <- 1
  }
  
  list(df.prior = df.prior, var.prior = fit$scale, var.post = var.post)
}

##################
## Mixed models ##
##################

setGeneric (
  name = "getBetaB",
  def = function(model,...){standardGeneric("getBetaB")}
)

.getBetaBMermod <- function(model) {
  betaB <- c(as.vector(lme4::getME(model,"beta")),as.vector(lme4::getME(model,"b")))
  names(betaB) <- c(colnames(lme4::getME(model,"X")),rownames(lme4::getME(model,"Zt")))
  betaB
}
setMethod("getBetaB", "lmerMod", .getBetaBMermod)

.getBetaBGlm <- function(model) 
  model$coefficients

setGeneric (
  name = "getVcovBetaBUnscaled",
  def = function(model,...){standardGeneric("getVcovBetaBUnscaled")}
)

setMethod("getBetaB", "lm", .getBetaBGlm)
setMethod("getBetaB", "glm", .getBetaBGlm)

.getVcovBetaBUnscaledMermod <- function(model){
  ## TODO speed up (see code GAM4)
  p <- ncol(lme4::getME(model,"X"))
  q <- nrow(lme4::getME(model,"Zt"))
  Ct <- rbind2(t(lme4::getME(model,"X")),lme4::getME(model,"Zt"))
  Ginv <- Matrix::solve(Matrix::tcrossprod(lme4::getME(model,"Lambda"))+Matrix::Diagonal(q,1e-18))
  vcovInv <- Matrix::tcrossprod(Ct)
  vcovInv[((p+1):(q+p)),((p+1):(q+p))] <- vcovInv[((p+1):(q+p)),((p+1):(q+p))]+Ginv

  Matrix::solve(vcovInv)
}

setMethod("getVcovBetaBUnscaled", "lmerMod", .getVcovBetaBUnscaledMermod)

.getVcovBetaBUnscaledGlm <- function(model)
  ## cov.scaled is scaled with the dispersion, "cov.scaled" is without the dispersion!
  ## MSqRob::getSigma is needed because regular "sigma" function can return "NaN" when sigma is very small!
  ## This might cause contrasts that can be estimated using summary() to be NA with our approach!
  summary(model)$cov.scaled/MSqRob::getSigma(model)^2

.getVcovBetaBUnscaledLm <- function(model)
  ## cov.scaled is scaled with the dispersion, "cov.scaled" is without the dispersion!
  ## MSqRob::getSigma is needed because regular "sigma" function can return "NaN" when sigma is very small!
  ## This might cause contrasts that can be estimated using summary() to be NA with our approach!
  summary(model)$cov.unscaled

setMethod("getVcovBetaBUnscaled", "glm", .getVcovBetaBUnscaledGlm)
setMethod("getVcovBetaBUnscaled", "lm", .getVcovBetaBUnscaledLm)

## Estimate pvalues contrasts
contrast_helper <- function(formula, msnset, contrast = NULL){
  ## Gives back the coefficients you can use to make contrasts with given the formula and dataset
  ## If a factor variable is specified (that is present in the formula) all the possible contrasts
  ## within this variable are returned
  contrast <- enquo(contrast) ;#contrast = quo(condition)
  df <- MSnSet2df(msnset)
  all_vars <- formula %>% terms %>% delete.response %>% all.vars
  names(all_vars) <- all_vars
  df[,all_vars] <- map2_dfr(all_vars,df[,all_vars],paste0)
  coefficients <- c("(Intercept)", df %>% select(all_vars) %>% unlist %>% unique %>% as.character)
  if (!is.null(rlang::quo_get_expr(contrast))) {
    c <- pull(df, !!contrast) %>% unique %>% sort %>% as.factor
    comp <- combn(c, 2, simplify = FALSE)
    ## condIds = map(comp,~which( coefficients %in% .x))
    ## L = rep(0,length(coefficients))
    ## L = sapply(condIds,function(x){L[x]=c(-1,1);L})
    ## rownames(L) = coefficients
    ## colnames(L) = map_chr(comp, ~paste(.x,collapse = "-"))
    condIds <- map(comp, ~which(coefficients %in% .x))
    L <- rep(0,nlevels(c))
    L <- sapply(comp,function(x){L[x]=c(-1,1);L})
    rownames(L) <- levels(c)
    colnames(L) <- map_chr(comp, ~paste(rev(.x),collapse = "-"))
    L
  } else coefficients
}

setGeneric (
  name= "getXLevels",
  def=function(model,...){standardGeneric("getXLevels")}
)

.getXLevelsGlm <- function(model)
  map2(names(model$xlevels), model$xlevels, paste0) %>% unlist

setMethod("getXLevels", "glm", .getXLevelsGlm)
setMethod("getXLevels", "lm", .getXLevelsGlm)

.getXLevelsMermod <- function(model)
  c(lme4::fixef(model) %>% names, lme4::getME(model,"flist") %>% map(levels) %>% unlist %>% unname)

setMethod("getXLevels", "lmerMod", .getXLevelsMermod)

contEst <- function(model, contrasts, var, df, lfc = 0){
  
  ## Contrasts should have a name
  if (is.null(colnames(contrasts))) stop("Contrast matrix should have column names.")
  
  betaB <- getBetaB(model)
  vcov <- getVcovBetaBUnscaled(model)
  coefficients <- names(betaB)
  id <- coefficients %in% rownames(contrasts)

  coefficients <- coefficients[id]
  vcov <- vcov[id,id]
  betaB <- betaB[id]

  xlevels <- getXLevels(model)
  id <- !apply(contrasts,2,function(x){any(x[!(rownames(contrasts) %in% xlevels)] !=0)})
  contrasts <- contrasts[coefficients, id, drop = FALSE]
  ## If no contrasts could be found, terminate
  if (is.null(colnames(contrasts))) return(new_tibble(list(), nrow = 0))

  se <- sqrt(diag(t(contrasts)%*%vcov%*%contrasts)*var)
  logFC <- (t(contrasts)%*%betaB)[,1]
  
  ### Addition to allow testing against another log FC (lfc)
  ### See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2654802/
  
  lfc <- abs(lfc)
  aest <- abs(logFC)
  
  Tval <- setNames(rep(0, length(logFC)),names(se))
  tstat.right <- (aest - lfc)/se
  tstat.left <- (aest + lfc)/se
  pval <- pt(tstat.right, df = df, lower.tail = FALSE) +
    pt(tstat.left, df = df, lower.tail = FALSE)
  tstat.right <- pmax(tstat.right, 0)
  fc.up <- (logFC >= lfc)
  fc.up[is.na(fc.up)] <- FALSE
  fc.down <- (logFC < -lfc)
  fc.down[is.na(fc.down)] <- FALSE
  Tval[fc.up] <- tstat.right[fc.up]
  Tval[fc.down] <- -tstat.right[fc.down]
  Tval[is.na(logFC)] <- NA
  
  
  new_tibble(list(contrast = colnames(contrasts),
                  logFC = logFC,
                  se = se,
                  t = Tval,
                  df = rep(df, length(se)),
                  pvalue = pval),
             nrow = ncol(contrasts))
}

do_lmerfit <- function(df, form, nIter = 10, tol = 1e-6, control = lme4::lmerControl(calc.derivs = FALSE)){
  fit <- lmer(form, data = df, control = control)
  ##Initialize SSE
  res <- resid(fit)
  ## sseOld=sum(res^2)
  sseOld <- fit@devcomp$cmp["pwrss"]
  while (nIter > 0){
    nIter = nIter-1
    fit@frame$`(weights)` <- MASS::psi.huber(res/(mad(res)))
    fit <- refit(fit)
    res <- resid(fit)
    ## sse=sum(res^2)
    sse <- fit@devcomp$cmp["pwrss"]
    if(abs(sseOld-sse)/sseOld <= tol) break
    sseOld <- sse
  }
  return(fit)
}

calculate_df <- function(df, model, vars){
  ## Get all the variables in the formula that are not defined in vars
  form <- attributes(model@frame)$formula
  vars_formula <- all.vars(form)
  vars_drop <- vars_formula[!vars_formula %in% vars]
  ## Sum of number of columns -1 of Zt mtrix of each random effect that does not involve a variable in vars_drop
  mq <- lme4::getME(model,"q_i")
  id <- !map_lgl(names(mq),~{any(stringr::str_detect(.x,vars_drop))})
  p <- sum(mq[id]) - sum(id)
  ## Sum of fixed effect parameters that do not involve a variable in vars_drop
  mx <- lme4::getME(model,"X")
  id <- !map_lgl(colnames(mx),~{any(stringr::str_detect(.x,vars_drop))})
  p <- p + sum(id)

  ## n is number of sample because 1 protein defined per sample
  n <- n_distinct(df$sample)
  n-p
}

## msnset = peptidesCPTAC
## fits mixed model on all proteins
do_mm <- function(formula, msnset, type_df, group_var = feature,
                 contrasts = NULL, lfc = 0, p.adjust.method = "BH", max_iter = 20L,
                 ## choose parallel_plan =sequential if you don't want parallelisation 
                 control = lme4::lmerControl(calc.derivs = FALSE)#, parallel_plan = multiprocess
                 , parallel = FALSE, n_cores = max(c(1,future::availableCores()-4))){

                   if(!(type_df %in% c("conservative", "traceHat"))) {
                     stop("Invalid input `type_df`.")
                   }
                   
  system.time({## can take a while
    if (parallel){
      cl <- makeClusterPSOCK(n_cores)
      plan(cluster, workers = cl)   
    }
    ## future::plan(parallel_plan, gc = TRUE)
    formula <- update(formula, expression ~ . )
    group_var <- enquo(group_var) # group_var = quo(protein)
    # target <- pData(msnset)
    df <- MSnSet2df(msnset)

    ## Glm adds variable name to levels in categorical (eg for contrast)
    ## lme4 doesn't do this for random effect, so add beforehand
    ## all_vars <- formula %>% terms %>% delete.response %>% all.vars
    df <- lme4:::findbars(formula) %>% map_chr(all.vars) %>%
      purrr::reduce(~{mutate_at(.x,.y,funs(paste0(.y,.)))}, .init=df)

    cat("Fitting mixed models\n")
    ## select only columns needed for fitting
    df_prot <- select(df, !!group_var, one_of(all.vars(formula))) %>%
      group_by(!!group_var) %>% nest %>%
      mutate(model = furrr::future_map(data,~{gc();try(do_lmerfit(.x, formula, nIter = max_iter,
                                                                  control = control), silent = TRUE)}))
    ## Return also failed ones afterward
    df_prot_failed <- filter(df_prot, map_lgl(model,~{class(.x) != "lmerMod"}))
    df_prot <- filter(df_prot, map_lgl(model, ~{class(.x)=="lmerMod"}))

    if(nrow(df_prot) == 0) {print("No models could be fitted"); return()}

    df_prot <- mutate(df_prot,
                      ## get trace hat df for squeezeVar
                      df = map_dbl(model, ~MSqRob::getDf(.x)),
                      sigma = map_dbl(model,~{MSqRob::getSigma(.x)}))
    ## df_protein = map_dbl(model,~{MSqRob::getDf(.x)}))
    ## Squeeze variance
    squeezeHlp <- squeezeVarRob(df_prot$sigma^2, df_prot$df, robust = TRUE) # MSqRob::squeezeVarRob
    df_prot <- mutate(ungroup(df_prot),
                      df_prior = squeezeHlp$df.prior,
                      var_prior = squeezeHlp$var.prior,
                      var_post = squeezeHlp$var.post,
                      sigma_post = sqrt(var_post),
                      df_post = df + df_prior)

    if (type_df == "conservative"){
    ## Calculate df on protein level, assumption is that there is only one protein value/run,
      df_prot <- mutate(df_prot,
                        df_protein = map2_dbl(data, model,~calculate_df(.x,.y, vars = colnames(pData(msnset)))))
    } else if (type_df == "traceHat"){
    # Alternative: MSqRob implementation with trace(Hat):
    df_prot <- df_prot %>% mutate(df_protein = df_post)
    }
   
    ## Calculate fold changes and p values for contrast
    cat("Estimating p-values contrasts\n")
    ## If only within shrinkage, use residual squeezed df, else use conservative df
    df_prot <- df_prot %>% 
      transmute(!!group_var,
                contrasts = furrr::future_pmap(list(model = model, contrasts = list(contrasts), var = var_post,
                                             df = df_protein, lfc = lfc), contEst))  %>%
      ## Calculate qvalues BH
      unnest(cols = contrasts) %>%
      group_by(contrast) %>%
      mutate(qvalue = p.adjust(pvalue, method = p.adjust.method)) %>%
      group_by(!!group_var) %>% nest(contrasts = c(contrast, logFC, se, t, df, pvalue, qvalue), .key = contrasts) %>%
      left_join(df_prot,.)
  }
  ) %>% print
  model <- bind_rows(df_prot,df_prot_failed)
  result <- model %>% select(!!group_var, contrasts) %>% filter(map_lgl(contrasts, ~{!is.null(.x)})) %>% unnest(cols = contrasts)
  if (parallel) stopCluster(cl)
  return(list(model = model, result = result))
}

##########################
## Quasibinomial models ##
##########################

do_count_groups <- function(msnset, group_var = protein, keep_fData_cols = NULL){
  ## Very general because feature(names) and sample(names) variable will be found in every msnset
  ## Can also be used for multiple rounds of normalization, e.g. first from features to peptides, then from peptides to proteins
  
  group_var <- enquo(group_var) #group_var <- quo(protein)
  dt <- MSnSet2df(msnset)
  
  if("group_count" %in% colnames(dt)){
    stop("Column names in the \"fData\" or \"pData\" slot of the \"msnset\" object cannot be named 
         \"group_count\". Please rename this column before running the analysis.")
  }
  
  dt_s <- filter(dt, !is.na(expression)) %>% group_by(sample, !!group_var) %>%
    summarise(count = n_distinct(feature)) %>% distinct
  ddd  <- spread(dt_s, sample, count) %>% ungroup
  exprs_data <- dplyr::select(ddd, - !!group_var) %>% as.matrix
  exprs_data[is.na(exprs_data)] <- 0
  # DO NOT DO THIS: names may not be ordered!: colnames(exprs_data) <- sampleNames(msnset)
  exprs_data <- exprs_data[ , sampleNames(msnset)]
  rownames(exprs_data) <- as.character(pull(ddd, !!group_var))
  
  fd1 <- dplyr::select(ddd, !!group_var)
  # Select the group variable and all variables you want to keep
  fd2 <- msnset %>% fData %>% dplyr::select(c(rlang::quo_text(group_var), keep_fData_cols)) %>% distinct(!!group_var, .keep_all = TRUE)
  
  if(nrow(fd1) != nrow(fd2 %>% distinct)){
    stop("Values in the \"group_var\" column can only correspond to a single value in the \"keep_fData_cols\" column.")
  }
  
  # fd <- fd2 %>% slice(match(unlist(fd1), !!group_var))
  fd <- fd2 %>% slice(match(unlist(fd1), pull(fd2, !!group_var)))
  
  # Either in a separate slot in MSnBase or added to the fData slot
  dt_p <- filter(dt, !is.na(expression)) %>% group_by(!!group_var) %>%
    summarise(group_count = n_distinct(feature)) %>% distinct
  
  fd$group_count <- dt_p$group_count
  
  # Extra safety check
  if(all(pull(dt_p, !!group_var) != pull(fd, !!group_var))){
    stop("Something went wrong with sorting the proteins. Please contact the creators of this package.")
  }
  
  fd <- as.data.frame(fd)
  rownames(fd) <- as.character(pull(ddd, !!group_var))
  
  pData(msnset)$numIds <- colSums(exprs_data)
  pData(msnset)$odds <- pData(msnset)$numIds/(sum(fd$group_count)-pData(msnset)$numIds)
  pData(msnset)$lodds <- log(pData(msnset)$odds)
  
  out <- MSnSet(exprs_data, fData = AnnotatedDataFrame(fd) , pData = pData(msnset))
  
  return(out)
}

stripGlmLR <- function(cm) {
  # cm$y = c()
  cm$model = c()
  
  # cm$residuals = c()
  # cm$fitted.values = c()
  # cm$effects = c()
  # cm$qr$qr = c()
  # cm$linear.predictors = c()
  # cm$weights = c()
  # cm$prior.weights = c()
  cm$data = c()
  
  
  # cm$family$variance = c()
  # cm$family$dev.resids = c()
  # cm$family$aic = c()
  # cm$family$validmu = c()
  # cm$family$simulate = c()
  # attr(cm$terms,".Environment") = c()
  # attr(cm$formula,".Environment") = c()
  
  cm
}

## msnset = peptidesCPTAC
## fits mixed model on all proteins
## default: limma with dispersion parameter smaller than 1 set to 1 (i.e. no underdispersion allowed)
do_glm <- function(formula = ~ condition + lab, msnset,  group_var = feature, familyfun = quasibinomial(link = "logit"),
                  contrasts = NULL, add_val = 0.1, contFun = "contEst", p.adjust.method = "BH", squeezeVar = TRUE, allow_underdispersion = FALSE, parallel_plan = multiprocess){
  ## choose parallel_plan = sequential if you don't want parallelisation 
  
  system.time({## can take a while
    
    plan(parallel_plan)
    
    e <- new.env(parent=environment(formula))
    assign("add_val", add_val, envir=e)
    environment(formula) <- e
    
    if(contFun == "contEst"){
      pvals <- "pvalue"
    } else if(contFun == "contEstQLT"){
      pvals <- c("pvalue.LRT", "pvalue.F")
    }
    
    contFun <- rlang::parse_quosure(contFun)
    
    formula <- update(formula, cbind(expression, group_count-expression) + add_val ~ . )
    group_var <- enquo(group_var) # group_var <- quo(protein)
    # target <- pData(msnset)
    df <- MSnSet2df(msnset)
    cat("Fitting glm models\n")
    df_prot <- select(df, !!group_var, one_of(all.vars(formula)[-3])) %>% # -3 is to remove "add_val"
      group_by(!!group_var) %>% nest %>%
      mutate(model = future_map(data,~{gc();try(glm(formula, data = .x, family=familyfun) %>% stripGlmLR, silent = TRUE)})) #future_map
    ## return also failed ones afterwards
    df_prot_failed <- filter(df_prot, map_lgl(model,~{!("glm" %in% class(.x))}))
    df_prot <- filter(df_prot, map_lgl(model,~{("glm" %in% class(.x))}))
    if(nrow(df_prot) == 0 ) {print("No models could be fitted"); return()}
    
    ## Calculate df for glm model
    df_prot <- df_prot %>% mutate(df_protein = map_dbl(model, ~MSqRob::getDf(.x)))
    
    ## Squeeze variance
    df_prot <- mutate(df_prot,
                      sigma = map_dbl(model,~{MSqRob::getSigma(.x)}),
                      df = map_dbl(model,~{MSqRob::getDf(.x)}))
    squeezeHlp <- squeezeVarRob(df_prot$sigma^2,df_prot$df, robust = TRUE, allow_underdispersion = allow_underdispersion) # MSqRob::squeezeVarRob
    df_prot <- mutate(ungroup(df_prot),
                      df_prior = squeezeHlp$df.prior,
                      var_prior = squeezeHlp$var.prior,
                      var_post = squeezeHlp$var.post,
                      sigma_post = sqrt(var_post),
                      df_post = df + df_prior)
    
    if(squeezeVar){
      vars <- quo(var_post)
      dfs <- quo(df_post)
    } else{
      df_prot <- mutate(df_prot, var = sigma^2)
      vars <- quo(var)
      dfs <- quo(df)
    }
    
    ## Calculate fold changes and p values for contrast
    cat("Estimating p-values contrasts\n")
    df_prot <- df_prot %>%
      transmute(!!group_var,
                contrasts = future_pmap(list(model = model, contrasts = list(contrasts), var = !!vars,
                                             df = !!dfs), !!contFun))  %>%
      ## There must be at least one contrast, otherwise "unnest" doesn't work...
      ## filter(unlist(lapply(contrasts, nrow))!=0) %>%
      filter(map_int(contrasts, nrow)!=0) %>%
      ## Calculate qvalues BH
      unnest(cols = contrasts) %>%
      group_by(contrast) %>%
      mutate_at(.vars = pvals, .funs = funs(qvalue = p.adjust(., method = p.adjust.method))) %>%
      group_by(!!group_var) %>% nest(contrasts = c(contrast, logFC, se, t, df, pvalue, qvalue),.key = contrasts) %>% #.key added for backwards compatibility
      left_join(df_prot,.)
  }) %>% print
  model <- bind_rows(df_prot, df_prot_failed)
  result <- model %>% dplyr::select(!!group_var,contrasts) %>% filter(map_lgl(contrasts,~{!is.null(.x)})) %>% unnest(cols = contrasts) %>% rename(logOR = logFC)
  return(list(model = model, result = result))
}


