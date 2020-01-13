### Functions needed for the Accelerated Failure Time (AFT) model
### Code obtained from https://www.degruyter.com/view/j/sagmb.2013.12.issue-6/sagmb-2013-0021/Appendix2_corrected.pdf ###

################################################################################ 
# FUNCTIONS 
################################################################################ 
################################################################################ 
# These functions test for covariates effects using an accelerated failure time model
# Then the presence of a point mass at 0 is tested with a LRT test that compares 
# the AFT model to the point mass mixture model.
# Parameter estimates from the AFT and mixture model are provided
# Data are assumed to come from a log normal distribution
# The survival package is required 
################################################################################
################################################################################ 
# Function Arguments
# data - vector of observations to analyze
# covar - matrix of covariates. Must be a matrix and must be numeric. 
################################################################################
################################################################################ 
# AnalyzeMixture is a wrapper function calling functions to fit the AFT model and
# mixture model and return results
# pm.alpha is the desired significance level for testing whether the AFT or mixture
# model fits best. Default is 0.2. The null distribution for this test is a 50:50
# mixture of a chi-square distribution with 1 df and a chi-square distribution
# with 0 df. 
################################################################################
# library(survival)
AnalyzeMixture <- function(data, covar, pm.alpha=0.2){
  censored <- AFT(data, covar)
  mixture <- Mixture(data, covar)
  pmTest <- -2*(censored$loglik[2]-mixture$FullLogl)
  Out <- ifelse(pmTest > 0.5*qchisq(1-pm.alpha, df=1), "Reject", "Dont Reject")
  Cov.Test <- c(Z = censored$table["covar","z"], P.value=censored$table["covar","p"])
  PM.Test <- list(Statistic=pmTest, Reject=Out)
  CensorEstimates <- c(CntlMean=censored$table["(Intercept)","Value"], CaseMean= sum(censored$table[1:2,"Value"]), STD=censored$scale)
  MixtureEstimates <- mixture$Estimates
  ans <- list(CovTest=Cov.Test, PMTest=PM.Test, CensorEstimates=CensorEstimates, MixtureEstimates=MixtureEstimates)
  return(ans) }
################################################################################ 
# Mixture
# This function formats the data for analysis with the mixture models,
# maximizes the likelihood, tests for covariate effects, returns the null and
# full logl and parameter estimates
# As currently written this function only processes 1 covariate with 2 levels
# that must be code (0,1). Can be modified to evaluate multiple covariates
# Values:
# Factor - name of covariate
# Statistic - Likelihood ratio chi-square statistic testing for covariate effect on tau and mu
# P.value - p-value of Statistic. Assumes a 2 df test
# FullLogl - log likelihood of full model
# NullLogl - log likeilhood of null model (tau and mu equal)
# Estimates - CntlMean = log normal mean of baseline group; CntlTau = proportion
# from the log normal distirbution in the baseline group; CaseMean = log normal
# mean in the treatment group; CaseTau = proportion from log normal in treatment
# group; STD = standard deviation (log scale) of the log normal distributions 
################################################################################
Mixture <- function(data, covar){
  # Create binary variable for point mass vs continuous component
  # pm = 0 if in point mass
  pm <- ifelse(data > 0,1,0)
  # FIT FULL MODEL
  cont <- data[data > 0]
  cont.cov <- covar[data > 0, ,drop=FALSE]
  cens <- data[data==0]
  cens.cov <- covar[data==0, , drop=FALSE]
  # Getting starting parameter values
  data[data==0] <- 0.5*min(cont)
  propTrunc <- pnorm(min(log(cont)),mean=mean(log(data)), sd=sd(log(data)))
  propAll <- length(cens)/length(data)
  # OneMinusTau = 1-tau which is the proportion of true zeros and tau is the proportion
  # from the continuous component
  OneMinusTau <- propAll-propTrunc
  # Can get OneMinusTau less than 1. Check for this and substitute one fifth of
  # the observed proportion missing (this is arbitrary)
  if (OneMinusTau < 0) OneMinusTau <- propAll/5
  Tau <- 1-OneMinusTau
  beta0 <- log(Tau/(1-Tau))
  linear <- lm(log(cont)~cont.cov)
  Params <- c(coefficients(linear), beta0, rep(0,ncol(covar)), summary(linear)$sigma)
  out <- try(optim(Params, mixtureLogL, method="L-BFGS-B",obs=cont, cen=cens, obs.covars=cont.cov, cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05), control=list(maxit=1000)),
             silent=TRUE)
  # Check for convergence error
  if (class(out)!="try-error"){ if (out$convergence==1){
    warning("Algorithm did not converge for full likelihood.")
    full.logl <- NA } else {
      
      full.logl <- out$value
      cntl.beta <- out$par[3]
      cntl.mean <- out$par[1]
      cntl.tau <- exp(out$par[3])/(1+exp(out$par[3])) 
      case.mean <- out$par[1]+out$par[2]
      case.beta <- out$par[3]+out$par[4]
      case.tau <- exp(case.beta)/(1+exp(case.beta)) 
      full.logl <- -1*out$value[1]
      std <- out$par[5]
    }
  } else {
    # If there is a boundary error then optimize without the lower bound specified 
    out <- nlminb(Params, mixtureLogL, obs=cont, cen=cens, obs.covars=cont.cov,
    cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05))
    cntl.mean <- out$par[1]
    cntl.beta <- out$par[3]
    cntl.tau <- exp(out$par[3])/(1+exp(out$par[3]))
    case.mean <- out$par[1]+out$par[2] 
    case.beta <- out$par[3]+out$par[4]
    case.tau <- exp(case.beta)/(1+exp(case.beta)) 
    full.logl <- -1*out$objective[1]
    std <- out$par[5]
  }
  Pvalues <- NULL
  Stats <- NULL
  cov.names <- colnames(covar) 
  if (is.null(colnames(covar))){
    cov.names <- as.character(seq(1,ncol(covar),1)) }
  # FIT NULL MODEL
  linear <- lm(log(cont)~1)
  ParamsN <- c(coefficients(linear), beta0, summary(linear)$sigma)
  outN <- try(optim(ParamsN, mixtureLogL, method="L-BFGS-B",obs=cont, cen=cens,
                    obs.covars=NULL,
                    cen.covars=NULL, lower=c(rep(-Inf,length(ParamsN)-1),0.05), control=list(maxit=1000)),
              silent=TRUE)
  if (class(outN)!="try-error"){ if (outN$convergence==1){
    warning("Algorithm did not converge for null likelihood.")
    null.logl <- NA } else {
      null.logl <- -outN$value }
  } else {
    # If there is a boundary error then optimize with nlminb
    outN <- nlminb(ParamsN, mixtureLogL, obs=cont, cen=cens, obs.covars=NULL, cen.covars=NULL, lower=c(rep(-Inf,length(ParamsN)-1),0.05))
    null.logl <- -outN$objective }
  
  X2 <- -2*(null.logl-full.logl) 
  p.value <- 1-pchisq(X2, df=2) 
  Pvalues <- c(Pvalues,p.value) 
  Stats <- c(Stats,X2)
  ans <- list(Factor=cov.names, Statistic=Stats, P.value=Pvalues, FullLogl=full.logl, NullLogl=null.logl, Estimates=c(CntlMean=as.numeric(cntl.mean), CntlTau=cntl.tau,
                                                                                                                      CaseMean=as.numeric(case.mean), CaseTau=case.tau, STD=std)) 
  return(ans)
}
################################################################################ 
# mixtureLogl is the log likelihood for the mixture model that is maximized
# params is a vector of the model parameters with the linear coefficients first
# followed by the logistic coefficients and last is the standard deviation of the
# the log normal distirbution
# This function allows evaluation of multiple covariates
# obs is a vector of the observed (non-zero values)
# cen is a vector of the censored or 0 values (all values are 0)
# obs.covars is a matrix n.obs X p of the covariates of the non-zero observations
# cen.covars is a n.cen X p matrix of the covariates of the zero observations 
################################################################################ 
mixtureLogL <- function(params, obs, cen, obs.covars, cen.covars){
d <- min(obs)
Obs.Covars <- cbind(rep(1,length(obs)), obs.covars) 
Cen.Covars <- cbind(rep(1,length(cen)), cen.covars)
n.param <- ncol(Obs.Covars)
betas <- 1:n.param
gammas <- c((n.param+1):(2*n.param))
sd.loc <- length(params)
mu.obs <- Obs.Covars%*%params[betas]
tau.obs <- 1/(1+exp(-1*(Obs.Covars%*%params[gammas]))) 
logl.obs <- sum(log(tau.obs*(exp((-1*(log(obs)-mu.obs)^2)/(2*(params[sd.loc]^2)))/(obs*sqrt(2*pi)*params[sd.loc]))))

# Changed to avoid error when fitting a model without missing data
if(length(cens) != 0){
mu.cen <- Cen.Covars%*%params[betas]
tau.cen <- 1/(1+exp(-1*(Cen.Covars%*%params[gammas]))) 
logl.cen <- sum(log((1-tau.cen)+tau.cen*pnorm((log(d)-mu.cen)/params[sd.loc]))) 
} else{
  logl.cen <- 0
}

logl <- logl.obs+logl.cen
return(-logl)
}

################################################################################ 
# AFT function is wrapper function for formatting and fitting accelerated
# failure time model with survreg function in survival package 
################################################################################
AFT <- function(data, covar){
  # Make event indicator
  # Z=0 censored, Z=1 event observed
  Z <- ifelse(data > 0,1,0)
  # Put in minimum value for censored values
  d <- min(data[data > 0])
  data[data==0] <- d
  # Create survival object
  surv <- Surv(data,Z,type="left")
  censor <- survreg(surv~covar, dist="lognormal") 
  ans <- summary(censor)
  return(ans)
}

## msnset = peptidesCPTAC
## fits AFT model on all proteins
do_AFT = function(formula, msnset, type_df, group_var = feature,
                 contrasts = NULL, lfc = 0, p.adjust.method = "BH",
                 ## choose parallel_plan =sequential if you don't want parallelisation 
                 #, parallel_plan = multiprocess
                 parallel = TRUE){
  
  if(!(type_df %in% c("conservative", "traceHat"))) {
    stop("Invalid input `type_df`.")
  }
  
  system.time({## can take a while
    if (parallel){
      cl <- makeClusterPSOCK(availableCores())
      plan(cluster, workers = cl)   
    }
    ## future::plan(parallel_plan,gc = TRUE)
    formula <- update(formula, expression ~ . )
    group_var <- enquo(group_var) # group_var = quo(protein)
    # target <- pData(msnset)
    df <- MSnSet2df(msnset, na.rm = FALSE)
    
    # The ATF model needs log-normal data and zeros for missing values
    df$expression <- 2^df$expression
    df$expression[is.na(df$expression)] <- 0
    
    ## Glm adds variable name to levels in catogorical (eg for contrast)
    ## lme4 doesnt do this for random effect, so add beforehand
    ## all_vars <- formula %>% terms %>% delete.response %>% all.vars
    # df = lme4:::findbars(formula) %>% map_chr(all.vars) %>%
    #   purrr::reduce(~{mutate_at(.x,.y,funs(paste0(.y,.)))}, .init=df)

    cat("Fitting AFT models\n")
    ## select only columns needed for fitting
    df_prot <- select(df, !!group_var, one_of(all.vars(formula))) %>%
      group_by(!!group_var) %>% nest %>%
      mutate(model = furrr::future_map(data,~{gc();try(do_AFTfit(.x, formula))}))
    ## Return also failed ones afterward
    df_prot_failed <- filter(df_prot, map_lgl(model,~{class(.x) != "lmerMod"}))
    df_prot <- filter(df_prot, map_lgl(model, ~{class(.x)=="lmerMod"}))
    
    if(nrow(df_prot) == 0) {print("No models could be fitted"); return()}
    
    df_prot <- mutate(df_prot,
                      ## get trace hat df for squeezeVar
                      df = map_dbl(model, ~getDf(.x)),
                      sigma = map_dbl(model,~{MSqRob::getSigma(.x)}))
    ## df_protein = map_dbl(model,~{MSqRob::getDf(.x)}))
    ## Squeeze variance
    squeezeHlp <- squeezeVarRob(df_prot$sigma^2, df_prot$df, robust = TRUE) # MSqRob::squeezeVarRob
    df_prot <- mutate(df_prot,
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
      unnest %>%
      group_by(contrast) %>%
      mutate(qvalue = p.adjust(pvalue, method = p.adjust.method)) %>%
      group_by(!!group_var) %>% nest(.key = contrasts) %>%
      left_join(df_prot,.)
  }
  ) %>% print
  model = bind_rows(df_prot,df_prot_failed)
  result = model %>% select(!!group_var, contrasts) %>% filter(map_lgl(contrasts,~{!is.null(.x)})) %>% unnest()
  if (parallel) stopCluster(cl)
  list(model = model, result = result)
}

df <- exprs
do_AFTfit <- function(df, form, pd){
  
  data <- unlist(df[1,])
  dataframe <- cbind(expression = data, pd)
  
  covar <- model.matrix(form, data = dataframe)[,2,drop=FALSE]
  
  # If there is no response variable in the formula, throw an error
  if(attr(terms(form),"response") != 1){
    stop("The formula should contain a response variable!")
  }
  # Now that we know that there is a response variable in our formula, we know that it is the first term.
  data <- df[, all.vars(form)[1], drop = TRUE]

  AnalyzeMixture(data = data, covar = covar)
  
}


