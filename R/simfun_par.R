### Scripts obtained from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4061266/ ###


my.simSigma <- function(p, pattern=c("random","decay"),diag=NULL, pCor=0, rho1=0, rho2=0, sd.rho2=0)
{    
  pattern=pattern[1]
  if (pattern=="random"){ 
    sigma = matrix(rho1,nrow=p,ncol=p)
    while (!is.positive.definite(sigma)){
      sigma = matrix(rho1,nrow=p,ncol=p)
      if (is.null(diag)){
        diag(sigma) <- 1
      } else {
        diag(sigma) = diag
      }
      if (pCor>0){
        for (i in 1:(pCor-1)) {
          for (j in (i+1):pCor) {
            rr <- rnorm(1, rho2, sd.rho2)
            if (abs(rr)>0.8) rr = sign(rr)*0.8
            sigma[i, j] = sigma[j,i] = rr*sqrt(diag(sigma)[i]*diag(sigma)[j])
          }}
      }
    }
  } else if (pattern=="decay"){
    sigma <- matrix(0,p,p)
    diag(sigma) <- diag
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        sigma[i,j] <- sigma[j,i] <- rho2^abs(i-j)*sqrt(sigma[i,i]*sigma[j,j])
      }
    }
  }
  
  return(sigma)
}


simX <- function(n, p, mu, rho, miss.prob=0.2, phi=0, missing=c("MCAR", "MAR", "MNAR"), mu.typ = c("const", "var"), s.typ=c("const", "mix", "decay"), pos=TRUE){
  
  mu.typ <- mu.typ[1]
  s.typ <- s.typ[1]
  missing <- missing[1]
  
  if (mu.typ=="const"){
    muvec <- rep(mu,p)
  } else if (mu.typ=="var" & pos==TRUE){
    muvec <- mu*sample(c(1,2),p, replace=T)
  } else if (mu.typ=="var" & pos==FALSE){
    muvec <- sample(c(0, mu[1], - mu[1]), p, replace=TRUE)
  }
  
  if (s.typ=="const") {
    sigma <- matrix(rho, p, p)
    diag(sigma) <- 1
  }  else if (s.typ=="mix"){
    sigma <- my.simSigma(p, pattern="random", diag=1,pCor=round(p/2), rho1=0,rho2=rho, sd.rho2=0.1)
  } else if (s.typ=="decay") {
    sigma <- my.simSigma(p, pattern="decay", diag=1, rho1=0,rho2=rho, sd.rho2=0)
  }
  
  Xcomp <- rmvnorm(n, muvec, sigma)
  
  if (pos)   Xcomp[Xcomp<0] <- 0
  
  if (missing=="MCAR"){
    prob <- rep(miss.prob, length(Xcomp))
  } else if (missing=="MAR"){  
    prob <- NULL
    for (i in 1:n){
      Xobs <- Xcomp[i, 1:2]  ## to make sure we observe these two proteins for sure
      pi <- exp(-(abs(Xobs[1])+abs(Xobs[2])) )  ## decide if the sample is an easy or hard one to observe
      probi <- rep(pi, p) ## set the probability of missing based on the above decision
      probi[1:2] <- 0  ## the two carefully observed proteins have zero probability of missing
      
      prob <- rbind(prob, probi)
      
    }
    prob <- as.matrix(prob)/mean(prob)*miss.prob
    
  } else if (missing=="MNAR" & pos==FALSE){               
    prob <- exp(-phi*(Xcomp)^2)
    theta <- 1#miss.prob/mean(prob)
    prob <- prob*theta
  } else if (missing=="MNAR" & pos==TRUE){
    prob <- exp(-phi*Xcomp) #### changed by pei
    prob <- prob/exp(-phi*min(mu)*0.5)    
  }
  prob[prob>1] <- 1
  prob[prob<0] <- 0
  prob[Xcomp==0] <- 1
  
  X <- Xcomp
  X[which(rbinom(n*p,1,prob)==1)] <- NA
  
  mean(is.na(X))
  range(prob)
  summary(prob)
  summary(as.vector(prob))
  
  
  return(list(X=X, Xcomp=Xcomp, mu=muvec, sigma=sigma, phi=phi, pos=pos))
}


my.solve <- function(X){
  if (!is.matrix(X))  X <- matrix(X, nrow=sqrt(length(X)))
  ss <- svd(X)
  Xinv <- ss$u%*%diag(1/ss$d, nrow=nrow(X), ncol=nrow(X))%*%t(ss$v)
  return(Xinv)
}


get.llik <- function(X, mu, S, phi, phi0){
  
  
  if (length(which(X<0))==0) pos<- TRUE else pos<- FALSE
  
  llik <- 0
  for (j in 1:nrow(X)){
    Xi <- X[j,]
    idxi <- which(!is.na(Xi))
    Xi <- Xi[idxi]
    Si <- as.matrix(S[idxi,idxi])
    Sii <- my.solve(Si)
    oo <-    - log(det(Si))  -(Xi-mu[idxi])%*%Sii%*%(Xi-mu[idxi]) 
    oo <- 0.5*oo
    nmis <- length(X[j,])-length(idxi)
    if (phi==0){
      vv <- 0
    } else {
      if (pos){
        ### new code //Pei
        vvoo <- 1-exp(-phi0-phi*Xi)
        vvoo[vvoo<=0] <- 0.001
        vv1=sum(log(vvoo))
        vv2 <- -phi0*nmis
        if (nmis==0){
          vv3 <- 0
        } else {
          Smi <- as.matrix(S[-idxi,-idxi])
          mu.mis <- mu[-idxi]+matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%(Xi-mu[idxi])
          S.mis <- Smi - matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%matrix(S[idxi,-idxi],ncol=nmis)
          vv3<- -sum(mu.mis)*phi+0.5*sum(S.mis)*phi^2 ### Pei: this formula is the simplified version of the previous vv3. 
        }
        vv <- vv1+vv2+vv3
      } else {
        vv1=sum(log(1-exp(-phi0-phi*Xi^2)))
        vv2 <- -phi0*nmis
        if (nmis==0){
          vv3 <- 0
        } else {
          Smi <- as.matrix(S[-idxi,-idxi])
          mu.mis <- mu[-idxi]+matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%(Xi-mu[idxi])
          S.mis <- Smi - matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%matrix(S[idxi,-idxi],ncol=nmis)
          Smis.inv <- my.solve(S.mis)
          Smii <- Smis.inv
          diag(Smis.inv) <- diag(Smis.inv)+2*phi
          A <- my.solve(Smis.inv)
          
          vv3<- 0.5*(log(det(A)) -log(det(S.mis)) + matrix(mu.mis,nrow=1)%*%(Smii%*%A%*%Smii-Smii)%*%matrix(mu.mis,ncol=1) ) ## haven't tested
        }
        vv <- vv1+vv2+vv3
        
      }
    }
    llik <- llik+ oo+vv
  }  
  
  pllik <- llik
  return(llik) 
}




PEM <- function(X, lambda=NULL, tol=0.001, maxIter=100, K=NULL,delta=1){
  p <- ncol(X)
  N <- nrow(X)
  if (is.null(K)) K= 5
  X.hat <- X
  mu.new <- colMeans(X,na.rm=T)
  Sigma.new <- cov(X, use="pairwise.complete")
  Sigma.new[is.na(Sigma.new)] <- 0
  
  diff <- 999
  
  iter <- 0
  
  if (is.null(lambda)) {
    Lambda <- find.lambda(Sigma.new,N=N,p=p,K=K,delta=delta)
  } else {
    Lambda <- lambda
  }
  Sigma.new <- N/(N+K)*Sigma.new*(N-1)/N + Lambda/(N+K)*diag(1, p, p)
  illik <- 999
  while(iter<maxIter & diff>tol){
    iter <- iter+1
    mu <- mu.new
    Sigma <- Sigma.new
    
    cov.add <- matrix(0, p,p)
    for (i in 1:nrow(X)){
      ii <- which(is.na(X[i,]))
      if (length(ii)>=1){
        Soo <- as.matrix(Sigma[-ii,-ii])
        pi <- nrow(Soo)
        X.hat[i,ii] <-  mu[ii]+Sigma[ii,-ii]%*%my.solve(Soo)%*%(X[i,-ii]-mu[-ii])    ## Key step
        cov.ii <- Sigma[ii,ii] - Sigma[ii,-ii]%*%my.solve(Soo)%*%Sigma[-ii,ii]
        cov.add[ii,ii] <- cov.add[ii,ii]+cov.ii
      }
    }
    
    mu.new <- colMeans(X.hat)
    Sig <- cov(X.hat)*(N-1)/N+cov.add/N
    if (is.null(lambda)) {
      Lambda <- find.lambda(Sig,N=N,p=p,K=K,delta=delta)
    } else {
      Lambda <- lambda
    }
    Sigma.new <- N/(N+K)*(Sig)  + Lambda/(N+K)*diag(1, p, p)
    #diff <- sum(abs(((mu.new-mu)/mu)[which(mu!=0)]))+ sum(abs(((Sigma.new-Sigma)/Sigma)[which(Sigma!=0)])) 
    #if (is.na(diff)) diff <-0  ## if the algorithm does not converge, we let it stop by setting diff=0   
    illik <- c(illik, get.llik(X, mu.new, Sigma.new,phi=0) - sum(diag(Lambda*my.solve(Sigma.new))) -K*log(det(Sigma.new)) ) 
    diff <- abs(illik[iter+1]-illik[iter])/abs(illik[iter])
    if (is.na(diff)) diff <- 0  ## if the algorithm does not converge, we let it stop by setting diff=0
    
  }
  
  return(list(mu=mu, Sigma=Sigma, illik=illik))
}



gets1 <- function(sigma, phi){
  p <- nrow(sigma)
  s1 <- my.solve(sigma)+diag(2*phi, p, p)
  return(s1)
}


get.bg <- function(sigma, mu, phi){
  p <- length(mu)
  s1 <- gets1(sigma, phi)
  s1inv <- my.solve(s1)
  
  ccen <- my.solve(sigma)
  A <- s1inv%*%ccen
  beta <- A%*%mu
  
  gamma <- s1inv 
  
  return(list(beta=beta, gamma=gamma))
}


find.lambda <- function(Sigma,N,p,K, delta){
  ffL <- function(lambda, Sigma, N, p, K){
    Sigma.new <- N/(N+K)*Sigma*(N-1)/N + lambda/(N+K)*diag(1, p, p)
    return(abs(min(as.double(eigen(N*Sigma.new)$value))))
  }
  
  Sigma2 = Sigma
  while (is.complex(eigen(N*Sigma2)$value)){
    delta = delta+1
    Sigma2 = N/(N+K)*Sigma*(N-1)/N + delta/(N+K)*diag(1, p, p)
  }
  Sigma=Sigma2
  
  oo <- -min(as.double(eigen(N*Sigma)$value))
  if (oo>0){
    lambda <- optimize(ffL, lower=0,upper=oo, Sigma=Sigma, N=N, p=p, K=K)$minimum+delta
  } else {
    lambda <- delta
  }
  return(lambda)
}



mnarPEM <- function(X, lambda=NULL, tol=0.001, maxIter=100, phi=1, K=NULL, pos=TRUE, delta=1,phi0){
  p <- ncol(X)
  N <- nrow(X)
  if (is.null(K)){
    K= 5
  } 
  X.hat <- X
  mu.new <- matrix(colMeans(X,na.rm=T),ncol=1)
  Sigma.new <- cov(X, use="pairwise.complete")
  Sigma.new[is.na(Sigma.new)] <- 0
  
  diff <- 999
  
  iter <- 0
  
  if (is.null(lambda)) {
    Lambda <- find.lambda(Sigma.new,N=N,p=p,K=K,delta=delta)
  } else {
    Lambda <- lambda
  }
  Sigma.new <- N/(N+K)*Sigma.new*(N-1)/N + Lambda/(N+K)*diag(1, p, p)
  illik <- 999
  while(iter<maxIter & diff>tol){
    iter <- iter+1
    mu <- mu.new
    Sigma <- Sigma.new
    
    cov.add <- matrix(0, p, p)
    for (i in 1:nrow(X)){
      ii <- which(is.na(X[i,]))
      if (length(ii)>=1){
        Soo <- as.matrix(Sigma[-ii,-ii])
        pi <- nrow(Soo)
        mu.mis <-  mu[ii]+Sigma[ii,-ii]%*%my.solve(Soo)%*%(X[i,-ii]-mu[-ii]) 
        mu.mis <- matrix(mu.mis,ncol=1)
        cov.mis <- Sigma[ii,ii] - Sigma[ii,-ii]%*%my.solve(Soo)%*%Sigma[-ii, ii]
        
        if (phi!=0 & pos==TRUE){
          
          X.hat[i,ii]<- mu.mis - phi*cov.mis%*%matrix(1,nrow=length(mu.mis),ncol=1)
          cov.ii <- cov.mis
          
        } else if (phi!=0 & pos==FALSE){
          
          oo <- get.bg(cov.mis, mu.mis, phi)             
          X.hat[i,ii] <- oo$beta
          cov.ii <- oo$gamma
          
        } else if (phi==0) {
          
          X.hat[i,ii] <- mu.mis
          cov.ii <- cov.mis 
          
        }     
        cov.add[ii,ii] <- cov.add[ii,ii] + cov.ii
      }
    }          
    
    mu.new <- colMeans(X.hat)
    Sig <- cov(X.hat)*(N-1)/N+cov.add/N
    if (is.null(lambda)) {
      Lambda <- find.lambda(Sig,N=N,p=p,K=K,delta=delta)
    } else {
      Lambda <- lambda
    }
    Sigma.new <- N/(N+K)*(Sig)  + Lambda/(N+K)*diag(1, p, p)
    
    illik <- c(illik, get.llik(X, mu.new, Sigma.new,phi=phi,phi0=phi0) - sum(diag(Lambda*my.solve(Sigma.new))) -K*log(det(Sigma.new)) ) 
    diff <- abs(illik[iter+1]-illik[iter])/abs(illik[iter])
    if (is.na(diff)|diff>10^5) {
      return(list(mu=rep(NA,p),Sigma=diag(1,p),illik=NA))
    }
  }
  
  return(list(mu=mu, Sigma=Sigma,illik=illik[length(illik)]))
}


comp.imp <- function(simDat){
  delta <- 5
  X <- simDat$X
  Xcomp <- simDat$Xcomp
  mu <- simDat$mu
  sigma <- simDat$sigma
  phi <- simDat$phi
  pos <- simDat$pos
  
  idx <- which(rowSums(!is.na(X))<=1)
  if (length(idx)>0){
    X <- X[-idx,]
    Xcomp <- Xcomp[-idx,]
  }
  
  p <- ncol(X)
  N <- nrow(X)
  k.max=min(P,N)-1
  
  ## complete estimate
  mu.comp <- colMeans(Xcomp)
  sigma.comp <- cov(Xcomp)
  
  k1=cv.kNNImpute(X, k.max=k.max)$k
  X.imp <-  impute.knn(X, k=k1, rowmax=1, colmax=1)$data
  
  
  if (phi==0) phi0T <- -log(mean(is.na(X))) else { phi0T <- 0}
  res.mnp <- mnarPEM(X, phi=phi, pos=pos, delta=delta,phi0=phi0T)
  mu.mnp <- res.mnp$mu
  sigma.mnp <- res.mnp$Sigma
  
  X.hat <- X
  mu <- mu.mnp
  Sigma <- sigma.mnp
  for (i in 1:nrow(X)){
    ii <- which(is.na(X[i,]))
    if (length(ii)>=1){
      Soo <- as.matrix(Sigma[-ii,-ii])
      pi <- nrow(Soo)
      mu.mis <-  mu[ii]+Sigma[ii,-ii]%*%my.solve(Soo)%*%(X[i,-ii]-mu[-ii]) 
      mu.mis <- matrix(mu.mis,ncol=1)
      cov.mis <- Sigma[ii,ii] - Sigma[ii,-ii]%*%my.solve(Soo)%*%Sigma[-ii, ii]
      
      if (phi!=0 & pos==TRUE){
        X.hat[i,ii]<- mu.mis - phi*cov.mis%*%matrix(1,nrow=length(mu.mis),ncol=1)
        cov.ii <- cov.mis   
      } else if (phi!=0 & pos==FALSE){
        oo <- get.bg(cov.mis, mu.mis, phi)             
        X.hat[i,ii] <- oo$beta
        cov.ii <- oo$gamma
      } 
    } }
  
  RMSE.knn <- sqrt(mean((Xcomp-X.imp)^2))/(max(Xcomp)-min(Xcomp))
  RMSE.mnp <- sqrt(mean((Xcomp-X.hat)^2))/(max(Xcomp)-min(Xcomp))
  return(c(RMSE.knn, RMSE.mnp))
}




comp.stat <- function(simDat,delta=1){
  
  X <- simDat$X
  Xcomp <- simDat$Xcomp
  mu <- simDat$mu
  sigma <- simDat$sigma
  phi <- simDat$phi
  pos <- simDat$pos
  
  idx <- which(rowSums(!is.na(X))<=1)
  if (length(idx)>0){
    X <- X[-idx,]
    Xcomp <- Xcomp[-idx,]
  }
  
  p <- ncol(X)
  N <- nrow(X)
  
  ## complete estimate
  mu.comp <- colMeans(Xcomp)
  sigma.comp <- cov(Xcomp)
  
  ## available estimate
  mu.av <- colMeans(X,na.rm=T)
  sigma.av <- cov(X, use="pairwise.complete")
  sigma.av[is.na(sigma.av)] <- 0
  
  ## available + penalization
  mu.avp <- mu.av
  sigma.avp <- sigma.av+diag(delta/N, p, p)
  
  ## KNNimputation + complete estimate
  k.max=min(p,N)-1
  k1=cv.kNNImpute(X,k.max=k.max)$k
  X.imp <-  impute.knn(X, k=k1, rowmax=1, colmax=1)$data
  
  mu.imp <- colMeans(X.imp)
  sigma.imp <- cov(X.imp)+diag(delta/N, p, p)
  
  
  ## EM
  if (p<N){
    res.em <- PEM(X, lambda=0, K=0)
    mu.em <- res.em$mu
    sigma.em <- res.em$Sigma
  } else {
    mu.em <- rep(NA,p)
    sigma.em <- matrix(0, p,p)
  }
  
  ## PEM, K=0
  #res.pemK0 <- PEM(X,K=0, delta=delta)
  #mu.pemK0 <- res.pemK0$mu
  #sigma.pemK0 <- res.pemK0$Sigma
  
  ## PEM
  res.pem <- PEM(X, delta=delta)
  mu.pem <- res.pem$mu
  sigma.pem <- res.pem$Sigma
  
  
  ## MnarPEM
  if (phi==0) phi0T <- -log(mean(is.na(X))) else { phi0T <- 0}
  res.mnp <- mnarPEM(X, phi=phi, pos=pos, delta=delta,phi0=phi0T)
  mu.mnp <- res.mnp$mu
  sigma.mnp <- res.mnp$Sigma
  
  
  ## get an AV-case-based estimate of phi first.
  
  if (sum(is.na(mu.mnp))>0) return(rep(NA,12)) else{
    
    ffc <- function(para, prob, x, pos){
      phi <- para[1]
      const <- para[2]
      if (pos){
        res <- sum((prob- const*exp(-phi*x))^2)
      } else {
        res <- sum((prob- const*exp(-phi*x^2))^2)
      }
      return(res)
    }
    b <- colMeans(is.na(X),na.rm=TRUE)
    a <- colMeans(X,na.rm=TRUE)
    oo <- optim(c(0.1,1), ffc,prob=b,x=a,pos=pos,lower=0,upper=1+as.integer(!pos), method="L-BFGS-B")
    phi.hat <- oo$par[1]
    phi.hat[phi.hat<0] <- 0
    phi0 <- -log(oo$par[2])
    
    
    if (pos){
      phi.vec <- phi.hat+(-5):5/100
      phi.vec <- phi.vec[phi.vec>=0]
    } else{
      if (phi.hat>0.1) {
        phi0 <- 0
        phi.vec <- (1:12)/10 
      } else {
        phi0 <- -log(mean(is.na(X)))
        phi.vec <- (0:5)/10
      }
      
    }
    
    illik <- NULL 
    for (phi in phi.vec){  
      illik <- rbind(illik, mnarPEM(X=X,phi=phi,pos=pos,delta=delta,phi0=phi0)$illik)  ## here you can change the plikelihood function
    }
    phi.hat <- phi.vec[which.max(illik)]
    print(phi.hat)
    
    res.mnp2 <- mnarPEM(X, phi=phi.hat, pos=pos, delta=delta,phi0=phi0)
    mu.mnp2 <- res.mnp2$mu
    sigma.mnp2 <- res.mnp2$Sigma
    
    
    mseC <- function(mm, ss){
      temp1 <- mean((mm-mu)^2)
      temp2 <- mean((ss-sigma)^2)
      #temp3 <- mean((my.solve(ss)-solve(sigma))^2)
      #temp4 <- is.positive.definite(ss)
      return(c(temp1, temp2) )
    }
    
    ##mse.c <- c(mseC(mu.comp, sigma.comp), mseC(mu.av, sigma.av), mseC(mu.avp, sigma.avp), mseC(mu.imp, sigma.imp), mseC(mu.em, sigma.em), 
    ##                       mseC(mu.pemK0, sigma.pemK0), mseC(mu.pem, sigma.pem), mseC(mu.mnp, sigma.mnp),  mseC(mu.mnp2, sigma.mnp2)   )
    
    mse.c <- c(mseC(mu.comp, sigma.comp), mseC(mu.avp, sigma.avp),  mseC(mu.imp, sigma.imp), mseC(mu.em, sigma.em), 
               mseC(mu.pem, sigma.pem), mseC(mu.mnp, sigma.mnp),  mseC(mu.mnp2, sigma.mnp2)  )
    
    
    
    res <- matrix(mse.c, nrow=2)
    
    rownames(res) <- c("RMSE.mu","RMSE.s")
    #colnames(res) <- c("COMP", "AV", "AVP", "IMP", "EM", "PEM.K0","PEM", "MnarPEM")
    
    return(res[,-1]/res[,1])
  }
}





par.result <- function(cl, p, n, mu, rho, miss.prob=0.2, phi=1, missing=c("MCAR", "MAR", "MNAR"), 
                       mu.typ=c("const", "var"), s.typ=c("const", "mix", "decay"), nsim=100, pos=TRUE, delta=delta){
  library(corpcor)
  library(mvtnorm)
  library(imputation)
  library(impute)
  
  ff <- function(i, p, n, mu, rho, miss.prob, phi, missing, mu.typ, s.typ, pos, delta){
    set.seed(i)
    simDat <- simX(p=p,n=n,mu=mu,rho=rho,miss.prob=miss.prob, phi=phi, missing=missing, mu.typ=mu.typ, s.typ=s.typ, pos=pos)
    idx <- which(colSums(!is.na(simDat$X))<=2)
    if (length(idx)>0) {
      for (j in idx){
        dd <- simDat$Xcomp[,j]
        zz <- order(-abs(dd))[1:2]
        simDat$X[zz,j] <- dd[zz]  
      }}
    
    mis <- mean(is.na(simDat$X))
    oo <- t(comp.stat(simDat, delta=delta))
    res <- c(matrix(oo, byrow=TRUE, nrow=1),mis)
    return(res)
  }
  
  res <- parLapply(cl, 1:nsim, ff, p=p, n=n, mu=mu,rho=rho,miss.prob=miss.prob, phi=phi, missing=missing[1], mu.typ=mu.typ[1], s.typ=s.typ[1], pos=pos,delta=delta)
  res <- matrix(unlist(res), byrow=TRUE, nrow=nsim)
  
  print(mean(res[,ncol(res)]))
  res <- res[,-ncol(res)]
  
  
  
  return(res)
}



par.compimp <- function(cl, p, n, mu, rho, miss.prob=0.2, phi=1, missing=c("MCAR", "MAR", "MNAR"), 
                        mu.typ=c("const", "var"), s.typ=c("const", "mix", "decay"), nsim=100, pos=TRUE, delta=delta){
  library(corpcor)
  library(mvtnorm)
  library(imputation)
  
  ff <- function(i, p, n, mu, rho, miss.prob, phi, missing, mu.typ, s.typ, pos, delta){
    set.seed(i)
    simDat <- simX(p=p,n=n,mu=mu,rho=rho,miss.prob=miss.prob, phi=phi, missing=missing, mu.typ=mu.typ, s.typ=s.typ, pos=pos)
    idx <- which(colSums(!is.na(simDat$X))<=2)
    if (length(idx)>0) {
      for (j in idx){
        dd <- simDat$Xcomp[,j]
        zz <- order(-abs(dd))[1:2]
        simDat$X[zz,j] <- dd[zz]  
      }}
    
    mis <- mean(is.na(simDat$X))
    oo <- comp.imp(simDat)
    return(oo)
  }
  
  res <- parLapply(cl, 1:nsim, ff, p=p, n=n, mu=mu,rho=rho,miss.prob=miss.prob, phi=phi, missing=missing[1], mu.typ=mu.typ[1], s.typ=s.typ[1], pos=pos,delta=delta)
  res <- matrix(unlist(res), byrow=TRUE, nrow=nsim)  
  
  return(res)
}

