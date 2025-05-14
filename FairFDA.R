
fourier.basis <- function(s, tseq) {
  nt <- length(tseq)
  basis.mat <- matrix(NA, nrow=nt, ncol=s)
  for (j in 1:s) {
    basis.mat[,j] <- sqrt(2)*cos(j*pi*tseq)
  }
  return(basis.mat)
}

classGen <- function(n, pa, p_ycona){
  class_a = rbinom(n, 1, pa) # generate A
  na_1 = sum(class_a==1); na_0 = sum(class_a==0)
  class_0y = rbinom(na_0, 1, p_ycona[1]) # generate Y given A=0
  class_1y = rbinom(na_1, 1, p_ycona[2]) # generate Y given A=1
  n_01 = sum(class_0y==1); n_00 = na_0 - n_01
  n_11 = sum(class_1y==1); n_10 = na_1 - n_11
  return(c(n_00, n_01, n_10, n_11))
}

dataGen <- function(n, pa, p_ycona, tseq, alpha, beta, type='gauss', const_snr0=1, const_snr1=1, s=50){
  
  class_n = classGen(n, pa, p_ycona)
  n_00 = class_n[1]; n_01 = class_n[2]
  n_10 = class_n[3]; n_11 = class_n[4]
  na_0 = n_00+n_01
  na_1 = n_10+n_11
  
  basis.mat <- fourier.basis(s,tseq)
  ### data under A=0
  if(type=='gauss'){
    z <- matrix(rnorm(na_0*s, 0, 1), na_0, s) # gaussian
  }
  if(type=='unif'){
    z <- matrix(runif(na_0*s, -sqrt(3), sqrt(3)), na_0, s) # uniform
  }
  D0 <- diag((1:s)^(-alpha/2))
  score0 <- z%*%D0 
  x0 <- tcrossprod(score0, basis.mat)
  x_00 = x0[1:n_00,,drop=FALSE]
  mu_01 = matrix((-1)^(1:s)*(1:s)^(-beta)*const_snr0, n_01, s, byrow = TRUE)
  x_01 = x0[(n_00+1):na_0,,drop=FALSE] + tcrossprod(mu_01, basis.mat)
  ### data under A=1
  if(type=='gauss'){
    z <- matrix(rnorm(na_1*s, 0, 1), na_1, s) # gaussian
  }
  if(type=='unif'){
    z <- matrix(runif(na_1*s, -sqrt(3), sqrt(3)), na_1, s) # uniform
  }
  D1 <- sqrt(2) * diag((1:s)^(-alpha/2)) ### different eigenvalues
  score1 <- z%*%D1
  x1 <- tcrossprod(score1, basis.mat)
  x_10 = x1[1:n_10,,drop=FALSE]
  mu_11 = matrix(sqrt(2)*(-1)^(1:s)*(1:s)^(-beta)*const_snr1, n_11, s, byrow = TRUE)
  x_11 = x1[(n_10+1):na_1,,drop=FALSE] + tcrossprod(mu_11, basis.mat)

  return(list(x_00=x_00, x_01=x_01, x_10=x_10, x_11=x_11))
  
}

sampleSplit <- function(x_00, x_01, x_10, x_11, split_ratio=0.5){
  
  n_00 = nrow(x_00); n_01 = nrow(x_01)
  n_10 = nrow(x_10); n_11 = nrow(x_11)
  
  idx = sample(1:n_00, ceiling(n_00*split_ratio), replace = FALSE)
  x_00_train = x_00[idx,,drop=FALSE]
  x_00_cal = x_00[-idx,,drop=FALSE]
  
  idx = sample(1:n_01, ceiling(n_01*split_ratio), replace = FALSE)
  x_01_train = x_01[idx,,drop=FALSE]
  x_01_cal = x_01[-idx,,drop=FALSE]
  
  idx = sample(1:n_10, ceiling(n_10*split_ratio), replace = FALSE)
  x_10_train = x_10[idx,,drop=FALSE]
  x_10_cal = x_10[-idx,,drop=FALSE]
  
  idx = sample(1:n_11, ceiling(n_11*split_ratio), replace = FALSE)
  x_11_train = x_11[idx,,drop=FALSE]
  x_11_cal = x_11[-idx,,drop=FALSE]
  
  return(list(x_00_train=x_00_train, x_00_cal=x_00_cal,
              x_01_train=x_01_train, x_01_cal=x_01_cal,
              x_10_train=x_10_train, x_10_cal=x_10_cal,
              x_11_train=x_11_train, x_11_cal=x_11_cal))
}


getEigenhat <- function(covfun, J){
  
  nt = nrow(covfun)
  t_by = 1 / (nt-1)
  eigenhat <- eigen(covfun)
  evalhat <- eigenhat$values*t_by
  evalhat <- evalhat[evalhat>0]
  if(length(evalhat)>=J){
    n_pc <- J
  }else{
    n_pc <- length(evalhat)
  } 
  efunhat <- eigenhat$vectors[,1:n_pc,drop=F]/sqrt(t_by)
  evalhat <- evalhat[1:n_pc]
  
  return(list(efunhat=efunhat, evalhat=evalhat))
  
}

getEst <- function(x_00, x_01, x_10, x_11, J){
  
  n_00 = nrow(x_00); n_01 = nrow(x_01)
  n_10 = nrow(x_10); n_11 = nrow(x_11)
  na_0 = n_00+n_01
  na_1 = n_10+n_11; n = na_0+na_1
  nt = ncol(x_00)
  t_by = 1/(nt-1) ### common design
  phat_00 = n_00/n; phat_01 = n_01/n
  phat_10 = n_10/n; phat_11 = n_11/n
  
  ### estimate under A=0
  muhat_00 = colMeans(x_00)
  muhat_01 = colMeans(x_01)
  covhat_0 = cov(x_00)*n_00/na_0 + cov(x_01)*n_01/na_0
  eigenhat = getEigenhat(covhat_0, J)
  evalhat_0 = eigenhat$evalhat
  efunhat_0 = eigenhat$efunhat
  thetahat_00 = crossprod(efunhat_0, muhat_00)*t_by
  thetahat_01 = crossprod(efunhat_0, muhat_01)*t_by
  
  ### estimate under A=1
  muhat_10 = colMeans(x_10)
  muhat_11 = colMeans(x_11)
  covhat_1 = cov(x_10)*n_10/na_1 + cov(x_11)*n_11/na_1
  eigenhat = getEigenhat(covhat_1, J)
  evalhat_1 = eigenhat$evalhat
  efunhat_1 = eigenhat$efunhat
  thetahat_10 = crossprod(efunhat_1, muhat_10)*t_by
  thetahat_11 = crossprod(efunhat_1, muhat_11)*t_by
  
  return(list(efunhat_0=efunhat_0, efunhat_1=efunhat_1,
              evalhat_0=evalhat_0, evalhat_1=evalhat_1,
              thetahat_00=thetahat_00, thetahat_01=thetahat_01,
              thetahat_10=thetahat_10, thetahat_11=thetahat_11,
              phat_00=phat_00, phat_01=phat_01,
              phat_10=phat_10, phat_11=phat_11))
  
}

### default parameters: s0=-1, s1=1, b0=b1=0 (DO)
### tau=0, without fairness constraint
### group-aware prediction
getPrehat <- function(x_00, x_01, x_10, x_11, est_list, tau=0, s0=-1, s1=1, b0=0, b1=0){
  
  n_00 = nrow(x_00); n_01 = nrow(x_01)
  n_10 = nrow(x_10); n_11 = nrow(x_11)
  nt = ncol(x_00)
  t_by = 1 / (nt-1)
  phat_00 = est_list$phat_00; phat_01 = est_list$phat_01
  phat_10 = est_list$phat_10; phat_11 = est_list$phat_11
  
  ### classify under A=0
  zetahat_00 = x_00 %*% est_list$efunhat_0 * t_by
  zetahat_01 = x_01 %*% est_list$efunhat_0 * t_by
  J = ncol(est_list$efunhat_0)
  tmp1 = (est_list$thetahat_01 - est_list$thetahat_00) / est_list$evalhat_0
  tmp2 = sum((est_list$thetahat_01 - est_list$thetahat_00)^2 / est_list$evalhat_0) / 2
  etahat_00 = exp((zetahat_00 - matrix(est_list$thetahat_00, n_00, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  etahat_01 = exp((zetahat_01 - matrix(est_list$thetahat_00, n_01, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  pre_00 =  ((phat_01-tau*s0)*etahat_00) > (phat_00+tau*b0)
  pre_01 =  ((phat_01-tau*s0)*etahat_01) > (phat_00+tau*b0)
  
  ### classify under A=1
  zetahat_10 = x_10 %*% est_list$efunhat_1 * t_by
  zetahat_11 = x_11 %*% est_list$efunhat_1 * t_by
  J = ncol(est_list$efunhat_1)
  tmp1 = (est_list$thetahat_11 - est_list$thetahat_10) / est_list$evalhat_1
  tmp2 = sum((est_list$thetahat_11 - est_list$thetahat_10)^2 / est_list$evalhat_1) / 2
  etahat_10 = exp((zetahat_10 - matrix(est_list$thetahat_10, n_10, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  etahat_11 = exp((zetahat_11 - matrix(est_list$thetahat_10, n_11, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  pre_10 =  ((phat_11-tau*s1)*etahat_10) > (phat_10+tau*b1)
  pre_11 =  ((phat_11-tau*s1)*etahat_11) > (phat_10+tau*b1)
  
  return(list(pre_00=pre_00, pre_01=pre_01,
              pre_10=pre_10, pre_11=pre_11))
  
}

getDhat <- function(x_00, x_01, x_10, x_11, est_list, tau, s0=-1, s1=1, b0=0, b1=0){
  
  n_00 = nrow(x_00); n_01 = nrow(x_01)
  n_10 = nrow(x_10); n_11 = nrow(x_11)
  nt = ncol(x_00)
  t_by = 1 / (nt-1)
  phat_00 = est_list$phat_00; phat_01 = est_list$phat_01
  phat_10 = est_list$phat_10; phat_11 = est_list$phat_11
  
  zetahat_10 = x_10 %*% est_list$efunhat_1 * t_by
  zetahat_11 = x_11 %*% est_list$efunhat_1 * t_by
  J = ncol(est_list$efunhat_1)
  tmp1 = (est_list$thetahat_11 - est_list$thetahat_10) / est_list$evalhat_1
  tmp2 = sum((est_list$thetahat_11 - est_list$thetahat_10)^2 / est_list$evalhat_1) / 2
  etahat_10 = exp((zetahat_10 - matrix(est_list$thetahat_10, n_10, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  etahat_11 = exp((zetahat_11 - matrix(est_list$thetahat_10, n_11, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  prehat_10 = ((phat_11-tau*s1)*etahat_10) > (phat_10+tau*b1)
  prehat_11 = ((phat_11-tau*s1)*etahat_11) > (phat_10+tau*b1)
  
  zetahat_00 = x_00 %*% est_list$efunhat_0 * t_by
  zetahat_01 = x_01 %*% est_list$efunhat_0 * t_by
  J = ncol(est_list$efunhat_0)
  tmp1 = (est_list$thetahat_01 - est_list$thetahat_00) / est_list$evalhat_0
  tmp2 = sum((est_list$thetahat_01 - est_list$thetahat_00)^2 / est_list$evalhat_0) / 2
  etahat_00 = exp((zetahat_00 - matrix(est_list$thetahat_00, n_00, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  etahat_01 = exp((zetahat_01 - matrix(est_list$thetahat_00, n_01, J, byrow=TRUE)) %*%  tmp1 - tmp2)
  prehat_00 = ((phat_01-tau*s0)*etahat_00) > (phat_00+tau*b0)
  prehat_01 = ((phat_01-tau*s0)*etahat_01) > (phat_00+tau*b0)
  
  Dhat = s0*mean(prehat_01) + b0*mean(prehat_00) + s1*mean(prehat_11) + b1*mean(prehat_10)
  
  return(Dhat)
  
}

getTauhat <- function(x_00, x_01, x_10, x_11, est_list, delta, eps=0, s0=-1, s1=1, b0=0, b1=0){
  
  # nt = ncol(x_01)
  delta = delta-eps
  phat_00 = est_list$phat_00
  phat_01 = est_list$phat_01
  phat_10 = est_list$phat_10
  phat_11 = est_list$phat_11
  
  tmp = getTauRange(phat_00, phat_01, phat_10, phat_11, s0=s0, s1=s1, b0=b0, b1=b1)
  tau_min = tmp[1]
  tau_max = tmp[2]
  
  Dhat0 = getDhat(x_00, x_01, x_10, x_11, est_list, tau=0, s0=s0, s1=s1, b0=b0, b1=b1)
  
  
  if(abs(Dhat0)<=delta){
    tauhat = 0
  } 
  if(Dhat0>delta){ ### tauhat>0
    tau_min = 0
    tau_old = tau_min
    tau_new = (tau_min+tau_max)/2
    while((tau_max-tau_min)>1e-6){
      Dhat_new = getDhat(x_00, x_01, x_10, x_11, est_list, tau=tau_new, s0=s0, s1=s1, b0=b0, b1=b1)
      tau_old = tau_new
      if(Dhat_new>delta){
        tau_min = tau_old
      }else{
        tau_max = tau_old
      }
      tau_new = (tau_min+tau_max)/2
    }
    tauhat = tau_old
  } 
  if(Dhat0<(-delta)){ ### tauhat<0
    tau_max = 0
    tau_old = tau_max
    tau_new = (tau_min+tau_max)/2
    while((tau_max-tau_min)>1e-6){
      Dhat_new = getDhat(x_00, x_01, x_10, x_11, est_list, tau=tau_new, s0=s0, s1=s1, b0=b0, b1=b1)
      tau_old = tau_new
      if(Dhat_new>(-delta)){
        tau_min = tau_old
      }else{
        tau_max = tau_old
      }
      tau_new = (tau_min+tau_max)/2
    }
    tauhat = tau_old
  }
  
  return(tauhat)
  
}

cvTune <- function(x_00, x_01, x_10, x_11, Jseq, nfold=5){
  nJ = length(Jseq)
  n_00 = nrow(x_00); n_01 = nrow(x_01)
  n_10 = nrow(x_10); n_11 = nrow(x_11)
  n_00_val = floor(n_00/nfold)
  n_01_val = floor(n_01/nfold)
  n_10_val = floor(n_10/nfold)
  n_11_val = floor(n_11/nfold)
  n_val = n_00_val + n_01_val + n_10_val + n_11_val
  ###
  error <- matrix(NA, nJ, nfold)
  for(i in 1:nfold){
    idx_00 = ((i-1)*n_00_val+1):(i*n_00_val)
    x_00_val = x_00[idx_00,,drop=FALSE]
    x_00_train = x_00[-idx_00,,drop=FALSE]
    idx_01 = ((i-1)*n_01_val+1):(i*n_01_val)
    x_01_val = x_01[idx_01,,drop=FALSE]
    x_01_train = x_01[-idx_01,,drop=FALSE]
    idx_10 = ((i-1)*n_10_val+1):(i*n_10_val)
    x_10_val = x_10[idx_10,,drop=FALSE]
    x_10_train = x_10[-idx_10,,drop=FALSE]
    idx_11 = ((i-1)*n_11_val+1):(i*n_11_val)
    x_11_val = x_11[idx_11,,drop=FALSE]
    x_11_train = x_11[-idx_11,,drop=FALSE]
    for(j in 1:nJ){
      est_list <- getEst(x_00_train, x_01_train, x_10_train, x_11_train, J=Jseq[j])
      pre_list <- getPrehat(x_00_val, x_01_val, x_10_val, x_11_val, est_list, tau=0)
      error[j, i] = (sum(pre_list$pre_00)+sum(pre_list$pre_10)+sum(1-pre_list$pre_01)+sum(1-pre_list$pre_11))/n_val
    }
  }
  Jhat = Jseq[which.min(rowMeans(error))]
  
  return(Jhat)
  
}


TuneCal <- function(x_00, x_01, x_10, x_11, J, delta, eps_seq, s0=-1, s1=1, b0=0, b1=0, rho=0.05){

  nseq = length(eps_seq)
  Dhat = matrix(NA, 100, nseq)
  for(j in 1:nseq){
    eps = eps_seq[j]
    for(i in 1:100){
      xx = sampleSplit(x_00, x_01, x_10, x_11, split_ratio=0.5)
      x_00_tmp = xx$x_00_train; x_00_test = xx$x_00_cal
      x_01_tmp = xx$x_01_train; x_01_test = xx$x_01_cal
      x_10_tmp = xx$x_10_train; x_10_test = xx$x_10_cal
      x_11_tmp = xx$x_11_train; x_11_test = xx$x_11_cal

      xx = sampleSplit(x_00_tmp, x_01_tmp, x_10_tmp, x_11_tmp, split_ratio=0.5)
      x_00_train = xx$x_00_train; x_00_cal = xx$x_00_cal
      x_01_train = xx$x_01_train; x_01_cal = xx$x_01_cal
      x_10_train = xx$x_10_train; x_10_cal = xx$x_10_cal
      x_11_train = xx$x_11_train; x_11_cal = xx$x_11_cal

      est_list = getEst(x_00_train, x_01_train, x_10_train, x_11_train, J=J)
      tauhat = getTauhat(x_00_cal, x_01_cal, x_10_cal, x_11_cal, est_list, delta, eps=eps, s0=s0, s1=s1, b0=b0, b1=b1)
      pre_list1 <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list, tau=tauhat, s0=s0, s1=s1, b0=b0, b1=b1)
      # Dhat1 = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list, tau=tauhat, s0=s0, s1=s1, b0=b0, b1=b1)

      est_list = getEst(x_00_cal, x_01_cal, x_10_cal, x_11_cal, J=J)
      tauhat = getTauhat(x_00_train, x_01_train, x_10_train, x_11_train, est_list, delta, eps=eps, s0=s0, s1=s1, b0=b0, b1=b1)
      pre_list2 <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list, tau=tauhat, s0=s0, s1=s1, b0=b0, b1=b1)
      # Dhat2 = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list, tau=tauhat, s0=s0, s1=s1, b0=b0, b1=b1)

      pre_00 = rbinom(nrow(x_00_test), 1, (pre_list1$pre_00+pre_list2$pre_00)/2)
      pre_01 = rbinom(nrow(x_01_test), 1, (pre_list1$pre_01+pre_list2$pre_01)/2)
      pre_10 = rbinom(nrow(x_10_test), 1, (pre_list1$pre_10+pre_list2$pre_10)/2)
      pre_11 = rbinom(nrow(x_11_test), 1, (pre_list1$pre_11+pre_list2$pre_11)/2)
      
      Dhat[i, j] <- s0*mean(pre_01) + b0*mean(pre_00) + s1*mean(pre_11) + b1*mean(pre_10)
    }
  }
  D_quan = apply(abs(Dhat), 2, quantile, probs=1-rho)
  if(any(D_quan<=delta)){
    epshat = eps_seq[min(which(D_quan<=delta))]
  }else{
    epshat = eps_seq[which.min(D_quan)]
  }

  return(epshat)

}


### explicit formula for homogeneous Gaussian processes
getoracleRisk <- function(oracle_list, tau=0, s0=-1, s1=1, b0=0, b1=0){
  snr0 = oracle_list$snr0; snr1 = oracle_list$snr1
  p_00 = oracle_list$p_00; p_01 = oracle_list$p_01
  p_10 = oracle_list$p_10; p_11 = oracle_list$p_11
  
  temp00 = -snr0/2 - log((p_00+tau*b0)/(p_01-tau*s0))/snr0
  temp01 = -snr0/2 + log((p_00+tau*b0)/(p_01-tau*s0))/snr0
  temp10 = -snr1/2 - log((p_10+tau*b1)/(p_11-tau*s1))/snr1
  temp11 = -snr1/2 + log((p_10+tau*b1)/(p_11-tau*s1))/snr1
  
  risk = p_00*pnorm(temp00) + p_01*pnorm(temp01) + p_10*pnorm(temp10) + p_11*pnorm(temp11)
  
  return(risk)
  
}

getoracleD <- function(oracle_list, tau=0, s0=-1, s1=1, b0=0, b1=0){
  
  snr0 = oracle_list$snr0; snr1 = oracle_list$snr1
  p_00 = oracle_list$p_00; p_01 = oracle_list$p_01
  p_10 = oracle_list$p_10; p_11 = oracle_list$p_11
  
  ### disparity measure
  if((p_11-tau*s1)>0){
    if((p_10+tau*b1)>0){
      temps1 = snr1/2 - log((p_10+tau*b1)/(p_11-tau*s1))/snr1
      tempb1 = -snr1/2 - log((p_10+tau*b1)/(p_11-tau*s1))/snr1
    }
    else{
      temps1 = Inf
      tempb1 = Inf
    }
  }
  if((p_11-tau*s1)<=0){
    if((p_10+tau*b1)>=0){
      temps1 = -Inf
      tempb1 = -Inf
    }
    else{
      if((p_11-tau*s1)==0){
        temps1 = Inf
        tempb1 = Inf
      }
      if((p_11-tau*s1)<0){
        temps1 = snr1/2 - log((p_10+tau*b1)/(p_11-tau*s1))/snr1
        tempb1 = -snr1/2 - log((p_10+tau*b1)/(p_11-tau*s1))/snr1
      }
    }
  }
  if((p_01-tau*s0)>0){
    if((p_00+tau*b0)>0){
      temps0 = snr0/2 - log((p_00+tau*b0)/(p_01-tau*s0))/snr0
      tempb0 = -snr0/2 - log((p_00+tau*b0)/(p_01-tau*s0))/snr0
    }
    else{
      temps0 = Inf
      tempb0 = Inf
    }
  }
  if((p_01-tau*s0)<=0){
    if((p_00+tau*b0)>=0){
      temps0 = -Inf
      tempb0 = -Inf
    }
    else{
      if((p_01-tau*s0)==0){
        temps0 = Inf
        tempb0 = Inf
      }
      else{
        temps0 = snr0/2 - log((p_00+tau*b0)/(p_01-tau*s0))/snr0
        tempb0 = -snr0/2 - log((p_00+tau*b0)/(p_01-tau*s0))/snr0
      }
    }
  }
  
  D = s1*pnorm(temps1) + s0*pnorm(temps0) + b1*pnorm(tempb1) + b0*pnorm(tempb0)
  
  return(D)
  
}

getTauRange <- function(p_00, p_01, p_10, p_11, s0=-1, s1=1, b0=0, b1=0){
  ### disparity measures: s1>=0, b1>=0, s0<=0, b0<=0
  
  if(b0<0){
    tau_max = min(p_11/s1, -p_00/b0)
  }
  if(b0==0){
    tau_max = p_11/s1
  }
  
  if(s0<0){
    tau_min = max(p_01/s0, -p_10/b1)
  }
  if(s0==0){
    tau_min = -p_10/b1
  }
  return(c(tau_min, tau_max))
}

getTaustar <- function(oracle_list, delta, s0=-1, s1=1, b0=0, b1=0){
  
  p_00 = oracle_list$p_00; p_01 = oracle_list$p_01
  p_10 = oracle_list$p_10; p_11 = oracle_list$p_11
  
  tmp = getTauRange(p_00, p_01, p_10, p_11, s0=s0, s1=s1, b0=b0, b1=b1)
  tau_min = tmp[1]
  tau_max = tmp[2]
  D0 = getoracleD(oracle_list, tau=0, s0=s0, s1=s1, b0=b0, b1=b1)
  
  if(abs(D0)<=delta){
    taustar = 0
  }
  if(D0>delta){ ### taustar>0
    tau_min = 0
    tau_old = tau_min
    tau_new = (tau_min+tau_max)/2
    while((tau_max-tau_min)>1e-6){
      D_new = getoracleD(oracle_list, tau=tau_new, s0=s0, s1=s1, b0=b0, b1=b1)
      tau_old = tau_new
      if(D_new>delta){
        tau_min = tau_old
      }else{
        tau_max = tau_old
      }
      tau_new = (tau_min+tau_max)/2
    }
    taustar = tau_old
  }
  if(D0<(-delta)){ ### taustar<0
    tau_max = 0
    tau_old = tau_max
    tau_new = (tau_min+tau_max)/2
    while((tau_max-tau_min)>1e-6){
      D_new = getoracleD(oracle_list, tau=tau_new, s0=s0, s1=s1, b0=b0, b1=b1)
      tau_old = tau_new
      if(D_new>-delta){
        tau_min = tau_old
      }else{
        tau_max = tau_old
      }
      tau_new = (tau_min+tau_max)/2
    }
    taustar = tau_old
  }
  
  return(taustar)
  
}
