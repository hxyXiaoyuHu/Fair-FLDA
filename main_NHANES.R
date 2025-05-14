
load("NHANES/savedemo_race.dat")
dim(savedemo)

# preprocessing following NHANES protocol
# keep intensity values between 1 and 1000 for each observation 
# remove subjects with no more than 100 observations
load("NHANES/remainingidrange1to1000.dat")
# length(remainingid) # 7014

load("NHANES/allquantilefunctionrange1to1000.dat")
# dim(allquantilefunction) # 7014 * 1001

n <- length(remainingid)
age <- rep(NA, n)
race <- rep(NA, n)
quantilefunction<-matrix(NA, n, ncol(allquantilefunction))
for (i in 1:n){
  tempid <- remainingid[i]
  index <- which(savedemo[,1]==tempid)
  age[i] <- savedemo[index, 2]
  # race: 1-Mexican American, 2-Other Hispanic, 3-Non-Hispanic White, 4-Non-Hispanic Black, 5-Other Non-Hispanic
  race[i] <- savedemo[index, 8]
  index2<-which(remainingid==tempid)
  quantilefunction[i,]<-allquantilefunction[index2,]
}

lgrid<-1001
grid<-seq(0, 1, length.out=lgrid) # quantile grid
idx_grid = which(grid %in% seq(0,1,length.out=101))

# 3252 instances
idx_00 = which(race==3 & age<=20)
idx_01 = which(race==3 & age>=50)
idx_10 = which(race==4 & age<=20)
idx_11 = which(race==4 & age>=50)

xx_00 = quantilefunction[idx_00, idx_grid]
xx_01 = quantilefunction[idx_01, idx_grid]
xx_10 = quantilefunction[idx_10, idx_grid]
xx_11 = quantilefunction[idx_11, idx_grid] 
n_00 = dim(xx_00)[1]
n_01 = dim(xx_01)[1]
n_10 = dim(xx_10)[1]
n_11 = dim(xx_11)[1]


library(foreach)
library(doParallel)
source('FairFDA.R')

rho = 0.05
Jseq = 1:10
delta_seq <- seq(0,0.50,by=0.02) # default calibration
# delta_seq <- seq(0.08,0.48,by=0.04) # tuning strategy for calibration
mc = 500 
# mc=100 # tuning 
ratio = 0.5 # training-test split
disparity_list = c('DO', 'PD', 'DD')
for(idx_disparity in 1:length(disparity_list)){
  disparity = disparity_list[idx_disparity]
  print(paste0("disparity: ", disparity, " ratio: ", ratio))
  ncores <- 2
  cl <- makeCluster(ncores)  
  registerDoParallel(cl)
  rec <- foreach(run = 1:mc) %dopar% {
    
    xx = sampleSplit(xx_00, xx_01, xx_10, xx_11, split_ratio=ratio)
    x_00 = xx$x_00_train; x_00_test = xx$x_00_cal
    x_01 = xx$x_01_train; x_01_test = xx$x_01_cal
    x_10 = xx$x_10_train; x_10_test = xx$x_10_cal
    x_11 = xx$x_11_train; x_11_test = xx$x_11_cal
    
    n_test <- nrow(x_00_test)+nrow(x_01_test)+nrow(x_10_test)+nrow(x_11_test)
    n <- nrow(x_00)+nrow(x_01)+nrow(x_10)+nrow(x_11)
    
    phat_00 = nrow(x_00)/n
    phat_01 = nrow(x_01)/n
    phat_10 = nrow(x_10)/n
    phat_11 = nrow(x_11)/n
    
    xx = sampleSplit(x_00, x_01, x_10, x_11, split_ratio=0.5)
    x_00_train = xx$x_00_train; x_00_cal = xx$x_00_cal
    x_01_train = xx$x_01_train; x_01_cal = xx$x_01_cal
    x_10_train = xx$x_10_train; x_10_cal = xx$x_10_cal
    x_11_train = xx$x_11_train; x_11_cal = xx$x_11_cal
    
    if(disparity == 'DO'){
      s0 = -1; s1 = 1; b0 = 0; b1 = 0 ### DO
    }
    if(disparity == 'PD'){
      s0 = 0; s1 = 0; b0 = -1; b1 = 1 ### PD
    }
    if(disparity == 'DD'){
      s0 = -phat_01/(phat_00+phat_01)
      s1 = phat_11/(phat_10+phat_11)
      b0 = -phat_00/(phat_00+phat_01)
      b1 = phat_10/(phat_10+phat_11)
    }
    
    error_flda <- rep(0, 3) # FLDA; either subset for training; whole dataset; 
    D_flda <- rep(0, 3) 
    error_fairflda <- matrix(0, length(delta_seq), 3) # Fair-FLDA; either subset for training; cross-fitting; 
    D_fairflda <- matrix(0, length(delta_seq), 3) 
    tauhat_fairflda <- matrix(0, length(delta_seq), 2)
    error_fairfldac <- matrix(0, length(delta_seq), 3) # Fair-FLDAc
    D_fairfldac <- matrix(0, length(delta_seq), 3) 
    tauhat_fairfldac <- matrix(0, length(delta_seq), 2)
    
    Jhat = cvTune(x_00, x_01, x_10, x_11, Jseq, nfold=5)
    est_list = getEst(x_00, x_01, x_10, x_11, J=Jhat)
    pre_list <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list, tau=0)
    error_flda[3] = (sum(pre_list$pre_00)+sum(pre_list$pre_10)+sum(1-pre_list$pre_01)+sum(1-pre_list$pre_11))/n_test
    D_flda[3] = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list, tau=0, s0=s0, s1=s1, b0=b0, b1=b1)
    
    ### sample splitting
    Jhat1 = cvTune(x_00_train, x_01_train, x_10_train, x_11_train, Jseq, nfold=5)
    est_list1 = getEst(x_00_train, x_01_train, x_10_train, x_11_train, J=Jhat1)
    pre_list <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list1, tau=0)
    error_flda[1] = (sum(pre_list$pre_00)+sum(pre_list$pre_10)+sum(1-pre_list$pre_01)+sum(1-pre_list$pre_11))/n_test
    D_flda[1] = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list1, tau=0, s0=s0, s1=s1, b0=b0, b1=b1)
    
    Jhat2 = cvTune(x_00_cal, x_01_cal, x_10_cal, x_11_cal, Jseq, nfold=5)
    est_list2 = getEst(x_00_cal, x_01_cal, x_10_cal, x_11_cal, J=Jhat2)
    pre_list <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list2, tau=0)
    error_flda[2] = (sum(pre_list$pre_00)+sum(pre_list$pre_10)+sum(1-pre_list$pre_01)+sum(1-pre_list$pre_11))/n_test
    D_flda[2] = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list2, tau=0, s0=s0, s1=s1, b0=b0, b1=b1)
    
    for(idx_delta in 1:length(delta_seq)){
      delta <- delta_seq[idx_delta]
      eps = min(sqrt(2*log(1/rho)/n), delta) # default choice of calibration parameters
      ### tuning strategy
      # eps_seq = exp(seq(log(sqrt(2*log(1/rho)/n)), log(delta), length.out=8))
      # eps = TuneCal(x_00, x_01, x_10, x_11, J=Jhat, delta, eps_seq, s0=s0, s1=s1, b0=b0, b1=b1, rho=rho)
      
      ### Fair-FLDA
      tauhat_fairflda[idx_delta, 1] = getTauhat(x_00_cal, x_01_cal, x_10_cal, x_11_cal, est_list1, delta, eps=0, s0=s0, s1=s1, b0=b0, b1=b1)
      pre_list1_fairflda <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list1, tau=tauhat_fairflda[idx_delta, 1], s0=s0, s1=s1, b0=b0, b1=b1)
      error_fairflda[idx_delta, 1] = (sum(pre_list1_fairflda$pre_00)+sum(pre_list1_fairflda$pre_10)+sum(1-pre_list1_fairflda$pre_01)+sum(1-pre_list1_fairflda$pre_11))/n_test
      D_fairflda[idx_delta, 1] = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list1, tau=tauhat_fairflda[idx_delta, 1], s0=s0, s1=s1, b0=b0, b1=b1)
      
      ### Fair-FLDAc (with calibration)
      tauhat_fairfldac[idx_delta, 1] = getTauhat(x_00_cal, x_01_cal, x_10_cal, x_11_cal, est_list1, delta, eps=eps, s0=s0, s1=s1, b0=b0, b1=b1)
      pre_list1_fairfldac <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list1, tau=tauhat_fairfldac[idx_delta, 1], s0=s0, s1=s1, b0=b0, b1=b1)
      error_fairfldac[idx_delta, 1] = (sum(pre_list1_fairfldac$pre_00)+sum(pre_list1_fairfldac$pre_10)+sum(1-pre_list1_fairfldac$pre_01)+sum(1-pre_list1_fairfldac$pre_11))/n_test
      D_fairfldac[idx_delta, 1] = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list1, tau=tauhat_fairfldac[idx_delta, 1], s0=s0, s1=s1, b0=b0, b1=b1)
      
      
      ### Fair-FLDA
      tauhat_fairflda[idx_delta, 2] = getTauhat(x_00_train, x_01_train, x_10_train, x_11_train, est_list2, delta, eps=0, s0=s0, s1=s1, b0=b0, b1=b1)
      pre_list2_fairflda <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list2, tau=tauhat_fairflda[idx_delta, 2], s0=s0, s1=s1, b0=b0, b1=b1)
      error_fairflda[idx_delta, 2] = (sum(pre_list2_fairflda$pre_00)+sum(pre_list2_fairflda$pre_10)+sum(1-pre_list2_fairflda$pre_01)+sum(1-pre_list2_fairflda$pre_11))/n_test
      D_fairflda[idx_delta, 2] = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list2, tau=tauhat_fairflda[idx_delta, 2], s0=s0, s1=s1, b0=b0, b1=b1)
      
      ### Fair-FLDAc (with calibration)
      tauhat_fairfldac[idx_delta, 2] = getTauhat(x_00_train, x_01_train, x_10_train, x_11_train, est_list2, delta, eps=eps, s0=s0, s1=s1, b0=b0, b1=b1)
      pre_list2_fairfldac <- getPrehat(x_00_test, x_01_test, x_10_test, x_11_test, est_list2, tau=tauhat_fairfldac[idx_delta, 2], s0=s0, s1=s1, b0=b0, b1=b1)
      error_fairfldac[idx_delta, 2] = (sum(pre_list2_fairfldac$pre_00)+sum(pre_list2_fairfldac$pre_10)+sum(1-pre_list2_fairfldac$pre_01)+sum(1-pre_list2_fairfldac$pre_11))/n_test
      D_fairfldac[idx_delta, 2] = getDhat(x_00_test, x_01_test, x_10_test, x_11_test, est_list2, tau=tauhat_fairfldac[idx_delta, 2], s0=s0, s1=s1, b0=b0, b1=b1)
      
      ### cross-fitting average
      pre_00 = rbinom(nrow(x_00_test), 1, (pre_list1_fairflda$pre_00+pre_list2_fairflda$pre_00)/2)
      pre_01 = rbinom(nrow(x_01_test), 1, (pre_list1_fairflda$pre_01+pre_list2_fairflda$pre_01)/2)
      pre_10 = rbinom(nrow(x_10_test), 1, (pre_list1_fairflda$pre_10+pre_list2_fairflda$pre_10)/2)
      pre_11 = rbinom(nrow(x_11_test), 1, (pre_list1_fairflda$pre_11+pre_list2_fairflda$pre_11)/2)
      error_fairflda[idx_delta, 3] <- (sum(pre_00)+sum(pre_10)+sum(1-pre_01)+sum(1-pre_11))/n_test
      D_fairflda[idx_delta, 3] <- s0*mean(pre_01) + b0*mean(pre_00) + s1*mean(pre_11) + b1*mean(pre_10)
      # if(disparity == 'DO'){
      #   D_fairflda[idx_delta, 3] <- mean(pre_11) - mean(pre_01)
      # }
      # if(disparity == 'PD'){
      #   D_fairflda[idx_delta, 3] <- mean(pre_10) - mean(pre_00)
      # }
      # if(disparity == 'DD'){
      #   D_fairflda[idx_delta, 3] <- mean(c(pre_10,pre_11)) - mean(c(pre_00, pre_01))
      # }
      
      pre_00 = rbinom(nrow(x_00_test), 1, (pre_list1_fairfldac$pre_00+pre_list2_fairfldac$pre_00)/2)
      pre_01 = rbinom(nrow(x_01_test), 1, (pre_list1_fairfldac$pre_01+pre_list2_fairfldac$pre_01)/2)
      pre_10 = rbinom(nrow(x_10_test), 1, (pre_list1_fairfldac$pre_10+pre_list2_fairfldac$pre_10)/2)
      pre_11 = rbinom(nrow(x_11_test), 1, (pre_list1_fairfldac$pre_11+pre_list2_fairfldac$pre_11)/2)
      error_fairfldac[idx_delta, 3] <- (sum(pre_00)+sum(pre_10)+sum(1-pre_01)+sum(1-pre_11))/n_test
      D_fairfldac[idx_delta, 3] <- s0*mean(pre_01) + b0*mean(pre_00) + s1*mean(pre_11) + b1*mean(pre_10)
      # if(disparity == 'DO'){
      #   D_fairfldac[idx_delta, 3] <- mean(pre_11) - mean(pre_01)
      # }
      # if(disparity == 'PD'){
      #   D_fairfldac[idx_delta, 3] <- mean(pre_10) - mean(pre_00)
      # }
      # if(disparity == 'DD'){
      #   D_fairfldac[idx_delta, 3] <- mean(c(pre_10,pre_11)) - mean(c(pre_00, pre_01))
      # }
      
    }
    return(list(error_flda=error_flda, error_fairflda=error_fairflda, error_fairfldac=error_fairfldac,
                D_flda=D_flda, D_fairflda=D_fairflda, D_fairfldac=D_fairfldac,
                tauhat_fairflda=tauhat_fairflda, tauhat_fairfldac=tauhat_fairfldac))
    
  }
  stopCluster(cl)
  save(rec, file=paste0("result_NHANES/rec_default_", disparity, "_", ratio, "_race_age_20_50.Rdata"))
}








