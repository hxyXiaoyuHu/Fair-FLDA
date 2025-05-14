
### @pa: class probability P(A=1)
### @p_ycona: a vector containing P(Y=1|A=0) and P(Y=1|A=1)
### @p_01: class probability P(A=0, Y=1)
### @x_{01}: data from the group (A=0, Y=1)

library(foreach)
library(doParallel)

source("FairFDA.R")

pa = 0.7
p_ycona = c(0.4, 0.7)
# p_ycona = c(0.5, 0.5)
p_01 = (1-pa)*p_ycona[1]; p_00 = (1-pa)*(1-p_ycona[1])
p_11 = pa*p_ycona[2]; p_10 = pa*(1-p_ycona[2])
n_test = 5000
tseq = seq(0,1,by=0.01)
mc = 500
s = 50
rho=0.05
alpha = 2
const_snr0 = 0.8
const_snr1 = 1
delta=0.05
Jseq = 1:10
delta_seq = seq(0, 0.5, by=0.025) # result_cv
n_seq = c(1000, 2000, 5000)
beta_seq = c(2, 1.5, 0.5)
disparity_list = c('DO', 'PD', 'DD')
data_type = 'gauss'
# data_type = 'unif'
for(idx_disparity in 1:3){
  disparity = disparity_list[idx_disparity]
  
  for(idx_n in 1:3){
    n = n_seq[idx_n]
    for(idx_beta in 1:3){
      beta = beta_seq[idx_beta]
      print(paste0('disparity: ', disparity, ' data_type: ', data_type, ' n: ', n, ' beta: ', beta))
      
      ncores <- 5
      cl <- makeCluster(ncores)  
      registerDoParallel(cl)
      rec <- foreach(run = 1:mc) %dopar% {
        
        error_flda <- rep(0, 3) # FLDA; either subset for training; whole dataset; 
        D_flda <- rep(0, 3) 
        error_fairflda <- matrix(0, length(delta_seq), 3) # Fair-FLDA; either subset for training; cross-fitting; 
        D_fairflda <- matrix(0, length(delta_seq), 3) 
        tauhat_fairflda <- matrix(0, length(delta_seq), 2)
        error_fairfldac <- matrix(0, length(delta_seq), 3) # Fair-FLDAc
        D_fairfldac <- matrix(0, length(delta_seq), 3) 
        tauhat_fairfldac <- matrix(0, length(delta_seq), 2)
        
        x = dataGen(n, pa, p_ycona, tseq, alpha, beta, type=data_type, const_snr0=const_snr0, const_snr1=const_snr1)
        x_00 = x$x_00; x_01 = x$x_01
        x_10 = x$x_10; x_11 = x$x_11
        
        phat_00 = nrow(x_00)/n
        phat_01 = nrow(x_01)/n
        phat_10 = nrow(x_10)/n
        phat_11 = nrow(x_11)/n
        
        x_test = dataGen(n_test, pa, p_ycona, tseq, alpha, beta, type=data_type, const_snr0=const_snr0, const_snr1=const_snr1)
        x_00_test = x_test$x_00; x_01_test = x_test$x_01
        x_10_test = x_test$x_10; x_11_test = x_test$x_11
        
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
          delta = delta_seq[idx_delta]
          eps = min(sqrt(2*log(1/rho)/n), delta) # calibration constant for Fair-FLDAc
          
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
      save(rec, file=paste0("result_cv_pa_", pa, "_pycona_", p_ycona[1], "_", p_ycona[2], "/rec", disparity, "_cv_", data_type, "_n_", n, "_alpha_", alpha, "_beta_", beta, "_constsnr0_", const_snr0, "_constsnr1_", const_snr1, ".Rdata"))
    }
  }
}




