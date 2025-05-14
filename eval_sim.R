library(xtable)
library(ggplot2)
source("FairFDA.R")

### model parameters
mc = 500
rho=0.05
s = 50
alpha = 2
const_snr0 = 0.8
const_snr1 = 1
pa = 0.7
p_ycona = c(0.4, 0.7)
p_01 = (1-pa)*p_ycona[1]; p_00 = (1-pa)*(1-p_ycona[1])
p_11 = pa*p_ycona[2]; p_10 = pa*(1-p_ycona[2])
beta_seq = c(2, 1.5, 0.5)
delta_seq = seq(0, 0.5, by=0.025)
### get oracle results
### snr0 \|\mu_{0,1}-\mu_{0,0}\|_{K_0}
disparity = 'DO'
if(disparity == 'DO'){
  s0 = -1; s1 = 1; b0 = 0; b1 = 0 ### DO
}
if(disparity == 'PD'){
  s0 = 0; s1 = 0; b0 = -1; b1 = 1 ### PD
}
if(disparity == 'DD'){
  s0 = -p_01/(p_00+p_01)
  s1 = p_11/(p_10+p_11)
  b0 = -p_00/(p_00+p_01)
  b1 = p_10/(p_10+p_11)
}

D_oracle = matrix(0, length(beta_seq), length(delta_seq))
error_oracle = matrix(0, length(beta_seq), length(delta_seq))
taustar_oracle = matrix(0, length(beta_seq), length(delta_seq))
for(i in 1:length(beta_seq)){
  beta = beta_seq[i]
  for(j in 1:length(delta_seq)){
    delta = delta_seq[j]
    snr0 = sqrt(sum((1:s)^(-2*beta)/(1:s)^(-alpha)))*const_snr0
    snr1 = sqrt(sum((1:s)^(-2*beta)/(1:s)^(-alpha)))*const_snr1
    oracle_list <- list(snr0=snr0, snr1=snr1, p_00=p_00, p_01=p_01, p_10=p_10, p_11=p_11)
    taustar <- getTaustar(oracle_list, delta=delta, s0=s0, s1=s1, b0=b0, b1=b1)
    taustar_oracle[i,j] = taustar
    D_oracle[i,j] = getoracleD(oracle_list, tau=taustar, s0=s0, s1=s1, b0=b0, b1=b1)
    error_oracle[i,j] = getoracleRisk(oracle_list, tau=taustar, s0=s0, s1=s1, b0=b0, b1=b1)
  }
}
D_oracle <- abs(D_oracle)
result_oracle <- NULL
for(i in 1:length(beta_seq)){
  result_oracle <- rbind(result_oracle, rbind(taustar_oracle[i,], D_oracle[i,], error_oracle[i,]))
}
# xtable(result_oracle, ncol=length(delta_seq), digits=3)

### evaluation
idx_beta = 2
beta = beta_seq[idx_beta]
data_type = 'gauss'
n_seq = c(1000, 2000, 5000)
for(idx_n in 1:length(n_seq)){
  n = n_seq[idx_n]

  load(paste0("result_cv_pa_", pa, "_pycona_", p_ycona[1], "_", p_ycona[2], "/rec", disparity, "_cv_", data_type, "_n_", n, "_alpha_", alpha, "_beta_", beta, "_constsnr0_", const_snr0, "_constsnr1_", const_snr1, ".Rdata"))
  
  error_flda <- matrix(0, mc, 3) 
  D_flda <- matrix(0, mc, 3) 
  error_fairflda <- array(0, dim=c(mc, length(delta_seq), 3))
  D_fairflda <- array(0, dim=c(mc, length(delta_seq), 3)) 
  error_fairfldac <- array(0, dim=c(mc, length(delta_seq), 3)) 
  D_fairfldac <- array(0, dim=c(mc, length(delta_seq), 3)) 
  for(i in 1:mc){
    error_flda[i,] = rec[[i]]$error_flda
    D_flda[i,] = rec[[i]]$D_flda
    error_fairflda[i,,] = rec[[i]]$error_fairflda
    D_fairflda[i,,] = rec[[i]]$D_fairflda
    error_fairfldac[i,,] = rec[[i]]$error_fairfldac
    D_fairfldac[i,,] = rec[[i]]$D_fairfldac
  }
  # result_errormean <- matrix(0, length(delta_seq), 9)
  # result_Dmean <- matrix(0, length(delta_seq), 9)
  # sd_errormean <- matrix(0, length(delta_seq), 9)
  # sd_Dmean <- matrix(0, length(delta_seq), 9)
  result_error <- matrix(0, length(delta_seq), 9)
  result_Dmedian <- matrix(0, length(delta_seq), 9)
  result_D95quan <- matrix(0, length(delta_seq), 9)
  for(idx_delta in 1:length(delta_seq)){
    delta = delta_seq[idx_delta]
    
    ### mean
    # result_errormean[idx_delta,] = c(colMeans(error_flda), colMeans(error_fairflda[, idx_delta, ]), colMeans(error_fairfldac[, idx_delta, ]))
    # result_Dmean[idx_delta,] = c(colMeans(abs(D_flda)), colMeans(abs(D_fairflda[, idx_delta, ])), colMeans(abs(D_fairfldac[, idx_delta, ])))
    # sd_errormean[idx_delta,] = c(apply(error_flda, 2, sd), apply(error_fairflda[, idx_delta, ], 2, sd), apply(error_fairfldac[, idx_delta, ], 2, sd))
    # sd_Dmean[idx_delta,] = c(apply(abs(D_flda),2,sd), apply(abs(D_fairflda[, idx_delta, ]),2,sd), apply(abs(D_fairfldac[, idx_delta, ]),2,sd))
    
    result_error[idx_delta,] = c(apply(error_flda, 2, median), apply(error_fairflda[, idx_delta, ], 2, median), apply(error_fairfldac[, idx_delta, ], 2, median))
    result_Dmedian[idx_delta,] = c(apply(abs(D_flda), 2, median), apply(abs(D_fairflda[, idx_delta, ]), 2, median), apply(abs(D_fairfldac[, idx_delta, ]), 2, median))
    result_D95quan[idx_delta,] = c(apply(abs(D_flda),2,quantile, probs=1-rho), apply(abs(D_fairflda[, idx_delta, ]),2,quantile, probs=1-rho), apply(abs(D_fairfldac[, idx_delta, ]),2,quantile, probs=1-rho))
  }
  
  ### cross-fitting average
  flda <- rbind(result_Dmedian[,3], result_D95quan[,3], result_error[,3])
  fair_flda <- rbind(result_Dmedian[,6], result_D95quan[,6], result_error[,6])
  fair_fldac <- rbind(result_Dmedian[,9], result_D95quan[,9], result_error[,9])
  result <- rbind(flda, fair_flda, fair_fldac)
  # xtable(result, nrow=9, ncol=length(delta_seq), digits=3)
  
  ######### figure illustration #####################
  
  ### delta-error plot
  error <- c(result[3,], result[6,], result[9,], error_oracle[idx_beta,])
  data <- cbind(rep(delta_seq, 4), error)
  colnames(data) <- c('delta', 'error')
  data = as.data.frame(data)
  data$group <- rep(c('A', 'B', 'C', 'D'), each=length(delta_seq))
  
  w<-10
  h<-9
  pdf_file <- paste0('figs/error_disparity_', disparity, '_cv_', data_type, '_n_', n, '_beta_', beta, '.pdf')
  pdf(file=pdf_file, width=w, height=h)
  
  
  p <- ggplot(data,aes(x=delta)) + theme(panel.grid.major =element_blank(),
                                         panel.grid.minor = element_blank(),
                                         panel.background = element_rect(colour = 'black', fill='white'),
                                         axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(hjust=0.5))
  
  p <- p + geom_line(data=data, aes(y=error, colour=group, linetype=group), lwd=2)
  p <- p + geom_point(data=data, aes(y=error, color=group, shape=group), size=4)
  
  p <- p + scale_x_continuous(limits = c(0, 0.5), breaks=c(0, 0.2, 0.4))
  p <- p + scale_y_continuous(limits = c(0.12, 0.25))
  # p <- p + scale_y_continuous(limits = c(0., 0.005), breaks=c(0, 0.002, 0.004))
  p <- p + scale_colour_manual(name='group', 
                               values=c('A'='#f47a00', 'B'='#00b0be', 'C'='#f45f74', 'D'='#d31f11'),
                               breaks=c('A', 'B', 'C'),
                               labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
  p <- p + scale_linetype_manual(name='group',
                                 values=c('A'='blank', 'B'='blank', 'C'='blank', 'D'='solid'),
                                 breaks=c('A', 'B', 'C'),
                                 labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
  p <- p + scale_shape_manual(name='group',
                              values=c('A'=16, 'B'=8, 'C'=2, 'D'=NA),
                              breaks=c('A', 'B', 'C'),
                              labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
  p <- p + labs(y = 'Error', x = expression(delta))
  p <- p + theme(axis.title = element_text(size=70)) + 
    theme(axis.text = element_text(size=70)) + 
    theme(legend.position = 'none')
  
  print(p)  
  dev.off()  
  
  for(idx_disparity_type in 1:2){
    disparity_type = c('median', '95quantile')[idx_disparity_type]
    
    if(disparity_type == 'median'){
      Dmedian <- c(result[1,], result[4,], result[7,], D_oracle[idx_beta,]) 
      if(disparity == 'DO'){
        # xlab_title = expression(bar(U)[DO])
        xlab_title = expression(U[DO*","*50])
      }
      if(disparity == 'PD'){
        # xlab_title = expression(bar(U)[PD])
        xlab_title = expression(U[PD*","*50])
      }
      if(disparity == 'DD'){
        # xlab_title = expression(bar(U)[DD])
        xlab_title = expression(U[DD*","*50])
      }
    }
    if(disparity_type == '95quantile'){
      Dmedian <- c(result[2,], result[5,], result[8,], D_oracle[idx_beta,]) # quantile disparity
      if(disparity == 'DO'){
        xlab_title <- expression(U[DO*","*95])
      }
      if(disparity == 'PD'){
        xlab_title <- expression(U[PD*","*95])
      }
      if(disparity == 'DD'){
        xlab_title <- expression(U[DD*","*95])
      }
    }
    error <- c(result[3,], result[6,], result[9,], error_oracle[idx_beta,])
    data <- cbind(Dmedian, error)
    colnames(data) <- c('UD', 'error')
    data <- as.data.frame(data)
    data$group <- rep(c('A', 'B', 'C', 'D'), each=length(delta_seq))

    w<-10
    h<-9
    pdf_file <- paste0('figs/tradeoff_', disparity_type, '_disparity_', disparity, '_cv_', data_type, '_n_', n, '_beta_', beta, '.pdf')
    pdf(file=pdf_file, width=w, height=h)

    p <- ggplot(data,aes(x=UD)) + theme(panel.grid.major =element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.background = element_rect(colour = 'black', fill='white'),
                                        axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(hjust=0.5))

    p <- p + geom_line(data=data, aes(y=error, colour=group, linetype=group), lwd=2)
    p <- p + geom_point(data=data, aes(y=error, color=group, shape=group), size=4)

    p <- p + scale_x_continuous(limits = c(0, 0.4), breaks=c(0, 0.15, 0.3))
    p <- p + scale_y_continuous(limits = c(0.12, 0.25))
    p <- p + scale_colour_manual(name='group',
                                values=c('A'='#f47a00', 'B'='#00b0be', 'C'='#f45f74', 'D'='#d31f11'),
                                breaks=c('A', 'B', 'C'),
                                labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
    p <- p + scale_linetype_manual(name='group',
                                  values=c('A'='blank', 'B'='blank', 'C'='blank', 'D'='solid'),
                                  breaks=c('A', 'B', 'C'),
                                  labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
    p <- p + scale_shape_manual(name='group',
                                values=c('A'=16, 'B'=8, 'C'=2, 'D'=NA),
                                breaks=c('A', 'B', 'C'),
                                labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
    # p <- p + ggtitle(paste0(expression(beta),'=', beta, ', n=', n)) + theme(plot.title = element_text(size=70))
    p <- p + labs(y = 'Error', x = xlab_title) +
      theme(axis.title = element_text(size=70)) +
      theme(axis.text = element_text(size=70)) +
      theme(legend.position = 'none')
    print(p)
    dev.off()
    
    
    ### disparity Q-Q plot
    if(disparity_type == 'median'){
      data = c(result[1,], result[4,], result[7,]) 
      if(disparity == 'DO'){
        ylab_title = expression(U[DO*","*50])
      }
      if(disparity == 'PD'){
        ylab_title = expression(U[PD*","*50])
      }
      if(disparity == 'DD'){
        ylab_title = expression(U[DD*","*50])
      }
    }
    if(disparity_type == '95quantile'){
      data = c(result[2,], result[5,], result[8,]) # quantile disparity
      main_title = expression(U[D*","*95])
      if(disparity == 'DO'){
        ylab_title = expression(U[DO*","*95])
      }
      if(disparity == 'PD'){
        ylab_title = expression(U[PD*","*95])
      }
      if(disparity == 'DD'){
        ylab_title = expression(U[DD*","*95])
      }
    }
    data = cbind(rep(delta_seq, 3), data)
    colnames(data) <- c('delta', 'UD')
    data = as.data.frame(data)
    data$group <- rep(c('A', 'B', 'C'), each=length(delta_seq))
    
    w<-10
    h<-9
    pdf_file <- paste0('figs/qq_', disparity_type, '_disparity_', disparity, '_cv_', data_type, '_n_', n, '_beta_', beta, '.pdf')
    pdf(file=pdf_file, width=w, height=h)
    
    
    p <- ggplot(data,aes(x=delta)) + theme(panel.grid.major =element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.background = element_rect(colour = 'black', fill='white'),
                                           axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(hjust=0.5))

    p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", size = 1)
    
    p <- p + geom_point(data=data, aes(y=UD, color=group, shape=group), size=4)
    
    p <- p + scale_x_continuous(limits = c(0, 0.5), breaks=c(0, 0.2, 0.4))
    p <- p + scale_y_continuous(limits = c(0, 0.5), breaks=c(0, 0.2, 0.4))
    # p <- p + scale_y_continuous(limits = c(0, 0.02), breaks=c(0, 0.01, 0.02))
    p <- p + scale_colour_manual(name='group', 
                                 values=c('A'='#f47a00', 'B'='#00b0be', 'C'='#f45f74'),
                                 breaks=c('A', 'B', 'C'),
                                 labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
    p <- p + scale_shape_manual(name='group',
                                values=c('A'=16, 'B'=8, 'C'=2),
                                breaks=c('A', 'B', 'C'),
                                labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
    # p <- p + ggtitle(paste0(expression(beta),'=', beta, ', n=', n)) + theme(plot.title = element_text(size=70))
    p <- p + labs(y = ylab_title, x = expression(delta))
    p <- p + theme(axis.title = element_text(size=70)) + 
      theme(axis.text = element_text(size=70)) + 
      theme(legend.position = 'none')
    
    print(p)  
    dev.off()  
    
  }
}









