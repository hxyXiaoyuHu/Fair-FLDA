library(xtable)
library(ggplot2)
source("FairFDA.R")

rho=0.05
ratio = 0.5
n = 1627 # training sample size

mc = 500
delta_seq <- seq(0,0.50,by=0.02)
# mc = 100 # tuning for calibration
# delta_seq <- seq(0.08,0.48,by=0.04) # tuning for calibration

disparity_list = c('DO', 'PD', 'DD')
for(idx_disparity in 1:3){
  disparity = disparity_list[idx_disparity]
  load(paste0("result_NHANES/rec_default_", disparity, "_", ratio, "_race_age_20_50.Rdata"))
  
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
  error <- c(result[3,], result[6,], result[9,])
  data <- cbind(rep(delta_seq, 3), error)
  colnames(data) <- c('delta', 'error')
  data = as.data.frame(data)
  data$group <- rep(c('A', 'B', 'C'), each=length(delta_seq))
  
  w<-10
  h<-9
  pdf_file <- paste0('figs/error_disparity_', disparity, "_", ratio, "_race_age_20_50.pdf")
  pdf(file=pdf_file, width=w, height=h)
  
  if(disparity=='DO'|disparity=='PD'){
    p <- ggplot(data,aes(x=delta)) + theme(panel.grid.major =element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.background = element_rect(colour = 'black', fill='white'),
                                           axis.line = element_line(colour = "black"),
                                           axis.text.x = element_blank())
  }
  if(disparity=='DD'){
    p <- ggplot(data,aes(x=delta)) + theme(panel.grid.major =element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.background = element_rect(colour = 'black', fill='white'),
                                           axis.line = element_line(colour = "black"))
  }
  
  p <- p + geom_point(data=data, aes(y=error, color=group, shape=group), size=4)
  
  p <- p + scale_x_continuous(limits = c(0, 0.5), breaks=c(0, 0.2, 0.4))
  p <- p + scale_y_continuous(limits = c(0.27, 0.32), breaks=round(seq(0.28, 0.32, length.out=3), 2))
  p <- p + scale_colour_manual(name='group', 
                               values=c('A'='#f47a00', 'B'='#00b0be', 'C'='#f45f74'),
                               breaks=c('A', 'B', 'C'),
                               labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))

  p <- p + scale_shape_manual(name='group',
                              values=c('A'=16, 'B'=8, 'C'=2),
                              breaks=c('A', 'B', 'C'),
                              labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
  
  if(disparity=='DO'){
    p <- p + ggtitle('Error') + theme(plot.title = element_text(size=70))+
      theme(plot.title = element_text(hjust=0.5))
  }
  p <- p + theme(axis.text = element_text(size=70)) + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(legend.position = 'none')
  
  print(p)  
  dev.off()  
  
  
  for(idx_disparity_type in 1:2){
    disparity_type = c('median', '95quantile')[idx_disparity_type]
    
    ### disparity Q-Q plot
    if(disparity_type == 'median'){
      
      data = c(result[1, ], result[4,], result[7,]) # median disparity
      main_title = expression(U[D*","*50])
    }
    if(disparity_type == '95quantile'){
      data = c(result[2, ], result[5,], result[8,]) # 95quantile disparity
      main_title = expression(U[D*","*95])
    }
    data = cbind(rep(delta_seq, 3), data)
    colnames(data) <- c('delta', 'UD')
    data = as.data.frame(data)
    data$group <- rep(c('A', 'B', 'C'), each=length(delta_seq))
    
    w<-10
    h<-9
    pdf_file <- paste0('figs/qq_', disparity_type, '_disparity_', disparity, "_", ratio, "_race_age_20_50.pdf")
    pdf(file=pdf_file, width=w, height=h)
    
    if(disparity=='DO'|disparity=='PD'){
      p <- ggplot(data,aes(x=delta)) + theme(panel.grid.major =element_blank(),
                                             panel.grid.minor = element_blank(),
                                             panel.background = element_rect(colour = 'black', fill='white'),
                                             axis.line = element_line(colour = "black"),
                                             axis.text.x = element_blank())
    }
    if(disparity=='DD'){
      p <- ggplot(data,aes(x=delta)) + theme(panel.grid.major =element_blank(),
                                             panel.grid.minor = element_blank(),
                                             panel.background = element_rect(colour = 'black', fill='white'),
                                             axis.line = element_line(colour = "black"))
    }
    
    p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", size = 1)

    p <- p + geom_point(data=data, aes(y=UD, color=group, shape=group), size=4)
    
    p <- p + scale_x_continuous(limits = c(0, 0.5), breaks=c(0, 0.2, 0.4))
    p <- p + scale_y_continuous(limits = c(0, 0.55), breaks = c(0.1, 0.3, 0.5))
    p <- p + scale_colour_manual(name='group', 
                                 values=c('A'='#f47a00', 'B'='#00b0be', 'C'='#f45f74'),
                                 breaks=c('A', 'B', 'C'),
                                 labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
    p <- p + scale_shape_manual(name='group',
                                values=c('A'=16, 'B'=8, 'C'=2),
                                breaks=c('A', 'B', 'C'),
                                labels=c('FLDA', 'Fair-FLDA', expression("Fair-FLDA"[c])))
    if(disparity=='DO'){
      p <- p + ggtitle(main_title) + theme(plot.title = element_text(size=70)) +
        theme(plot.title = element_text(hjust=0.5))
    }
    p <- p + theme(axis.text = element_text(size=70)) + 
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      theme(legend.position = 'none')
    print(p)  
    dev.off()  
  }
}
