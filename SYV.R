##############################################
######## SYV (continuous) Application ########
##############################################
library(ggplot2)
library(directlabels)
library(cowplot)

############################################
# Constant treatment effect
############################################

# Function of sample size calculations for longitudinal IRGT trials, when a nominal type I error rate is fixed
# at 5%, with continuous outcomes under the canonical identity link function, assuming no interaction between
# time and intervention
conSIMULATE <- function(K,t,phi,b,alpha,gamma=0.15){
  qc <- K/(K+1)
  lambda_t4 <- 1+(K-1)*alpha[1]+(t-1)*(K-1)*alpha[2]+(t-1)*alpha[3]
  lambda_c2 <- 1+(t-1)*alpha[3]
  sigma2 <- phi/t*(lambda_c2/qc+lambda_t4/qc)
  I <- (qnorm(0.05/2)+qnorm(gamma))^2/b^2*sigma2
  I <- ceiling(I)
  if (I%%(1+K) != 0){
    I <- I + ((K+1)-I%%(K+1))
  }
  It <- (qt(0.05/2, df=I-2)+qt(gamma, I-2))^2/b^2*sigma2
  while (I<It){
    I <- I + (K+1)
    It <- (qt(0.05/2, df=I-2)+qt(gamma, I-2))^2/b^2*sigma2
  }
  t_power <- pt(qt(0.05/2, df=I-2)+abs(b)*sqrt(I/sigma2), df=I-2)
  N <- 2*K*I/(K+1)
  return(data.frame(b=b, alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, phi=phi, a_power=t_power))
}

# Function of power calculations for longitudinal IRGT trials, with continuous outcomes under the canonical
# identity link function, assuming no interaction between time and intervention
conPOWER <- function(N,I,K,t,phi,b,alpha){
  qc <- K/(K+1)
  lambda_t4 <- 1+(K-1)*alpha[1]+(t-1)*(K-1)*alpha[2]+(t-1)*alpha[3]
  lambda_c2 <- 1+(t-1)*alpha[3]
  sigma2 <- phi/t*(lambda_c2/qc+lambda_t4/qc)
  
  t_power <- pt(qt(0.05/2, df=I-2)+abs(b)*sqrt(I/sigma2), df=I-2)
  return(data.frame(alpha1=alpha[2], alpha2=alpha[3], t_test=t_power))
}

# Function of power contours for longitudinal IRGT trials, with with continuous outcomes under the canonical
# identity link function, assuming no interaction between time and intervention
conPLOT <- function(N,I,K,t,phi,b,a0,a1_range,a2_range,pc=0.5,randomization=4){
  con_plot<-NULL
  syv<-NULL
  for (a1 in seq(a1_range[1], a1_range[2], 0.001)){
    for (a2 in seq(a2_range[1], a2_range[2], 0.0005)){
      syv<-rbind(syv, conPOWER(N,I,K,t,phi,b,c(a0, a1, a2)))
    }
  }
  fig<-ggplot()  +
    theme_bw() +
    ggtitle(bquote(alpha[0] == ~.(a0))) +
    xlab(expression(alpha[1])) +
    ylab(expression(alpha[2])) +
    ylim(c(a2_range[1], a2_range[2])) +
    scale_x_continuous(breaks=seq(a1_range[1], a1_range[2], 0.01)) +
    stat_contour(data = syv, aes(x = alpha1, y = alpha2, z = t_test, colour = ..level..), 
                 breaks = round(quantile(syv$t_test, seq(0, 1, 0.1)), 2), size = 1) +
    scale_color_continuous(name = "Power") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification=c(1, 0), legend.position=c(1, 0))
  con_plot<-direct.label(fig, "bottom.pieces")
  return(con_plot)
}

# ---------------------------------------------------------------------------------


### results
# Effect size of 0.3 standard deviation (SD)
ss <- conSIMULATE(8, 3, 1, 0.3, c(0.04, 0.05, 0.8)) # N = 416, I = 234, K = 8, power = 0.85128
ss
conPOWER(416, 234, 8, 3, 1, 0.3, c(0.04, 0.05, 0.8))

### power contours
# Effect size of 0.3 standard deviation (SD) and N = 416, I = 234, K = 8
syv_1_1 <- conPLOT(416, 234, 8, 3, 1, 0.3, 0.03, c(0.02, 0.1), c(0.7, 0.9))
syv_1_2 <- conPLOT(416, 234, 8, 3, 1, 0.3, 0.04, c(0.02, 0.1), c(0.7, 0.9))
syv_1_3 <- conPLOT(416, 234, 8, 3, 1, 0.3, 0.05, c(0.02, 0.1), c(0.7, 0.9))
Figure_5 <- plot_grid(syv_1_1, syv_1_2, syv_1_3, ncol=3)
ggsave("Figure_syv.pdf", height = 5, width = 15)



############################################
# Linear time by treatment interaction
############################################

# Function of sample size calculations for longitudinal IRGT trials, when a nominal type I error rate is fixed
# at 5%, with continuous outcomes under the canonical identity link function, assuming a linear time by treatment
# interaction
conSIMULATE <- function(K,t,phi,b,alpha,gamma=0.15){
  qc <- K/(K+1)
  lambda_t3 <- 1+(K-1)*(alpha[1]-alpha[2])-alpha[3]
  lambda_t4 <- 1+(K-1)*alpha[1]+(t-1)*(K-1)*alpha[2]+(t-1)*alpha[3]
  lambda_c1 <- 1-alpha[3]
  lambda_c2 <- 1+(t-1)*alpha[3]
  
  V3 <- 12*phi/(t*(t^2-1))*(lambda_c1/qc+lambda_t3/qc)
  V4 <- phi/t*(lambda_c2/qc+lambda_t4/qc)
  Va <- V3*(t+1)^2/4 + V4
  Vb <- -V3*(t+1)/2
  Vd <- V3
  Vsub_inv <- solve(matrix(c(Va,Vb,Vb,Vd),2,2))
  I <- K+1
  F_power <- 0
  while (F_power<(1-gamma)){
    I <- I + (K+1)
    lambda <- I*(b%*%Vsub_inv%*%b)
    fc <- qf(1-0.05, df1=2, df2=I-3)
    F_power <- 1-pf(fc, df1=2, df2=I-3, ncp=lambda)
  }
  N <- 2*K*I/(K+1)
  return(data.frame(b=b[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, phi=phi, a_power=F_power))
}

# ---------------------------------------------------------------------------------

### results
b3 <- c(0.3, 0.1)
ss <- conSIMULATE(8, 3, 1, b3, c(0.04, 0.05, 0.8)) # N = 128, I = 72, K = 8, power = 0.8574815
ss



############################################
# Categorical time by treatment interaction
############################################

# Function of sample size calculations for longitudinal IRGT trials, when a nominal type I error rate is fixed
# at 5%, with continuous outcomes under the canonical identity link function, assuming categorical time by treatment
# interactions
conSIMULATE <- function(K,t,phi,b,alpha,gamma=0.15){
  qc <- K/(K+1)
  lambda_t3 <- 1+(K-1)*(alpha[1]-alpha[2])-alpha[3]
  lambda_t4 <- 1+(K-1)*alpha[1]+(t-1)*(K-1)*alpha[2]+(t-1)*alpha[3]
  lambda_c1 <- 1-alpha[3]
  lambda_c2 <- 1+(t-1)*alpha[3]
  Va <- phi*((lambda_c1+lambda_t3)/qc)
  Vb <- phi/t*((lambda_c2-lambda_c1)/qc+(lambda_t4-lambda_t3)/qc)
  Vc <- 1/Va
  Vd <- -Vb/(Va*(Va+t*Vb))
  Vsub_inv <- Vc*diag(t) + Vd*matrix(1,t,t)
  I <- K+1
  F_power <- 0
  while (F_power<(1-gamma)){
    I <- I + (K+1)
    lambda <- I*(b%*%Vsub_inv%*%b)
    fc <- qf(1-0.05, df1=t, df2=I-t-1)
    F_power <- 1-pf(fc, df1=t, df2=I-t-1, ncp=lambda)
  }
  N <- 2*K*I/(K+1)
  return(data.frame(b=b[t], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, phi=phi, a_power=F_power))
}

# ---------------------------------------------------------------------------------

### results
b2j <- c(0.5, 0.3, 0.1)
ss <- conSIMULATE(8, 3, 1, b2j, c(0.04, 0.05, 0.8)) # N = 96, I = 54, K = 8, power = 0.8620326
ss


