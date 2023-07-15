############################################################################################
############################ Design continuous simulations #################################
############################################################################################
# N: number of individuals
# I: number of clusters
# K: number of individuals per cluster in the intervention arm
# t: number of time points
# phi: dispersion parameter
# b: effect size
# alpha - correlations
#         alpha[1] - within-period correlation
#         alpha[2] - inter-period correlation
#         alpha[3] - within-individual correlation
# qc: percentage of clusters assigned to the control arm
# gamma: type II error rate; power = 1-gamma
############################################################################################
library(openxlsx)

# Functions of sample size calculations for longitudinal IRGTs, when a nominal type I error
# rate is fixed at 5%, with continuous outcomes under the canonical identity link function
# (The following functions differ from the actual simulation program.
# These are only used to calculate the power analytically.)

###################################################
######### (1) No time effect #########
###################################################
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

# Analytical power
aPower <- NULL
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.03, 0.015, 0.2))) # N=180, I=99
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.1, 0.05, 0.2))) # N=240, I=132
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.01, 0.005, 0.4))) # N=200, I=110
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.03, 0.015, 0.2))) # N=160, I=88
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.1, 0.05, 0.2))) # N=220, I=121
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.01, 0.005, 0.4))) # N=180, I=99



###################################################
######### (2) Linear time effect #########
# without interaction between time and intervention
###################################################
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
  It <- (qt(0.05/2, df=I-3)+qt(gamma, I-3))^2/b^2*sigma2
  while (I<It){
    I <- I + (K+1)
    It <- (qt(0.05/2, df=I-3)+qt(gamma, I-3))^2/b^2*sigma2
  }
  t_power <- pt(qt(0.05/2, df=I-3)+abs(b)*sqrt(I/sigma2), df=I-3)
  N <- 2*K*I/(K+1)
  return(data.frame(b=b, alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, phi=phi, a_power=t_power))
}

# Analytical power
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.03, 0.015, 0.2))) # N=180, I=99
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.1, 0.05, 0.2))) # N=240, I=132
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.01, 0.005, 0.4))) # N=200, I=110
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.03, 0.015, 0.2))) # N=160, I=88
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.1, 0.05, 0.2))) # N=220, I=121
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.01, 0.005, 0.4))) # N=180, I=99



###################################################
######### (3) Categorical time effect #########
# without interaction between time and intervention
###################################################
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

# Analytical power
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.03, 0.015, 0.2))) # N=180, I=99
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.1, 0.05, 0.2))) # N=240, I=132
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, 0.35, c(0.01, 0.005, 0.4))) # N=200, I=110
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.03, 0.015, 0.2))) # N=160, I=88
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.1, 0.05, 0.2))) # N=220, I=121
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, 0.35, c(0.01, 0.005, 0.4))) # N=180, I=99



###################################################
######### (4) Linear time effect #########
# with interaction between time and intervention
# overall test
###################################################
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

# Analytical power
b3 <- c(0.18, 0.1)
b4 <- c(0.15, 0.07)
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, b3, c(0.03, 0.015, 0.2))) # N=180, I=99
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, b3, c(0.1, 0.05, 0.2))) # N=220, I=121
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, b3, c(0.01, 0.005, 0.4))) # N=180, I=99
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, b4, c(0.03, 0.015, 0.2))) # N=200, I=110
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, b4, c(0.1, 0.05, 0.2))) # N=260, I=143
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, b4, c(0.01, 0.005, 0.4))) # N=220, I=121



###################################################
######### (5) Categorical time effect #########
# witht interaction between time and intervention
# global test
###################################################
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

# Analytical power
b2j_in <- 0.3*c(0.6^1, 0.6^2, 0.6^3, 0.6^4)
b2j <- cumsum(b2j_in)
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, b2j[1:3], c(0.03, 0.015, 0.2))) # N=340, I=187
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, b2j[1:3], c(0.1, 0.05, 0.2))) # N=460, I=253
aPower <- rbind(aPower, conSIMULATE(10, 3, 1, b2j[1:3], c(0.01, 0.005, 0.4))) # N=360, I=198
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, b2j[1:4], c(0.03, 0.015, 0.2))) # N=260, I=143
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, b2j[1:4], c(0.1, 0.05, 0.2))) # N=360, I=198
aPower <- rbind(aPower, conSIMULATE(10, 4, 1, b2j[1:4], c(0.01, 0.005, 0.4))) # N=280, I=154

write.xlsx(aPower, file = "conResults/pred_power_con.xlsx", row.names = FALSE)


