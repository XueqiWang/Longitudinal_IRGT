############################################################################################
############################## Design binary simulations ###################################
############################################################################################
# N: number of individuals
# I: number of clusters
# K: number of individuals per cluster in the intervention arm
# t: number of time points
# beta: regression parameters
# alpha - correlations
#         alpha[1] - within-period correlation
#         alpha[2] - inter-period correlation
#         alpha[3] - within-individual correlation
# qc: percentage of clusters assigned to the control arm
# gamma: type II error rate; power = 1-gamma
############################################################################################
library(openxlsx)

# Functions of sample size calculations for longitudinal IRGTs, when a nominal type I error
# rate is fixed at 5%, with binary outcomes under the canonical identity link function
# (The following functions differ from the actual simulation program.
# These are only used to calculate the power analytically.)

###################################################
######### (1) No time effect #########
###################################################
binSIMULATE <- function(K,t,beta,alpha,gamma=0.15){
  gmu0 <- c(c(1,0)%*%beta)
  mu0 <- plogis(gmu0)
  
  gmu1 <- c(c(1,1)%*%beta)
  mu1 <- plogis(gmu1)
  
  b <- beta[2]
  qc <- K/(K+1)
  lambda_t4 <- 1+(K-1)*alpha[1]+(t-1)*(K-1)*alpha[2]+(t-1)*alpha[3]
  lambda_c2 <- 1+(t-1)*alpha[3]
  sigma2 <- 1/t*(lambda_c2/(qc*mu0*(1-mu0))+lambda_t4/(qc*mu1*(1-mu1)))
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
  return(data.frame(b=beta[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, a_power=t_power))
}

binSIMULATE <- function(K,t,beta,alpha,gamma=0.15){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  alpha2 <- alpha[3]
  
  R0 <- (1-alpha2)*diag(t) + alpha2*matrix(1,t,t)
  invR0 <- solve(R0)
  X0 <- cbind(matrix(1,t,1), matrix(0,t,1))
  gmu0 <- c(X0%*%beta)
  mu0 <- plogis(gmu0)
  W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
  
  bm1 <- (1-alpha0+alpha1-alpha2)*diag(t*K)
  bm2 <- (alpha2-alpha1)*kronecker(matrix(1,t,t), diag(K))
  bm3 <- (alpha0-alpha1)*kronecker(diag(t), matrix(1,K,K))
  bm4 <- alpha1*matrix(1,t*K,t*K)
  R1 <- bm1+bm2+bm3+bm4
  invR1 <- solve(R1)
  X1 <- cbind(matrix(1,t*K,1), matrix(1,t*K,1))
  gmu1 <- c(X1%*%beta)
  mu1 <- plogis(gmu1)
  W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
  
  I <- 0
  t_power <- 0
  while (t_power<(1-gamma)){
    I <- I + (K+1)
    N <- 2*K*I/(K+1)
    Omega <- (N/2)*(t(X0)%*%W0%*%X0) + (N/20)*(t(X1)%*%W1%*%X1)
    vardelta <- solve(Omega)[2,2]
    t_power <- pt(qt(0.05/2, df=I-2)+abs(beta[2])/sqrt(vardelta), df=I-2)
  }
  
  return(data.frame(b=beta[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, a_power=t_power))
}


# Analytical power
aPower <- NULL
beta <- c(0, log(2))
aPower <- rbind(aPower, binSIMULATE(10, 3, beta, c(0.03, 0.015, 0.2))) # N=200, I=110
aPower <- rbind(aPower, binSIMULATE(10, 3, beta, c(0.1, 0.05, 0.2))) # N=260, I=143
aPower <- rbind(aPower, binSIMULATE(10, 3, beta, c(0.01, 0.005, 0.4))) # N=220, I=121
aPower <- rbind(aPower, binSIMULATE(10, 4, beta, c(0.03, 0.015, 0.2))) # N=160, I=88
aPower <- rbind(aPower, binSIMULATE(10, 4, beta, c(0.1, 0.05, 0.2))) # N=240, I=132
aPower <- rbind(aPower, binSIMULATE(10, 4, beta, c(0.01, 0.005, 0.4))) # N=200, I=110



###################################################
######### (2) Linear time effect #########
# without interaction between time and intervention
###################################################
binSIMULATE <- function(K,t,beta,alpha,gamma=0.15){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  alpha2 <- alpha[3]
  
  R0 <- (1-alpha2)*diag(t) + alpha2*matrix(1,t,t)
  invR0 <- solve(R0)
  X0 <- cbind(matrix(1,t,1), seq(1,t,1), matrix(0,t,1))
  gmu0 <- c(X0%*%beta)
  mu0 <- plogis(gmu0)
  W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
  
  bm1 <- (1-alpha0+alpha1-alpha2)*diag(t*K)
  bm2 <- (alpha2-alpha1)*kronecker(matrix(1,t,t), diag(K))
  bm3 <- (alpha0-alpha1)*kronecker(diag(t), matrix(1,K,K))
  bm4 <- alpha1*matrix(1,t*K,t*K)
  R1 <- bm1+bm2+bm3+bm4
  invR1 <- solve(R1)
  X1 <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(1,t,1)), matrix(1,K,1))
  gmu1 <- c(X1%*%beta)
  mu1 <- plogis(gmu1)
  W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
  
  I <- 0
  t_power <- 0
  while (t_power<(1-gamma)){
    I <- I + (K+1)
    N <- 2*K*I/(K+1)
    Omega <- (N/2)*(t(X0)%*%W0%*%X0) + (N/20)*(t(X1)%*%W1%*%X1)
    vardelta <- solve(Omega)[3,3]
    t_power <- pt(qt(0.05/2, df=I-2)+abs(beta[3])/sqrt(vardelta), df=I-2)
  }
  
  return(data.frame(b=beta[3], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, a_power=t_power))
}

# Analytical power
beta <- c(log(0.7), log(1.2), log(1.9))
aPower <- rbind(aPower, binSIMULATE(10, 3, beta, c(0.03, 0.015, 0.2))) # N=220, I=121
aPower <- rbind(aPower, binSIMULATE(10, 3, beta, c(0.1, 0.05, 0.2))) # N=300, I=165
aPower <- rbind(aPower, binSIMULATE(10, 3, beta, c(0.01, 0.005, 0.4))) # N=240, I=132
aPower <- rbind(aPower, binSIMULATE(10, 4, beta, c(0.03, 0.015, 0.2))) # N=200, I=110
aPower <- rbind(aPower, binSIMULATE(10, 4, beta, c(0.1, 0.05, 0.2))) # N=280, I=154
aPower <- rbind(aPower, binSIMULATE(10, 4, beta, c(0.01, 0.005, 0.4))) # N=240, I=132



###################################################
######### (3) Categorical time effect #########
# without interaction between time and intervention
###################################################
binSIMULATE <- function(K,t,beta,alpha,gamma=0.15){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  alpha2 <- alpha[3]
  
  R0 <- (1-alpha2)*diag(t) + alpha2*matrix(1,t,t)
  invR0 <- solve(R0)
  X0 <- cbind(diag(t), matrix(0,t,1))
  gmu0 <- c(X0%*%beta)
  mu0 <- plogis(gmu0)
  W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
  
  bm1 <- (1-alpha0+alpha1-alpha2)*diag(t*K)
  bm2 <- (alpha2-alpha1)*kronecker(matrix(1,t,t), diag(K))
  bm3 <- (alpha0-alpha1)*kronecker(diag(t), matrix(1,K,K))
  bm4 <- alpha1*matrix(1,t*K,t*K)
  R1 <- bm1+bm2+bm3+bm4
  invR1 <- solve(R1)
  X1 <- kronecker(cbind(diag(t), matrix(1,t,1)), matrix(1,K,1))
  gmu1 <- c(X1%*%beta)
  mu1 <- plogis(gmu1)
  W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
  
  I <- 0
  t_power <- 0
  while (t_power<(1-gamma)){
    I <- I + (K+1)
    N <- 2*K*I/(K+1)
    Omega <- (N/2)*(t(X0)%*%W0%*%X0) + (N/20)*(t(X1)%*%W1%*%X1)
    vardelta <- solve(Omega)[(t+1),(t+1)]
    t_power <- pt(qt(0.05/2, df=I-2)+abs(beta[t+1])/sqrt(vardelta), df=I-2)
  }
  
  return(data.frame(b=beta[t+1], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, a_power=t_power))
}

# Analytical power
bj_in <- 0.1*c(0.5^1, 0.5^2, 0.5^3, 0.5^4)
bj <- cumsum(bj_in)
beta3 <- c(bj[1:3], log(2))
beta4 <- c(bj, log(2))
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.03, 0.015, 0.2))) # N=200, I=110
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.1, 0.05, 0.2))) # N=260, I=143
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.01, 0.005, 0.4))) # N=220, I=121
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.03, 0.015, 0.2))) # N=180, I=99
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.1, 0.05, 0.2))) # N=240, I=132
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.01, 0.005, 0.4))) # N=200, I=110



###################################################
######### (4) Linear time effect #########
# with interaction between time and intervention
# global test
###################################################
binSIMULATE <- function(K,t,beta,alpha,gamma=0.15){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  alpha2 <- alpha[3]
  
  R0 <- (1-alpha2)*diag(t) + alpha2*matrix(1,t,t)
  invR0 <- solve(R0)
  X0 <- cbind(matrix(1,t,1), seq(1,t,1), matrix(0,t,1), matrix(0,t,1))
  gmu0 <- c(X0%*%beta)
  mu0 <- plogis(gmu0)
  W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
  
  bm1 <- (1-alpha0+alpha1-alpha2)*diag(t*K)
  bm2 <- (alpha2-alpha1)*kronecker(matrix(1,t,t), diag(K))
  bm3 <- (alpha0-alpha1)*kronecker(diag(t), matrix(1,K,K))
  bm4 <- alpha1*matrix(1,t*K,t*K)
  R1 <- bm1+bm2+bm3+bm4
  invR1 <- solve(R1)
  X1 <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(1,t,1), seq(1,t,1)), matrix(1,K,1))
  gmu1 <- c(X1%*%beta)
  mu1 <- plogis(gmu1)
  W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
  
  I <- 0
  F_power <- 0
  while (F_power<(1-gamma)){
    I <- I + (K+1)
    N <- 2*K*I/(K+1)
    Omega <- (N/2)*(t(X0)%*%W0%*%X0) + (N/20)*(t(X1)%*%W1%*%X1)
    vardelta <- solve(Omega)[3:4, 3:4]
    lambda <- (beta[3:4])%*%(solve(vardelta))%*%(beta[3:4])
    fc <- qf(1-0.05, df1=2, df2=I-3)
    F_power <- 1-pf(fc, df1=2, df2=I-3, ncp=lambda)
  }
  
  return(data.frame(b=beta[4], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, a_power=F_power))
}

# Analytical power
beta3 <- c(log(0.2), log(1.005), log(1.1), log(1.4))
beta4 <- c(log(0.2), log(1.005), log(1.1), log(1.3))
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.03, 0.015, 0.2))) # N=200, I=110
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.1, 0.05, 0.2))) # N=240, I=132
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.01, 0.005, 0.4))) # N=200, I=110
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.03, 0.015, 0.2))) # N=180, I=99
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.1, 0.05, 0.2))) # N=220, I=121
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.01, 0.005, 0.4))) # N=180, I=99



###################################################
######### (5) Categorical time effect #########
# witht interaction between time and intervention
###################################################
binSIMULATE <- function(K,t,beta,alpha,gamma=0.15){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  alpha2 <- alpha[3]
  
  R0 <- (1-alpha2)*diag(t) + alpha2*matrix(1,t,t)
  invR0 <- solve(R0)
  X0 <- cbind(diag(t), matrix(0,t,t))
  gmu0 <- c(X0%*%beta)
  mu0 <- plogis(gmu0)
  W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
  
  bm1 <- (1-alpha0+alpha1-alpha2)*diag(t*K)
  bm2 <- (alpha2-alpha1)*kronecker(matrix(1,t,t), diag(K))
  bm3 <- (alpha0-alpha1)*kronecker(diag(t), matrix(1,K,K))
  bm4 <- alpha1*matrix(1,t*K,t*K)
  R1 <- bm1+bm2+bm3+bm4
  invR1 <- solve(R1)
  X1 <- kronecker(cbind(diag(t), diag(t)), matrix(1,K,1))
  gmu1 <- c(X1%*%beta)
  mu1 <- plogis(gmu1)
  W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
  
  I <- 0
  F_power <- 0
  while (F_power<(1-gamma)){
    I <- I + (K+1)
    N <- 2*K*I/(K+1)
    Omega <- (N/2)*(t(X0)%*%W0%*%X0) + (N/20)*(t(X1)%*%W1%*%X1)
    vardelta <- solve(Omega)[(t+1):(2*t), (t+1):(2*t)]
    lambda <- (beta[(t+1):(2*t)])%*%(solve(vardelta))%*%(beta[(t+1):(2*t)])
    fc <- qf(1-0.05, df1=t, df2=I-t-1)
    F_power <- 1-pf(fc, df1=t, df2=I-t-1, ncp=lambda)
  }
  
  return(data.frame(b=beta[2*t], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    N=N, I=I, K=K, t=t, a_power=F_power))
}

# Analytical power
b1j_in <- 0.1*c(0.5^1, 0.5^2, 0.5^3, 0.5^4)
b1j <- cumsum(b1j_in)
b2j_in <- 0.65*c(0.6^1, 0.6^2, 0.6^3, 0.6^4)
b2j <- cumsum(b2j_in)
beta3 <- c(b1j[1:3], b2j[1:3])
beta4 <- c(b1j, b2j)
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.03, 0.015, 0.2))) # N=320, I=176
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.1, 0.05, 0.2))) # N=440, I=242
aPower <- rbind(aPower, binSIMULATE(10, 3, beta3, c(0.01, 0.005, 0.4))) # N=340, I=187
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.03, 0.015, 0.2))) # N=260, I=143
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.1, 0.05, 0.2))) # N=340, I=187
aPower <- rbind(aPower, binSIMULATE(10, 4, beta4, c(0.01, 0.005, 0.4))) # N=280, I=154

write.xlsx(aPower, file = "binResults/pred_power_bin.xlsx", row.names = FALSE)


