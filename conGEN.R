#############################################################
# Simulate correlated continuous outcomes
# in a lomgitudinal IRGT

# Dec 2021

# INPUT
# N: number of individuals
# I: number of clusters
# K: number of individuals per cluster in the intervention arm
# t: number of time points
# phi: dispersion parameter
# beta: vector of regression parameters
# alpha: correlations
#       alpha[1] - within-period correlation
#       alpha[2] - inter-period correlation
#       alpha[3] - within-individual correlation
# model: time effect assumption of mean model
#        1 - no time effect
#        21 - categorical time effect without interaction between time and intervention
#        22 - categorical time effect with interaction between time and intervention
#        31 - linear time effect without interaction between time and intervention
#        32 - linear time effect with interaction between time and intervention
#############################################################

conGEN <- function(N,I,K,t,phi,beta,alpha,model){
  require(mvtnorm)
  
  ########################################################
  # Create block exchangeable correlation matrix
  # for the intervention arm
  ########################################################
  
  bxch <- function(alpha){
    alpha0 <- alpha[1]
    alpha1 <- alpha[2]
    alpha2 <- alpha[3]
    bm1 <- (1 - alpha0 + alpha1 - alpha2) * diag(t*K)
    bm2 <- (alpha2 - alpha1) * kronecker(matrix(1,t,t), diag(K))
    bm3 <- (alpha0 - alpha1) * kronecker(diag(t), matrix(1,K,K))
    bm4 <- alpha1 * matrix(1,t*K,t*K)
    r <- bm1+bm2+bm3+bm4
    return(r)
  }
  
  ########################################################
  # Create simple exchangeable correlation matrix
  # for the control arm
  ########################################################
  
  sxch <- function(alpha){
    alpha2 <- alpha[3]
    sm1 <- (1 - alpha2) * diag(t)
    sm2 <- alpha2 * matrix(1,t,t)
    r <- sm1+sm2
    return(r)
  }
  
  ########################################################
  # returns variance matrix of Gaussian variables
  # with dispersion parameter phi and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,phi){
    return(phi*r)
  }
  
  # Create X matrix
  if (model==1){ # no time effect
    x_c <- cbind(matrix(1,t,1), matrix(0,t,1))
    x_t <- cbind(matrix(1,t*K,1), matrix(1,t*K,1))
  } else if (model==2){ # linear time effect without interaction between time and intervention
    x_c <- cbind(matrix(1,t,1), seq(1,t,1), matrix(0,t,1))
    x_t <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(1,t,1)), matrix(1,K,1))
  } else if (model==3){ # categorical time effect without interaction between time and intervention
    x_c <- cbind(diag(t), matrix(0,t,1))
    x_t <- kronecker(cbind(diag(t), matrix(1,t,1)), matrix(1,K,1))
  } else if (model==4){ # linear time effect with interaction between time and intervention
    x_c <- cbind(matrix(1,t,1), seq(1,t,1), matrix(0,t,1), matrix(0,t,1))
    x_t <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(1,t,1), seq(1,t,1)), matrix(1,K,1))
  } else if (model==5){ # categorical time effect with interaction between time and intervention
    x_c <- cbind(diag(t), matrix(0,t,t))
    x_t <- kronecker(cbind(diag(t), diag(t)), matrix(1,K,1))
  }
  
  # Simulate correlated Gaussian outcomes for the control arm
  y <- NULL
  I_c <- N/2
  u <- x_c%*%beta
  r <- sxch(alpha)
  v <- cor2var(r,phi)
  y <- c(y, c(t(rmvnorm(I_c,u,v))))
  
  # Simulate correlated Gaussian outcomes for the intervention arm
  I_t <- N/(2*K)
  u <- x_t%*%beta
  r <- bxch(alpha)
  v <- cor2var(r,phi)
  y <- c(y, c(t(rmvnorm(I_t,u,v))))
  
  # Return simulated data vector
  return(y)
}


