##########################################################################
# Simulate correlated binary outcomes

# Ref: Qaqish, B. F. (2003). A family of multivariate binary distributions 
# for simulating correlated binary variables. Biometrika 90, 455-463.

# INPUT
# N: number of individuals
# I: number of clusters
# K: number of individuals per cluster in the intervention arm
# t: number of time points
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
##########################################################################

binGEN<-function(N,I,K,t,beta,alpha,model){
  
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
  # a[1:n, 1:n] is the input covariance matrix of Y[1:n].
  # Returns  b[1:n,1:n] such that b[1:t-1, t] are the 
  # slopes for regression of y[t] on y[1:t-1], for t=2:n.
  # Diagonals and lower half of b[,] are copied from a[,].
  # a[,] is assumed +ve definite symmetric, not checked.
  ########################################################
  
  allreg<-function(a){
    n<-nrow(a)
    b<-a
    for(t in 2:n){
      t1<-t-1
      gt<-a[1:t1,1:t1]
      st<-a[1:t1,t]
      bt<-solve(gt,st)
      b[1:t1,t]<-bt
    }
    return(b)
  }
  
  ########################################################
  # returns variance matrix of binary variables with mean
  # vector u[] and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,u){
    s<-diag(sqrt(u*(1-u)))
    return(s%*%r%*%s)
  }
  
  ########################################################
  # r[1:n, 1:n] is the corr mtx
  # u[1:n] is the mean of a binary vector
  # checks that pairwise corrs are in-range for the given u[]
  # only upper half of r[,] is checked.
  # return 0 if ok
  # return 1 if out of range
  ########################################################
  
  chkbinc<-function(r,u){
    n<-length(u)
    s<-sqrt(u*(1-u))
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        uij<-u[i]*u[j]+r[i,j]*s[i]*s[j]
        ok<-((uij <= min(u[i], u[j])) & (uij >= max(0, u[i]+u[j]-1)))
        if(!ok) {return(1)}
      }
    }
    return(0)
  }
  
  ########################################################
  # Multivariate Binary Simulation by Linear Regression.
  # Simulate a single vector.
  # Returns a simulated binary random vector y[1:n] with mean 
  # u[1:n] and regression coefs matrix b[1:n,1:n] (obtained 
  # by calling allreg() above).
  # y[] and u[] are column vectors.
  # Returns -1 if the cond. linear family not reproducible
  ########################################################
  
  mbslr1<-function(b,u){
    n<-nrow(b)
    y<-rep(-1,n)
    y[1]<-rbinom(1,1,u[1])
    for(i in 2:n){
      i1<-i-1
      r<-y[1:i1]-u[1:i1]              # residuals
      ci<-u[i]+sum(r*b[1:i1,i])       # cond.mean
      if(ci < 0 | ci > 1){
        stop(paste("mbslr1: ERROR:",ci))
        return(-1)
      }
      y[i]<-rbinom(1,1,ci)
    }
    return(y)
  }
  
  # Create X matrix and mu
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
  
  # Simulate correlated binary outcomes
  y<-NULL
  
  for (i in 1:(N/2)){
    gmu_c <- c(x_c%*%beta)
    u <- plogis(gmu_c)
    r <- sxch(alpha)
    v <- cor2var(r,u)                   # v <- cov matrix
    oor <- chkbinc(r,u)                 # corr out of range?
    if(oor){
      stop("ERROR: Corr out of range for given mean")
    }
    b <- allreg(v)     # prepare coeffs
    y <- c(y,mbslr1(b,u))         # simulate data matrix
  }
  
  for (i in (N/2+1):I){
    gmu_t <- c(x_t%*%beta)
    u <- plogis(gmu_t)
    r <- bxch(alpha)
    v <- cor2var(r,u)                   # v <- cov matrix
    oor <- chkbinc(r,u)                 # corr out of range?
    if(oor){
      stop("ERROR: Corr out of range for given mean")
    }
    b <- allreg(v)     # prepare coeffs
    y <- c(y,mbslr1(b,u))         # simulate data matrix
  }
  
  # Return simulated data matrix
  return(y)
}

