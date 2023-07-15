############################################################################################
# GEE analysis of cluster randomized trials with binary outcomes
# Adapted from the "PRENTICEcorrJSP" SAS macro of Preisser.

# The following estimating equations are used for estimating correlation parameters:
# MAEE: Matrix-adjusted estimating equations of Preisser et al., 2008

# Ref. Preisser, J. S., Lu, B., and Qaqish, B. F. (2008). Finite sample adjustments 
# in estimating equations and covariance estimators for intracluster correlations. 
# Statistics in Medicine 27, 5764-5785.

# The following variance estimators are output from the program:
# (1) MB: Model-based estimator (only available for marginal mean parameters)
# (2) BC0: Uncorrected sandwich estimator that extends Liang and Zeger, 1986
# (3) BC1: Bias-corrected sandwich estimator that extends Kauermann and Carroll, 2001
# (4) BC2: Bias-corrected sandwich estimator that extends Mancl and DeRouen, 2001
# (5) BC3: Bias-corrected sandwich estimator that extends Fay and Graubard, 2001
# (6) BC4: Bias-corrected sandwich estimator that extends Morel, Bokossa, and Neerchal, 2003

# INPUT
# y: The binary outcome variable
# X: Marginal mean covariates (design matrix including intercept)
# id: Cluster identifier
# n: Vector of cluster sample sizes
# Z: Marginal correlation covariates (design matrix for correlation)
# maxiter: Maximum number of iterations for Fisher Scoring
# epsilon: Tolerance for convergence
# printrange: Print details of range violations? (Should initially check) (YES, NO)
#             Range is always checked, however. Range is always checked and printed
#             when variance computed.       
# shrink: THETA or ALPHA. Used to tune step sizes in case of non-convergence.
# makevone: Sets VEE = 1 if YES instead of w_ijk (if NO) from Appendix B. 
#           See also Ref. Sharples and Breslow (1992).

# Note that all inputs are required.
# ID's should be integers ranging from 1 to the number of clusters. Data should be sorted by
# ID before calling contMAEE.  Z matrix should be all pairs j < k for each cluster i, ie
# (1,2) (1,3).... (1,n)... (2,3).... (n-1, n) for each cluster.
############################################################################################

binMAEE=function(y, X, id, n, Z, maxiter, epsilon, 
                 printrange, shrink, makevone){
  require(MASS)
  
  #####################################################################################
  # MODULE: BEGINEND
  # Creates two vectors that have the start and end points for each cluster
  
  # INPUT
  # n: Vector of cluster sample sizes
  
  # OUTPUT
  # first: Vector with starting row for cluster i
  # last: Vector with ending row for cluster i
  #####################################################################################
  
  BEGINEND=function(n){
    last=cumsum(n)
    first=last-n+1
    return(cbind(first,last))
  }
  
  #####################################################################################
  # Module: IS_POS_DEF
  # A = symmetric matrix
  # returns 1 if A is positive definite
  # 0 otherwise
  #####################################################################################
  
  is_pos_def=function(A){
    return(min(eigen(A)$values)>1e-13)
  }
  
  #####################################################################################
  # MODULE: GETROWB
  # Generates a row of the -E(d(corr)/dbeta) (ie Q) matrix
  
  # INPUT
  # mu: marginal means for cluster i
  # gamma: pairwise correlation for outcome j and k for cluster i
  # j: indicator for mean 
  # k: indicator for mean 
  # X: covariate matrix for cluster i
  # y: response vector for cluster i
  
  # OUTPUT
  # Row of -E(d(corr)/dbeta) (ie Q) matrix
  #####################################################################################
  
  GETROWB=function(mu,gamma,j,k,X,y){
    row=(gamma/2)*((1-2*mu[j])*X[j,]+(1-2*mu[k])*X[k,]+
                     2*mu[j]*(1-mu[j])*X[j,]/(y[j]-mu[j])+
                     2*mu[k]*(1-mu[k])*X[k,]/(y[k]-mu[k]))
    return(row)
  }
  
  #####################################################################################
  # MODULE: CREATEA
  # Creates residual for beta estimating equation, (Y - mu)
  
  # INPUT
  # mu: Vector of n_i marginal means
  # y: Outcome vector for ith cluster
  
  # OUTPUT
  # residuals for beta estimating equation
  #####################################################################################
  
  CREATEA=function(mu,y){
    return(y-mu)
  }
  
  #####################################################################################
  # MODULE: CREATEC
  # Creates derivative matrix for beta estimating equation, dmu/dbeta
  
  # INPUT
  # X: Covariate matrix for cluster i
  # mu: Vector of n_i marginal means
  
  # OUTPUT
  # derivative matrix for beta estimating equation
  #####################################################################################
  
  CREATEC=function(X,mu){
    return(X*(mu*(1-mu)))
  }
  
  #####################################################################################
  # MODULE: CREATEB
  # Creates covariance matrix for beta estimating equation, var(Y)
  
  # INPUT
  # mu: Vector of n_i marginal means
  # n: Sample size (scalar) for cluster i
  # gamma: Vector of rho_jk between outcome j and k for cluster i
  
  # OUTPUT
  # covariance matrix for beta estimating equation
  #####################################################################################
  
  CREATEB=function(mu,gamma,n){
    B=diag(mu*(1-mu))
    l=1
    for(j in 1:(n-1)){
      for(k in (j+1):n){
        B[j,k]<-sqrt(mu[j]*mu[k]*(1-mu[j])*(1-mu[k]))*gamma[l]
        l=l+1
      }
    }
    B[lower.tri(B)] = t(B)[lower.tri(B)]
    return(B)
  }
  
  #####################################################################################
  # MODULE: SCORE
  # Generates the score matrix for each cluster and approximate information
  # to be used to estimate parameters and generate standard errors
  
  # INPUT
  # Ustarold: initial values for information matrix
  # beta: Vector of marginal mean parameters
  # alpha: Vector of marginal correlation parameters
  # y: The 0/1 binary outcome
  # X: Marginal mean covariates
  # Z: Marginal correlation covariates
  # n: Vector of cluster sample sizes
  # p: Number of marginal mean parameters
  # q: Number of marginal correlation parameters
  # flag: Performs an Eigenanalysis of B to see if positive definite. Only
  #       called when computing the variance at the end (0 = no, 1 = yes). Prints
  #       warning for each cluster violation.
  # rangeflag: Checks to see if correlation is within range based on marginal means
  #           (0 = in range, 1 = out of range). See Prentice (1988).
  # VEEFLAG:  Checks if all variances for alpha e.e. are positive, terminates if not
  # NPSDFLAG: Checks if B is positive definite, terminates if not
  # NPSDADJFLAG: Checks if the alpha adjustment factor matrix is positive definite, terminates if not
  
  # OUTPUT
  # U: Score vector
  # UUtran: Sum of U_i*U_i` across all clusters
  # Ustar: Approximate information matrix
  #####################################################################################
  
  SCORE=function(Ustarold, beta, alpha, y, X, Z, n, p,
                 q,flag, rangeflag, VEEFLAG, NPSDFLAG, NPSDADJFLAG){
    U=rep(0,p+q)
    UUtran=Ustar=matrix(0,p+q,p+q)
    naiveold=ginv(Ustarold[1:p,1:p]) # needed for Hi1 below
    
    locx=BEGINEND(n)
    locz=BEGINEND(choose(n,2))
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      Z_c = Z[locz[i,1]:locz[i,2],]
      
      U_c=rep(0,p+q)
      Ustar_c=matrix(0,p+q,p+q)
      mu_c=1/(1+exp(c(-X_c%*%beta)))
      gamma_c=c(Z_c%*%alpha)
      
      VEE=R=rep(0,choose(n[i],2))
      DB=matrix(0,choose(n[i],2),p)
      
      C=CREATEC(X_c,mu_c)
      B=CREATEB(mu_c,gamma_c,n[i])
      A=CREATEA(mu_c,y_c)
      
      INVB=ginv(B)
      CtinvB=t(C)%*%INVB
      Hi1=C%*%naiveold%*%CtinvB
      
      SQVARFUN=sqrt(mu_c*(1-mu_c))
      INVSQVAR=1/SQVARFUN
      CT=t(C)
      omega=C%*%naiveold%*%CT
      vminomega=B-omega
      psd_vmin=is_pos_def(vminomega)
      mineig=min(eigen(vminomega)$values)
      if(psd_vmin==1){
        # if(mineig<1e-12){cat("vminomega with psd_min =",psd_min,"mininum eigenvalue is",mineig,"\n")}
        # Ci=B%*%ginv(vminomega)
        Ci=diag(INVSQVAR)%*%(B%*%ginv(vminomega))%*%diag(SQVARFUN)
        RX=(y_c-mu_c)*INVSQVAR
        Gi=tcrossprod(RX)
      } else{
        NPSDADJFLAG=1
        stop("vminomega is not positive definite")
      }
      
      # RANGE CHECKS
      l=1
      for(j in 1:(n[i]-1)){
        for(k in (j+1):n[i]){
          if((gamma_c[l] >= min(sqrt((mu_c[j]*(1-mu_c[k]))/(mu_c[k]*(1-mu_c[j]))),
                                sqrt((mu_c[k]*(1-mu_c[j]))/(mu_c[j]*(1-mu_c[k]))))) | 
             (gamma_c[l] <= max(-sqrt((mu_c[j]*mu_c[k])/((1-mu_c[j])*(1-mu_c[k]))),
                                -sqrt(((1-mu_c[j])*(1-mu_c[k]))/(mu_c[j]*mu_c[k])))) 
             & (flag==0)){
            rangeflag=1
            if(printrange=="YES"){
              cat("Range Violation Detected for Cluster",i,"and Pair",j,k,"\n")
            }
            break
          }
          if((gamma_c[l] >= min(sqrt((mu_c[j]*(1-mu_c[k]))/(mu_c[k]*(1-mu_c[j]))),
                                sqrt((mu_c[k]*(1-mu_c[j]))/(mu_c[j]*(1-mu_c[k]))))) | 
             (gamma_c[l] <= max(-sqrt((mu_c[j]*mu_c[k])/((1-mu_c[j])*(1-mu_c[k]))),
                                -sqrt(((1-mu_c[j])*(1-mu_c[k]))/(mu_c[j]*mu_c[k])))) 
             & (flag==1)){
            cat("Last Update Pushes Parameters Out of Range.","\n")
            cat("Range Violation Detected for Cluster",i,"and Pair",j,k,"\n")
          }
          VEE[l]<-1+((1-2*mu_c[j])*(1-2*mu_c[k])*gamma_c[l]) / 
            sqrt(mu_c[j]*(1-mu_c[j])*mu_c[k]*(1-mu_c[k]))-gamma_c[l]^2
          
          # insert check that variance is nonnegative
          if(VEE[l]<=0){
            VEEFLAG=1
            stop("Variance of correlation parameter is negative")
          }
          
          # DB[l,]<-GETROWB(mu_c,gamma_c[l],j,k,X_c,y_c)
          
          # Matrix-based multiplicative correction, (I - H_i)^{-1}
          R[l]<-Ci[j, ]%*%Gi[ ,k]-gamma_c[l]  
          
          l=l+1
        }
      }
      
      if(makevone=="YES"){
        VEE=rep(1,choose(n[i],2))
      }
      
      # Check for pos def of B;
      # For the simulations only, the check is done at every iteration;
      # score module is exited if not p.s.d., otherwise whole program would crash;
      if(min(eigen(B)$values)<=0){
        NPSDFLAG=1
        stop(paste("Var(Y) of Cluster",i,"is not Positive-Definite;","Joint Distribution Does Not Exist and Program terminates"))
      }
      
      U_c[1:p]=t(C)%*%INVB%*%A
      U_c[(p+1):(p+q)]=t(Z_c)%*%(R/VEE)  # this is specific to the linear structure
      UUtran_c=tcrossprod(U_c)
      Ustar_c[1:p,1:p]=t(C)%*%INVB%*%C
      
      # Next line is commented out to give updating equations (17) in Prentice (1988);
      # A more detailed procedure is given by un-commenting to get procedure based upon (15) of Prentice;
      # Ustar_c[(p+1):(p+q),1:p]=t(Z_c)%*%(DB/VEE)
      Ustar_c[(p+1):(p+q),(p+1):(p+q)]=t(Z_c)%*%(Z_c/VEE)
      
      U=U+U_c
      UUtran=UUtran+UUtran_c
      Ustar=Ustar+Ustar_c
    }
    rangeflag=0
    return(list(U=U,UUtran=UUtran,Ustar=Ustar,flag=flag,rangeflag=rangeflag,VEEFLAG=VEEFLAG,
                NPSDFLAG=NPSDFLAG,NPSDADJFLAG=NPSDADJFLAG))
  }
  
  #####################################################################################
  # MODULE: INITBETA
  # Generates initial values for beta. Approximates logistic regression
  # using Newton's method
  
  # INPUT
  # y: The 0/1 binary outcomer4
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes
  
  # OUTPUT
  # beta: Vector of marginal mean parameters
  # Ustar: approximate score vector
  #####################################################################################
  
  INITBETA=function(y,X,n){
    z=y+y-1
    beta=solve(t(X)%*%X,t(X)%*%z)
    
    for(i in 1:2){
      u=c(X%*%beta)
      u=1/(1+exp(-u))
      v=u*(1-u)
      z=t(X)%*%(y-u)
      Ustar=t(X)%*%(X*v)
      d=solve(Ustar,z)
      beta=beta+d
    }
    return(list(beta=c(beta),Ustar=Ustar))
  }
  
  
  #####################################################################################
  # MODULE: INVBIG
  # compute (A - mm`)^{-1}c without performing the inverse directly
  # INPUT
  # ainvc: inverse of matrix A times vector c
  # ainvm: inverse of matrix A times matrix (with low no. of columns) M
  # M: matrix of eigen column vectors m1,m2, ..
  # c: vector 
  # start: of do loop
  # end:   of do loop, rank of X
  
  # OUTPUT
  # ainvc: inverse of matrix A times vector c
  #####################################################################################
  
  INVBIG=function(ainvc,ainvm,m,c,start,end){
    for(i in start:end){
      b=ainvm[,i]
      bt=t(b)
      btm=bt%*%m
      btmi=btm[,i]
      gam=1-btmi
      bg=b/gam
      ainvc=ainvc+bg%*%(bt%*%c)
      if(i<end){
        ainvm=ainvm+bg%*%btm
      }
    }
    return(ainvc)
  }
  
  #####################################################################################
  # MODULE: MAKEVAR
  # Creates covariance matrix of beta and alpha.
  
  # INPUT
  # beta: Vector of marginal mean parameters
  # alpha: Vector of marginal correlation parameters
  # Ustarold: initial values for information matrix
  # y: The 0/1 binary outcome
  # X: Marginal mean covariates
  # Z: Marginal correlation covariates
  # n: Vector of cluster sample sizes
  # p: Number of marginal mean parameters
  # q: Number of marginal correlation parameters
  # VEEFLAG:  Checks if all variances for alpha e.e. are positive, terminates if not
  # NPSDFLAG: Checks if B is positive definite, terminates if not
  # NPSDADJFLAG: Checks if the alpha adjustment factor matrix is positive definite, terminates if not
  
  # OUTPUT
  # robust: Robust covariance matrix for beta and alpha
  # naive: Naive (Model-Based) covariance matrix for beta
  # varMD: bias-corrected variance by Mancl and Derouen (2001)
  # varKC: bias-corrected variance by Kauermann and Carroll (2001)
  # varFG: bias-corrected variance by Fay and Graubard (2001)
  # varMBN: bias-corrected variance by Morel, Bokossa, and Neerchal (2003)
  #####################################################################################
  
  MAKEVAR=function(Ustarold, beta, alpha, y, X, Z, n, p, q,
                   VEEFLAG,ROBFLAG,NPSDADJFLAG){
    
    SCORE_RES=SCORE(Ustarold, beta, alpha, y, X, Z, n, p, q, 
                    flag=1, rangeflag=0,VEEFLAG,NPSDFLAG,NPSDADJFLAG)
    U=SCORE_RES$U
    UUtran=SCORE_RES$UUtran
    Ustar=SCORE_RES$Ustar
    flag=SCORE_RES$flag
    rangeflag=SCORE_RES$rangeflag
    VEEFLAG=SCORE_RES$VEEFLAG
    NPSDFLAG=SCORE_RES$NPSDFLAG
    NPSDADJFLAG=SCORE_RES$NPSDADJFLAG
    
    naive=ginv(Ustar[1:p,1:p])
    naivealp=ginv(Ustar[(p+1):(p+q),(p+1):(p+q)])
    
    # new commands to compute INV(I - H1)
    eigenRES1=eigen(naive)
    evals1=eigenRES1$values
    evecs1=eigenRES1$vectors
    sqrevals1=sqrt(evals1)
    sqe1=evecs1%*%diag(sqrevals1)
    
    # new commands to compute INV(I - H2)
    eigenRES2=eigen(naivealp)
    evals2=eigenRES2$values
    evecs2=eigenRES2$vectors
    sqrevals2=sqrt(evals2)
    sqe2=evecs2%*%diag(sqrevals2)
    
    # Bias-corrected variance
    Ustar_c_array=UUtran_c_array=array(0,c(p+q,p+q,length(n)))
    UUtran=UUbc=UUbc2=UUbc3=Ustar=inustar=matrix(0,p+q,p+q)
    
    locx=BEGINEND(n)
    locz=BEGINEND(choose(n,2))
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      mu_c=1/(1+exp(c(-X_c%*%beta)))
      
      U_i=U_c=rep(0,p+q)
      Ustar_c=matrix(0,p+q,p+q)
      Z_c=Z[locz[i,1]:locz[i,2],]
      gamma_c=c(Z_c%*%alpha)
      
      # commands for beta
      C=CREATEC(X_c,mu_c)
      B=CREATEB(mu_c,gamma_c,n[i])
      A=CREATEA(mu_c,y_c)
      INVB=ginv(B)
      U_i[1:p]=t(C)%*%INVB%*%A
      
      CtinvB=t(C)%*%INVB 
      Hi1=C%*%naive%*%CtinvB
      
      # commands for generalized inverse - beta
      ai1=INVB
      mm1=C%*%sqe1
      ai1A=ai1%*%A
      ai1m1=ai1%*%mm1
      ai1A=INVBIG(ai1A,ai1m1,mm1,A,1,p)
      U_c[1:p]=t(C)%*%ai1A
      
      # commands for alpha
      VEE=R=rep(0,choose(n[i],2))
      DB=matrix(0,choose(n[i],2),p)
      RX=rep(0,n[i])
      
      # MAEE
      SQVARFUN=sqrt(mu_c*(1-mu_c))
      INVSQVAR=1/SQVARFUN
      CT=t(C)
      omega=C%*%naive%*%CT
      vminomega=B-omega
      psd_vmin=is_pos_def(vminomega)
      mineig=min(eigen(vminomega)$values)
      if(psd_vmin==1){
        # if(mineig<1e-12){cat("vminomega with psd_min =",psd_min,"mininum eigenvalue is",mineig,"\n")}
        # Ci=B%*%ginv(vminomega)
        Ci=diag(INVSQVAR)%*%(B%*%ginv(vminomega))%*%diag(SQVARFUN)
        RX=(y_c-mu_c)*INVSQVAR
        Gi=tcrossprod(RX)
      } else{
        NPSDADJFLAG=1
        stop("vminomega is not positive definite")
      }
      
      l=1
      for(j in 1:(n[i]-1)){
        for(k in (j+1):n[i]){
          VEE[l]<-1 + ((1-2*mu_c[j])*(1-2*mu_c[k])*gamma_c[l]) / 
            sqrt(mu_c[j]*(1-mu_c[j])*mu_c[k]*(1-mu_c[k])) - gamma_c[l]^2
          
          DB[l,]<-GETROWB(mu_c, gamma_c[l], j, k, X_c, y_c)
          
          # Matrix-based multiplicative correction, (I - H_i)^{-1}
          R[l]<-Ci[j, ]%*%Gi[ ,k]-gamma_c[l]  
          
          l=l+1
        }
      }
      
      U_i[(p+1):(p+q)]=t(Z_c)%*%(R/VEE)
      mm2=Z_c%*%sqe2
      ai2R=R/VEE
      ai2m2=mm2/VEE
      ai2R=INVBIG(ai2R,ai2m2,mm2,R,1,q)
      U_c[(p+1):(p+q)]<-t(Z_c)%*%ai2R
      
      Ustar_c[1:p,1:p]=t(C)%*%INVB%*%C
      Ustar_c[(p+1):(p+q),1:p]=t(Z_c)%*%(DB/VEE)
      Ustar_c[(p+1):(p+q),(p+1):(p+q)]=t(Z_c)%*%(Z_c/VEE)
      Ustar=Ustar+Ustar_c
      
      UUtran_c=tcrossprod(U_i)
      UUtran=UUtran+UUtran_c
      UUbc_c=tcrossprod(U_c)
      UUbc=UUbc+UUbc_c
      UUbc_ic=tcrossprod(U_c,U_i)
      UUbc2=UUbc2+UUbc_ic
      
      Ustar_c_array[,,i]=Ustar_c
      UUtran_c_array[,,i]=UUtran_c
    }
    
    inustar[1:p,1:p]=ginv(Ustar[1:p,1:p])
    inustar[(p+1):(p+q),(p+1):(p+q)]=ginv(Ustar[(p+1):(p+q),(p+1):(p+q)])
    inustar[(p+1):(p+q),1:p]=-inustar[(p+1):(p+q),(p+1):(p+q)]%*%Ustar[(p+1):(p+q),1:p]%*%inustar[1:p,1:p]
    
    # the minus sign above is crucial, esp. for large correlation;    
    inustartr=t(inustar)
    
    # calculating adjustment factor for BC3
    for(i in 1:length(n)){      
      Hi=diag(1/sqrt(1-pmin(0.75,c(diag(Ustar_c_array[,,i]%*%inustar)))))
      UUbc3=UUbc3+Hi%*%UUtran_c_array[,,i]%*%Hi
    }
    
    # BC0 or usual Sandwich estimator of Prentice (1988);     
    robust=inustar%*%UUtran%*%inustartr
    
    # BC1 or Variance estimator that extends Kauermann and Carroll (2001);
    varKC=inustar%*%(UUbc2+t(UUbc2))%*%inustartr/2
    
    # BC2 or Variance estimator that extends Mancl and DeRouen (2001);
    varMD=inustar%*%UUbc%*%inustartr
    
    # BC3 or Variance estimator that extends Fay and Graubard (2001);
    varFG=inustar%*%UUbc3%*%inustartr
    
    # BC4 or Variance estimator that extends Morel, Bokossa, and Neerchal (2003);
    varMBN=matrix(0,p+q,p+q)
    varMBN[1:p,1:p]=(sum(n)-1)*length(n)/((sum(n)-p)*(length(n)-1))*robust[1:p,1:p]+
      (min(0.5, p/(length(n)-p))*max(1,sum(diag(inustar[1:p,1:p]%*%((sum(n)-1)*length(n)/((sum(n)-p)*(length(n)-1))*UUtran[1:p,1:p])))/p))*inustartr[1:p,1:p]
    f_alpha=sum(sapply(n, function(x) choose(x,2)))
    varMBN[(p+1):(p+q),(p+1):(p+q)]=(f_alpha-1)*length(n)/((f_alpha-q)*(length(n)-1))*robust[(p+1):(p+q),(p+1):(p+q)]+
      (min(0.5, q/(length(n)-q))*max(1,sum(diag(inustar[(p+1):(p+q),(p+1):(p+q)]%*%((f_alpha-1)*length(n)/((f_alpha-q)*(length(n)-1))*UUtran[(p+1):(p+q),(p+1):(p+q)])))/q))*inustartr[(p+1):(p+q),(p+1):(p+q)]
    
    naive=inustar[1:p,1:p]
    
    if(min(diag(robust))<=0){ROBFLAG = 1}
    if(min(diag(varMD))<=0){ROBFLAG = 1}
    if(min(diag(varKC))<=0){ROBFLAG = 1}
    if(min(diag(varFG))<=0){ROBFLAG = 1}
    if(min(diag(varMBN))<=0){ROBFLAG = 1}
    
    return(list(robust=robust,naive=naive,varMD=varMD,varKC=varKC,varFG=varFG,varMBN=varMBN,
                VEEFLAG=VEEFLAG,ROBFLAG=ROBFLAG,NPSDADJFLAG=NPSDADJFLAG))
  }
  
  #####################################################################################
  # MODULE: FITPRENTICE
  # Performs the Prentice (1988) Method
  
  # INPUT
  # y: The 0/1 binary outcome
  # X: Marginal mean covariates
  # Z: Marginal correlation covariates
  # n: Vector of cluster sample sizes
  # maxiter: Max number of iterations
  # epsilon: Tolerence for convergence
  # VEEFLAG: THE ALGORITHM terminated due to nonpositive variance in weights
  # SINGFLAG: THE ALGORITHM terminated due to singular MB covariance matrix
  # ROBFLAG: THE ALGORITHM terminated due to singular robust variance matrix
  # VEEFLAG:  Checks if all variances for alpha e.e. are positive, terminates if not
  # NPSDFLAG: Checks if B is positive definite, terminates if not
  # NPSDADJFLAG: Checks if the alpha adjustment factor matrix is positive definite, terminates if not
  
  # OUTPUT
  # beta: A p x 1 vector of marginal mean parameter estimates  
  # alpha: A q x 1 vector of marginal correlation parameter estimates  
  # robust: Robust covariance matrix for beta and alpha
  # naive: Naive (Model-Based) covariance matrix for beta
  # varMD: Bias-corrected variance due to Mancl and DeRouen (2001)
  # varKC: Bias-corrected variance due to Kauermann and Carroll (2001)
  # varFG: Bias-corrected variance due to Fay and Graubard (2001)
  # varMBN: bias-corrected variance by Morel, Bokossa, and Neerchal (2003)
  # niter: Number of iterations required for convergence
  # converge: Did the algorithm converge (0 = no, 1 = yes)
  #####################################################################################
  
  FITPRENTICE=function(y, X, Z, n, maxiter, epsilon, VEEFLAG, SINGFLAG, ROBFLAG,
                       ALPFLAG, NPSDFLAG, NPSDADJFLAG){
    p=ncol(X)
    q=ncol(Z)
    delta=rep(2*epsilon,p+q)
    max_modi=20
    converge=0
    
    rangeflag=0
    alpha=rep(0.01,q)
    INITRES=INITBETA(y,X,n)
    beta=INITRES$beta
    Ustar=INITRES$Ustar
    
    niter=1
    while((niter<=maxiter) & (max(abs(delta))>epsilon)){
      n_modi=0
      
      SINGFLAG=0
      ALPFLAG=0
      repeat{
        Ustarold=Ustar
        
        NPSDFLAG=0
        NPSDADJFLAG=0
        SCORE_RES=SCORE(Ustarold, beta, alpha, y, X, Z, n, p,
                        q,flag=0, rangeflag, VEEFLAG, NPSDFLAG, NPSDADJFLAG)
        U=SCORE_RES$U
        UUtran=SCORE_RES$UUtran
        Ustar=SCORE_RES$Ustar
        rangeflag=SCORE_RES$rangeflag
        VEEFLAG=SCORE_RES$VEEFLAG
        
        if (VEEFLAG==1){
          stop("VEEFLAG=1 indicates program terminated due to division by zero in variance")
        }
        if(rangeflag==1){
          if(shrink=="THETA"){
            if(niter==1){
              alpha=rep(0,q)
            } else{
              theta=theta-(0.5)^(n_modi+1)*delta
              beta=theta[1:p]
              alpha=theta[(p+1):(p+q)]
            }
          } else if(shrink=="ALPHA"){
            if(niter==1){
              alpha=rep(0,q)
            } else{
              alpha=0.95*alpha
            }
          }
          n_modi=n_modi+1
          if(printrange=="YES"){
            cat("Iteration",niter,"and Shrink Number",n_modi,"\n")
          }
        }
        if((n_modi>max_modi) | (rangeflag==0)){break}
      }
      
      if(n_modi>max_modi){
        if(printrange=="YES"){
          cat("n_modi too great, more than 20 shrinks") 
        }
        ALPFLAG=1
      }
      theta=c(beta,alpha)
      
      psdustar=is_pos_def(Ustar)
      mineig=min(eigen(Ustar)$values)
      if(psdustar==TRUE){
        delta=solve(Ustar,U)
        theta=theta+delta
        beta=theta[1:p]
        alpha=theta[(p+1):(p+q)]
        converge=(max(abs(delta))<=epsilon)
      } else{
        SINGFLAG=1
      }
      niter<-niter+1
    }
    
    Ustarold=Ustar
    if(converge==1){
      MAKEVAR_RES=MAKEVAR(Ustarold, beta, alpha, y, X, Z, n, p, q, VEEFLAG, ROBFLAG, NPSDADJFLAG)
      robust=MAKEVAR_RES$robust
      naive=MAKEVAR_RES$naive
      varMD=MAKEVAR_RES$varMD
      varKC=MAKEVAR_RES$varKC
      varFG=MAKEVAR_RES$varFG
      varMBN=MAKEVAR_RES$varMBN
      VEEFLAG=MAKEVAR_RES$VEEFLAG
      ROBFLAG=MAKEVAR_RES$ROBFLAG
      NPSDADJFLAG=MAKEVAR_RES$NPSDADJFLAG
    }
    return(list(beta=beta,alpha=alpha,robust=robust,naive=naive,varMD=varMD,
                varKC=varKC,varFG=varFG,varMBN=varMBN,niter=niter,converge=converge,
                VEEFLAG=VEEFLAG,SINGFLAG=SINGFLAG,ROBFLAG=ROBFLAG,ALPFLAG=ALPFLAG, 
                NPSDFLAG=NPSDFLAG,NPSDADJFLAG=NPSDADJFLAG))
  }
  
  #####################################################################################
  # MODULE: RESULTS
  # Creates printed output to screen of parameters and other information
  
  # INPUT
  # beta: Vector of marginal mean parameters
  # alpha: Vector of marginal correlation Parameters  
  # robust: Robust covariance matrix for beta and alpha
  # naive: Naive (Model-Based) covariance matrix for beta 
  # niter: Number of iterations until convergence
  # n: Vector of cluster sample sizes
  
  # OUTPUT
  # To Screen
  #####################################################################################
  RESULTS=function(beta, alpha, robust, naive, varMD, varKC, varFG, varMBN, niter, n){
    p=length(beta)
    q=length(alpha)
    K=length(n)
    df=K-p
    
    #cat("GEE for correlated binary data","\n",
    #    "Number of Clusters:",K,"\n",
    #    "Maximum Cluster Size:",max(n),"\n",
    #    "Minimum Cluster Size:",min(n),"\n",
    #    "Number of Iterations:",niter,"\n")
    
    beta_numbers=as.matrix(seq(1:p))-1
    bSE=sqrt(diag(naive))
    bSEBC0=sqrt(diag(robust[1:p,1:p]))
    bSEBC1=sqrt(diag(varKC[1:p,1:p]))
    bSEBC2=sqrt(diag(varMD[1:p,1:p]))
    bSEBC3=sqrt(diag(varFG[1:p,1:p]))
    bSEBC4=sqrt(diag(varMBN[1:p,1:p]))
    
    alpha_numbers=as.matrix(seq(1:q))-1
    aSEBC0=sqrt(diag(robust[(p+1):(p+q),(p+1):(p+q)]))
    aSEBC1=sqrt(diag(varKC[(p+1):(p+q),(p+1):(p+q)]))
    aSEBC2=sqrt(diag(varMD[(p+1):(p+q),(p+1):(p+q)]))
    aSEBC3=sqrt(diag(varFG[(p+1):(p+q),(p+1):(p+q)]))
    aSEBC4=sqrt(diag(varMBN[(p+1):(p+q),(p+1):(p+q)]))
    
    outbeta=cbind(beta_numbers,beta,bSE,bSEBC0,bSEBC1,bSEBC2,bSEBC3,bSEBC4)
    outalpha=cbind(alpha_numbers,alpha,aSEBC0,aSEBC1,aSEBC2,aSEBC3,aSEBC4)
    colnames(outbeta)<-c("Beta","Estimate","MB-stderr","BC0-stderr","BC1-stderr","BC2-stderr","BC3-stderr","BC4-stderr")
    colnames(outalpha)<-c("Alpha","Estimate","BC0-stderr","BC1-stderr","BC2-stderr","BC3-stderr","BC4-stderr")
    
    #cat("Results for marginal mean parameters \n")
    #print(outbeta)
    #cat("\n")
    
    #cat("Results for correlation parameters \n")
    #print(outalpha)
    #cat("\n")
    return(list(outbeta=outbeta[,-1],outalpha=outalpha[,-1]))
  }
  
  # reasons for non-results are identified and tallied
  VEEFLAG=0
  SINGFLAG=0
  CONVFLAG=0 # omiited in the arguments for prentice fit
  ROBFLAG=0
  ALPFLAG=0
  NPSDFLAG=0
  NPSDADJFLAG=0
  
  # Fit the Prentice Model
  PRENTICE_RES=FITPRENTICE(y, X, Z, n, maxiter, epsilon, VEEFLAG, SINGFLAG, ROBFLAG,
                           ALPFLAG, NPSDFLAG, NPSDADJFLAG)
  beta=PRENTICE_RES$beta
  alpha=PRENTICE_RES$alpha
  robust=PRENTICE_RES$robust
  naive=PRENTICE_RES$naive
  varMD=PRENTICE_RES$varMD
  varKC=PRENTICE_RES$varKC
  varFG=PRENTICE_RES$varFG
  varMBN=PRENTICE_RES$varMBN
  niter=PRENTICE_RES$niter
  converge=PRENTICE_RES$converge
  VEEFLAG=PRENTICE_RES$VEEFLAG
  SINGFLAG=PRENTICE_RES$SINGFLAG
  ROBFLAG=PRENTICE_RES$ROBFLAG
  ALPFLAG=PRENTICE_RES$ALPFLAG
  NPSDFLAG=PRENTICE_RES$NPSDFLAG
  NPSDADJFLAG=PRENTICE_RES$NPSDADJFLAG
  
  # Final Results
  if(converge==0 & VEEFLAG==0 & SINGFLAG==0 & NPSDFLAG==0 & NPSDADJFLAG==0){CONVFLAG=1}
  if(converge==1 & ROBFLAG==0){
    RESULTS(beta, alpha, robust, naive, varMD, varKC, varFG, varMBN, niter, n)
  }
}

