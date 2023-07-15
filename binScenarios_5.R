############################################
# Simulation program for binary outcomes
# Longitudinal IRGTs
# Categorical time effect -
# with interaction between time and intervention
############################################
# N: number of individuals
# I: number of clusters
# K: number of individuals per cluster in the intervention arm
# t: number of time points
# beta: vector of regression parameters
# alpha: correlations
#       alpha[1] - within-period correlation
#       alpha[2] - inter-period correlation
#       alpha[3] - within-individual correlation
############################################

###############################
# GEE with working independence
###############################

binSIMULATE<-function(N,I,K,t,beta,alpha,model){
  
  # Source program code
  source("binGEN.R")
  source("geesmv_new.R")
  
  require(pracma)
  require(gee)
  require(geesmv)
  
  # Store results
  nrep <- 1000
  nerror <- 0
  p <- length(beta)
  results <- matrix(NA,nrep,79)
  colnames(results) <- c("nSIM","T","MB","ROB","KC","MD","AVG","N","I","p","error",
                         rep("beta", 4), rep("MB_cov", 16), rep("ROB_cov", 16), rep("KC_cov", 16), rep("MD_cov", 16))
  results[,8] <- N
  results[,9] <- I
  results[,10] <- p
  results[,11] <- 0
  results[,2] <- t
  
  # Create X matrix
  I_c <- N/2
  I_t <- N/(2*K)
  x_c <- cbind(diag(t), matrix(0,t,t))
  x_t <- kronecker(cbind(diag(t), diag(t)), matrix(1,K,1))
  X <- rbind(repmat(x_c, I_c, 1), repmat(x_t, I_t, 1))
  
  # Create id
  id <- c(rep(1:I_c, each=t), rep((I_c+1):(I_c+I_t), each=t*K))
  
  # For Loop
  j <- 0
  for(i in 1:nrep){
    # Generate outcome and check; store number of simulated data sets
    passGEN <- 1
    nSIM <- 0
    while (passGEN==1){
      nSIM <- nSIM+1
      j <- j+1
      set.seed(j)
      y <- try(binGEN(N,I,K,t,beta,alpha,model),silent=TRUE)
      if (class(y)[1]=="try-error"){passGEN <- 1; next}
      passGEN <- 0
    }
    results[i,1]<-nSIM
    
    # GEE program
    mydata <- data.frame(cbind(y,X,id))
    outfit <- try(gee(y~-1+X,id=id,data=mydata,family=binomial(link="logit"),corstr="independence"),silent=TRUE)
    if (class(outfit)[1]=="try-error" | is.null(outfit)==1){
      nerror <- nerror+1
      results[i,11] <- 1
      next
    }
    
    coef_2j <- as.numeric(outfit$coefficients)[(p/2+1):p]
    naive_2j <- outfit$naive.variance[(p/2+1):p,(p/2+1):p]
    robust_2j <- outfit$robust.variance[(p/2+1):p,(p/2+1):p]
    
    results[i, 12:(12+t-1)]<-coef_2j
    
    var <- try(coef_2j%*%solve(naive_2j)%*%coef_2j,silent=TRUE)
    if(class(var)[1]=="try-error" | is.null(var)==1){
      nerror <- nerror+1
      results[i,11] <- 1
    } else {
      results[i,3] <- var
      results[i, 16:(16+t^2-1)] <- as.numeric(naive_2j)
    }
    
    var <- try(coef_2j%*%solve(robust_2j)%*%coef_2j,silent=TRUE)
    if(class(var)[1]=="try-error" | is.null(var)==1){
      nerror <- nerror+1
      results[i,11] <- 1
    } else {
      results[i,4] <- var
      results[i, 32:(32+t^2-1)] <- as.numeric(robust_2j)
    }
    
    outfit <- try(GEE.var.kc.new(y~-1+X,id=id,data=mydata,family=binomial,corstr="independence"),silent=TRUE)
    if (class(outfit)[1]=="try-error" | is.null(outfit)==1){
      nerror <- nerror+1
      results[i,11] <- 1
      kc.new_2j <- NA
    } else {
      kc.new_2j <- outfit$cov.beta[(p/2+1):p,(p/2+1):p]
      var <- try(coef_2j%*%solve(kc.new_2j)%*%coef_2j,silent=TRUE)
      if(class(var)[1]=="try-error" | is.null(var)==1){
        nerror <- nerror+1
        results[i,11] <- 1
      } else {
        results[i,5] <- var
        results[i, 48:(48+t^2-1)] <- as.numeric(kc.new_2j)
      }
    }
    
    outfit <- try(GEE.var.md.new(y~-1+X,id=id,data=mydata,family=binomial,corstr="independence"),silent=TRUE)
    if (class(outfit)[1]=="try-error" | is.null(outfit)==1){
      nerror <- nerror+1
      results[i,11] <- 1
      md.new_2j <- NA
    } else {
      md.new_2j <- outfit$cov.beta[(p/2+1):p,(p/2+1):p]
      var <- try(coef_2j%*%solve(md.new_2j)%*%coef_2j,silent=TRUE)
      if(class(var)[1]=="try-error" | is.null(var)==1){
        nerror <- nerror+1
        results[i,11] <- 1
      } else {
        results[i,6] <- var
        results[i, 64:(64+t^2-1)] <- as.numeric(md.new_2j)
      }
    }
    
    if ((is.na(kc.new_2j) + is.na(md.new_2j))[1] == 0){
      avg_2j <- (kc.new_2j+md.new_2j)/2
      diag(avg_2j) <- (sqrt(diag(kc.new_2j)) + sqrt(diag(md.new_2j)))^2/4
      var <- try(coef_2j%*%solve(avg_2j)%*%coef_2j,silent=TRUE)
      if(class(var)[1]=="try-error" | is.null(var)==1){
        nerror <- nerror+1
        results[i,11] <- 1
      } else {
        results[i,7] <- var
      }
    } else {
      results[i,11] <- 1
    }
    
    # Loop control
    if(i %% 20 == 0){print(paste("Iteration",i,"NERROR =",nerror))}
  }

  # save data sets and clear workspace
  if (beta[p]==0){
    name<-"Size"
  } else{
    name<-"Power"}
  save(results,file=paste0(name,"_beta_ind.RData"))
}

# ---------------------------------------------------------------------------------

# Power

b1j_in <- 0.1*c(0.5^1, 0.5^2, 0.5^3, 0.5^4)
b1j <- cumsum(b1j_in)
b2j_in <- 0.65*c(0.6^1, 0.6^2, 0.6^3, 0.6^4)
b2j <- cumsum(b2j_in)
beta3 <- c(b1j[1:3], b2j[1:3])
beta4 <- c(b1j, b2j)

# (1)
setwd('binResults/binRData/5model/1row')
binSIMULATE(320, 176, 10, 3, beta3, c(0.03, 0.015, 0.2), 5) # N=320, I=176

# (2)
setwd('../2row')
binSIMULATE(440, 242, 10, 3, beta3, c(0.1, 0.05, 0.2), 5) # N=440, I=242

# (3)
setwd('../3row')
binSIMULATE(340, 187, 10, 3, beta3, c(0.01, 0.005, 0.4), 5) # N=340, I=187

# (4)
setwd('../4row')
binSIMULATE(260, 143, 10, 4, beta4, c(0.03, 0.015, 0.2), 5) # N=260, I=143

# (5)
setwd('../5row')
binSIMULATE(340, 187, 10, 4, beta4, c(0.1, 0.05, 0.2), 5) # N=340, I=187

# (6)
setwd('../6row')
binSIMULATE(280, 154, 10, 4, beta4, c(0.01, 0.005, 0.4), 5) # N=280, I=154

# Size

b1j_in <- 0.1*c(0.5^1, 0.5^2, 0.5^3, 0.5^4)
b1j <- cumsum(b1j_in)
b2j <- rep(0, 4)
beta3 <- c(b1j[1:3], b2j[1:3])
beta4 <- c(b1j, b2j)

# (1)
setwd('../1row')
binSIMULATE(320, 176, 10, 3, beta3, c(0.03, 0.015, 0.2), 5) # N=320, I=176

# (2)
setwd('../2row')
binSIMULATE(440, 242, 10, 3, beta3, c(0.1, 0.05, 0.2), 5) # N=440, I=242

# (3)
setwd('../3row')
binSIMULATE(340, 187, 10, 3, beta3, c(0.01, 0.005, 0.4), 5) # N=340, I=187

# (4)
setwd('../4row')
binSIMULATE(260, 143, 10, 4, beta4, c(0.03, 0.015, 0.2), 5) # N=260, I=143

# (5)
setwd('../5row')
binSIMULATE(340, 187, 10, 4, beta4, c(0.1, 0.05, 0.2), 5) # N=340, I=187

# (6)
setwd('../6row')
binSIMULATE(280, 154, 10, 4, beta4, c(0.01, 0.005, 0.4), 5) # N=280, I=154



###############################
# GEE with MAEE
###############################

binSIMULATE<-function(N,I,K,t,beta,alpha,model){
  
  # Source program code
  source("binGEN.R")
  source("binMAEE_F.R")
  
  require(pracma)
  require(gee)
  
  # Store results
  nrep <- 1000
  nerror <- 0
  p <- length(beta)
  results <- matrix(NA,nrep,79)
  colnames(results) <- c("nSIM","T","MB","ROB","KC","MD","AVG","N","I","p","error",
                         rep("beta", 4), rep("MB_cov", 16), rep("ROB_cov", 16), rep("KC_cov", 16), rep("MD_cov", 16))
  results[,8] <- N
  results[,9] <- I
  results[,10] <- p
  results[,11] <- 0
  results[,2] <- t
  
  results0 <- matrix(NA,nrep,5)
  results1 <- matrix(NA,nrep,5)
  results2 <- matrix(NA,nrep,5)
  colnames(results0) <- c("alpha0","ROB","KC","MD","error")
  colnames(results1) <- c("alpha1","ROB","KC","MD","error")
  colnames(results2) <- c("alpha2","ROB","KC","MD","error")
  
  # Create X matrix
  I_c <- N/2
  I_t <- N/(2*K)
  x_c <- cbind(diag(t), matrix(0,t,t))
  x_t <- kronecker(cbind(diag(t), diag(t)), matrix(1,K,1))
  X <- rbind(repmat(x_c, I_c, 1), repmat(x_t, I_t, 1))
  
  # Create Z matrix
  CREATEZ<-function(t,K){
    alpha0_pos <- 1
    alpha1_pos <- 2
    alpha2_pos <- 3
    
    sm1 <- (1 - alpha2_pos) * diag(t)
    sm2 <- alpha2_pos * matrix(1,t,t)
    POS_c <- sm1+sm2
    
    bm1 <- (1 - alpha0_pos + alpha1_pos - alpha2_pos) * diag(t*K)
    bm2 <- (alpha2_pos - alpha1_pos) * kronecker(matrix(1,t,t), diag(K))
    bm3 <- (alpha0_pos - alpha1_pos) * kronecker(diag(t), matrix(1,K,K))
    bm4 <- alpha1_pos * matrix(1,t*K,t*K)
    POS_t <- bm1+bm2+bm3+bm4
    
    zrow<-diag(3)
    z_c<-NULL
    z_t<-NULL
    for(jp in 1:(t-1)){
      for(tp in (jp+1):(t)){
        z_c<-rbind(z_c,zrow[POS_c[jp,tp],])
      }
    }
    Z<-z_c
    for(i in 1:(I_c-1)){
      Z<-rbind(Z,z_c)
    }
    
    for(jp in 1:(t*K-1)){
      for(tp in (jp+1):(t*K)){
        z_t<-rbind(z_t,zrow[POS_t[jp,tp],])
      }
    }
    for(i in (I_c+1):(I_c+I_t)){
      Z<-rbind(Z,z_t)
    }
    
    rm(z_c)
    rm(z_t)
    return(Z)
  }
  Z<-CREATEZ(t,K) # large matrix
  
  # Create id
  clsize <- c(rep(t, I_c), rep(t*K, I_t))
  id <- c(rep(1:I_c, each=t), rep((I_c+1):(I_c+I_t), each=t*K))
  
  # For Loop
  j <- 0
  for(i in 1:nrep){
    # Generate outcome and check; store number of simulated data sets
    passGEN <- 1
    nSIM <- 0
    while (passGEN==1){
      nSIM <- nSIM+1
      j <- j+1
      set.seed(j)
      y <- try(binGEN(N,I,K,t,beta,alpha,model),silent=TRUE)
      if (class(y)[1]=="try-error"){passGEN <- 1; next}
      passGEN <- 0
    }
    results[i,1]<-nSIM
    
    # MAEE program
    out=try(binMAEE(y=y,X=X,id=id,n=clsize,Z=Z,maxiter=25,epsilon=0.001,printrange="NO",shrink="ALPHA",makevone="NO"),silent=TRUE)
    if(class(out)[1]=="try-error" | is.null(out)==1){
      nerror<-nerror+1
      results[i,11]<-1
      next
    }
    
    # Save results
    coef_2j <- out$outbeta[(p/2+1):p,1]
    
    results[i, 12:(12+t-1)]<-coef_2j
    
    var <- try(coef_2j%*%solve(out$outbVAR[(p/2+1):p,(p/2+1):p])%*%coef_2j,silent=TRUE)
    if(class(var)[1]=="try-error" | is.null(var)==1){
      nerror <- nerror+1
      results[i,11] <- 1
    } else {
      results[i,3] <- var
      results[i, 16:(16+t^2-1)] <- as.numeric(out$outbVAR[(p/2+1):p,(p/2+1):p])
    }
    
    var <- try(coef_2j%*%solve(out$outbVARBC0[(p/2+1):p,(p/2+1):p])%*%coef_2j,silent=TRUE)
    if(class(var)[1]=="try-error" | is.null(var)==1){
      nerror <- nerror+1
      results[i,11] <- 1
    } else {
      results[i,4] <- var
      results[i, 32:(32+t^2-1)] <- as.numeric(out$outbVARBC0[(p/2+1):p,(p/2+1):p])
    }
    
    var <- try(coef_2j%*%solve(out$outbVARBC1[(p/2+1):p,(p/2+1):p])%*%coef_2j,silent=TRUE)
    if(class(var)[1]=="try-error" | is.null(var)==1){
      nerror <- nerror+1
      results[i,11] <- 1
    } else {
      results[i,5] <- var
      results[i, 48:(48+t^2-1)] <- as.numeric(out$outbVARBC1[(p/2+1):p,(p/2+1):p])
    }
    
    var <- try(coef_2j%*%solve(out$outbVARBC2[(p/2+1):p,(p/2+1):p])%*%coef_2j,silent=TRUE)
    if(class(var)[1]=="try-error" | is.null(var)==1){
      nerror <- nerror+1
      results[i,11] <- 1
    } else {
      results[i,6] <- var
      results[i, 64:(64+t^2-1)] <- as.numeric(out$outbVARBC2[(p/2+1):p,(p/2+1):p])
    }
    
    VARAVG <- (out$outbVARBC1[(p/2+1):p,(p/2+1):p]+out$outbVARBC2[(p/2+1):p,(p/2+1):p])/2
    diag(VARAVG) <- (sqrt(diag(out$outbVARBC1[(p/2+1):p,(p/2+1):p])) + sqrt(diag(out$outbVARBC2[(p/2+1):p,(p/2+1):p])))^2/4
    var <- try(coef_2j%*%solve(VARAVG)%*%coef_2j,silent=TRUE)
    if(class(var)[1]=="try-error" | is.null(var)==1){
      nerror <- nerror+1
      results[i,11] <- 1
    } else {
      results[i,7] <- var
    }
    
    results0[i,(1:4)] <- out$outalpha[1,1:4]
    results1[i,(1:4)] <- out$outalpha[2,1:4]
    results2[i,(1:4)] <- out$outalpha[3,1:4]
    
    # Loop control
    if(i %% 20 == 0){print(paste("Iteration",i,"NERROR =",nerror))}
  }
  
  results0[,5]<-results[,11]
  results1[,5]<-results[,11]
  results2[,5]<-results[,11]
  
  # save data sets and clear workspace
  if (beta[p]==0){
    name<-"Size"
  } else{
    name<-"Power"}
  save(results,file=paste0(name,"_beta_MAEE.RData"))
  save(results0,file=paste0(name,"_alpha0_MAEE.RData"))
  save(results1,file=paste0(name,"_alpha1_MAEE.RData"))
  save(results2,file=paste0(name,"_alpha2_MAEE.RData"))
}

# ---------------------------------------------------------------------------------

b1j_in <- 0.1*c(0.5^1, 0.5^2, 0.5^3, 0.5^4)
b1j <- cumsum(b1j_in)
b2j_in <- 0.65*c(0.6^1, 0.6^2, 0.6^3, 0.6^4)
b2j <- cumsum(b2j_in)
beta3 <- c(b1j[1:3], b2j[1:3])
beta4 <- c(b1j, b2j)

# (1)
setwd('../1row')
binSIMULATE(320, 176, 10, 3, beta3, c(0.03, 0.015, 0.2), 5) # N=320, I=176

# (2)
setwd('../2row')
binSIMULATE(440, 242, 10, 3, beta3, c(0.1, 0.05, 0.2), 5) # N=440, I=242

# (3)
setwd('../3row')
binSIMULATE(340, 187, 10, 3, beta3, c(0.01, 0.005, 0.4), 5) # N=340, I=187

# (4)
setwd('../4row')
binSIMULATE(260, 143, 10, 4, beta4, c(0.03, 0.015, 0.2), 5) # N=260, I=143

# (5)
setwd('../5row')
binSIMULATE(340, 187, 10, 4, beta4, c(0.1, 0.05, 0.2), 5) # N=340, I=187

# (6)
setwd('../6row')
binSIMULATE(280, 154, 10, 4, beta4, c(0.01, 0.005, 0.4), 5) # N=280, I=154

# Size

b1j_in <- 0.1*c(0.5^1, 0.5^2, 0.5^3, 0.5^4)
b1j <- cumsum(b1j_in)
b2j <- rep(0, 4)
beta3 <- c(b1j[1:3], b2j[1:3])
beta4 <- c(b1j, b2j)

# (1)
setwd('../1row')
binSIMULATE(320, 176, 10, 3, beta3, c(0.03, 0.015, 0.2), 5) # N=320, I=176

# (2)
setwd('../2row')
binSIMULATE(440, 242, 10, 3, beta3, c(0.1, 0.05, 0.2), 5) # N=440, I=242

# (3)
setwd('../3row')
binSIMULATE(340, 187, 10, 3, beta3, c(0.01, 0.005, 0.4), 5) # N=340, I=187

# (4)
setwd('../4row')
binSIMULATE(260, 143, 10, 4, beta4, c(0.03, 0.015, 0.2), 5) # N=260, I=143

# (5)
setwd('../5row')
binSIMULATE(340, 187, 10, 4, beta4, c(0.1, 0.05, 0.2), 5) # N=340, I=187

# (6)
setwd('../6row')
binSIMULATE(280, 154, 10, 4, beta4, c(0.01, 0.005, 0.4), 5) # N=280, I=154

setwd('../../../..')


