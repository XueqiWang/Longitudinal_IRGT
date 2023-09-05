########################################################################################################################
# Function of sample size calculations for longitudinal IRGT trials with a continuous or binary outcome

# August 2023

# INPUT
# Kt: number of individuals in the treatment arm
# Kc: number of individuals in the control arm
# t: number of time points
# family: "gaussian", or "binomial"
# link: "identity", "logit", or "log"
#       (default: link = "identity" for family = "gaussian";
#                 link = "logit" for family = "binomial")
# model: "NT-CTE", "LT-CTE", "CT-CTE", "LT-TI", or "CT-TI"
#       ("NT-CTE" - no-time constant treatment effect mode
#        "LT-CTE" - linear-time constant treatment effect model
#        "CT-CTE" - Categorical-time constant treatment effect mode
#        "LT-TI" - linear time by treatment interaction mode
#        "CT-TI" - categorical time by treatment interaction model
#        default: model = "NT-CTE")
# b: effect size if family="gaussian", or vector of b if family="binomial"
# phit: dispersion parameter for the treatment arm
# phic: dispersion parameter for the control arm
# alpha - correlations in the treatment arm
#         alpha[1] - within-period correlation in the treatment arm
#         alpha[2] - inter-period correlation in the treatment arm
#         alpha[3] - within-individual correlation in the treatment arm
# rho - correlations in the control arm
#         rho[1] - within-period correlation in the control arm
#         rho[2] - inter-period correlation in the control arm
#         rho[3] - within-individual correlation in the control arm
# qc: proportion of clusters assigned to the control arm
# size: Type I error rate
# gamma: type II error rate; power = 1-gamma
########################################################################################################################

lirgt_ss <- function(Kt, Kc, t, family="gaussian", link=NULL, model = "NT-CTE", b, phit=1, phic=1, alpha, rho, qc=1/2, size=0.05, gamma=0.15){
  
  require(MASS)
  
  # link
  if (is.null(link)){
    if (family=="gaussian"){
      link <- "identity"
    } else if (family=="binomial"){
      link <- "logit"
    }
  }
  
  # elements
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  alpha2 <- alpha[3]
  rho0 <- rho[1]
  rho1 <- rho[2]
  rho2 <- rho[3]
  Im <- as.numeric(strsplit(attr(fractions(qc),"fracs"), "/")[[1]][2])
  
  # continuous outcome
  if (family=="gaussian"){
    
    # eigenvalues of the block exchangeable correlation matrix for the treatment arm
    lambda_t3 <- 1+(Kt-1)*(alpha0-alpha1)-alpha2
    lambda_t4 <- 1+(Kt-1)*alpha0+(t-1)*(Kt-1)*alpha1+(t-1)*alpha2
    
    # eigenvalues of the block exchangeable correlation matrix for the control arm
    lambda_c3 <- 1+(Kc-1)*(rho0-rho1)-rho2
    lambda_c4 <- 1+(Kc-1)*rho0+(t-1)*(Kc-1)*rho1+(t-1)*rho2
    
    # required number of clusters
    if (model=="NT-CTE"){    # model  1
      sigma2 <- 1/t*(phic*lambda_c4/(qc*Kc)+phit*lambda_t4/((1-qc)*Kt))
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        pred_power <- pt(qt(size/2, df=I-2)+abs(b)*sqrt(I/sigma2), df=I-2)
      }
    }
    
    else if (model=="LT-CTE"){    # model  2
      sigma2 <- 1/t*(phic*lambda_c4/(qc*Kc)+phit*lambda_t4/((1-qc)*Kt))
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        pred_power <- pt(qt(size/2, df=I-3)+abs(b)*sqrt(I/sigma2), df=I-3)
      }
    }
    
    else if (model=="CT-CTE"){    # model  3
      sigma2 <- 1/t*(phic*lambda_c4/(qc*Kc)+phit*lambda_t4/((1-qc)*Kt))
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        pred_power <- pt(qt(size/2, df=I-2)+abs(b)*sqrt(I/sigma2), df=I-2)
      }
    }
    
    else if (model=="LT-TI"){    # model  4
      V3 <- (phic*lambda_c3/(qc*Kc)+phit*lambda_t3/((1-qc)*Kt))*12/(t*(t^2-1))
      V4 <- (phic*lambda_c4/(qc*Kc)+phit*lambda_t4/((1-qc)*Kt))/t
      Va <- V3*(t+1)^2/4 + V4
      Vb <- -V3*(t+1)/2
      Vd <- V3
      Vsub_inv <- solve(matrix(c(Va,Vb,Vb,Vd),2,2))
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        lambda <- I*(b%*%Vsub_inv%*%b)
        fc <- qf(1-size, df1=2, df2=I-3)
        pred_power <- 1-pf(fc, df1=2, df2=I-3, ncp=lambda)
      }
    }
    
    else if (model=="CT-TI"){    # model  5
      Va <- phic*lambda_c3/(qc*Kc)+phit*lambda_t3/((1-qc)*Kt)
      Vb <- (phic*(lambda_c4-lambda_c3)/(qc*Kc)+phit*(lambda_t4-lambda_t3)/((1-qc)*Kt))/t
      Vc <- 1/Va
      Vd <- -Vb/(Va*(Va+t*Vb))
      Vsub_inv <- Vc*diag(t) + Vd*matrix(1,t,t)
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        lambda <- I*(b%*%Vsub_inv%*%b)
        fc <- qf(1-size, df1=t, df2=I-t-1)
        pred_power <- 1-pf(fc, df1=t, df2=I-t-1, ncp=lambda)
      }
    }
  }
  
  # binary outcome
  else if (family=="binomial"){
    
    # block exchangeable correlation matrix for the treatment arm
    bmt1 <- (1-alpha0+alpha1-alpha2)*diag(t*Kt)
    bmt2 <- (alpha2-alpha1)*kronecker(matrix(1,t,t), diag(Kt))
    bmt3 <- (alpha0-alpha1)*kronecker(diag(t), matrix(1,Kt,Kt))
    bmt4 <- alpha1*matrix(1,t*Kt,t*Kt)
    R1 <- bmt1+bmt2+bmt3+bmt4
    invR1 <- solve(R1)
    
    # block exchangeable correlation matrix for the control arm
    bmc1 <- (1-rho0+rho1-rho2)*diag(t*Kc)
    bmc2 <- (rho2-rho1)*kronecker(matrix(1,t,t), diag(Kc))
    bmc3 <- (rho0-rho1)*kronecker(diag(t), matrix(1,Kc,Kc))
    bmc4 <- rho1*matrix(1,t*Kc,t*Kc)
    R0 <- bmc1+bmc2+bmc3+bmc4
    invR0 <- solve(R0)
    
    # sample size calculation
    if (model=="NT-CTE"){    # model  1
      X0 <- cbind(matrix(1,t*Kc,1), matrix(0,t*Kc,1))
      X1 <- cbind(matrix(1,t*Kt,1), matrix(1,t*Kt,1))
      
      if (link=="logit"){    # binomial with logit link
        gmu0 <- c(X0%*%b)
        mu0 <- plogis(gmu0)
        W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- plogis(gmu1)
        W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
      } else if (link=="identity"){    # binomial with identity link
        mu0 <- c(X0%*%b)
        W0 <- diag(sqrt(1/(mu0*(1-mu0))))%*%invR0%*%diag(sqrt(1/(mu0*(1-mu0))))
        mu1 <- c(X1%*%b)
        W1 <- diag(sqrt(1/(mu1*(1-mu1))))%*%invR1%*%diag(sqrt(1/(mu1*(1-mu1))))
      } else if (link=="log"){    # binomial with log link
        gmu0 <- c(X0%*%b)
        mu0 <- exp(gmu0)
        W0 <- diag(sqrt(mu0/(1-mu0)))%*%invR0%*%diag(sqrt(mu0/(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- exp(gmu1)
        W1 <- diag(sqrt(mu1/(1-mu1)))%*%invR1%*%diag(sqrt(mu1/(1-mu1)))
      }
      
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        Omega <- I*qc*(t(X0)%*%W0%*%X0) + I*(1-qc)*(t(X1)%*%W1%*%X1)
        vardelta <- solve(Omega)[2,2]
        pred_power <- pt(qt(size/2, df=I-2)+abs(b[2])/sqrt(vardelta), df=I-2)
      }
    }
    
    else if (model=="LT-CTE"){    # model  2
      X0 <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(0,t,1)), matrix(1,Kc,1))
      X1 <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(1,t,1)), matrix(1,Kt,1))
      
      if (link=="logit"){    # binomial with logit link
        gmu0 <- c(X0%*%b)
        mu0 <- plogis(gmu0)
        W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- plogis(gmu1)
        W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
      } else if (link=="identity"){    # binomial with identity link
        mu0 <- c(X0%*%b)
        W0 <- diag(sqrt(1/(mu0*(1-mu0))))%*%invR0%*%diag(sqrt(1/(mu0*(1-mu0))))
        mu1 <- c(X1%*%b)
        W1 <- diag(sqrt(1/(mu1*(1-mu1))))%*%invR1%*%diag(sqrt(1/(mu1*(1-mu1))))
      } else if (link=="log"){    # binomial with log link
        gmu0 <- c(X0%*%b)
        mu0 <- exp(gmu0)
        W0 <- diag(sqrt(mu0/(1-mu0)))%*%invR0%*%diag(sqrt(mu0/(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- exp(gmu1)
        W1 <- diag(sqrt(mu1/(1-mu1)))%*%invR1%*%diag(sqrt(mu1/(1-mu1)))
      }
      
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        Omega <- I*qc*(t(X0)%*%W0%*%X0) + I*(1-qc)*(t(X1)%*%W1%*%X1)
        vardelta <- solve(Omega)[3,3]
        pred_power <- pt(qt(size/2, df=I-2)+abs(b[3])/sqrt(vardelta), df=I-2)
      }
    }
    
    else if (model=="CT-CTE"){    # model  3
      X0 <- kronecker(cbind(diag(t), matrix(0,t,1)), matrix(1,Kc,1))
      X1 <- kronecker(cbind(diag(t), matrix(1,t,1)), matrix(1,Kt,1))
      
      if (link=="logit"){    # binomial with logit link
        gmu0 <- c(X0%*%b)
        mu0 <- plogis(gmu0)
        W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- plogis(gmu1)
        W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
      } else if (link=="identity"){    # binomial with identity link
        mu0 <- c(X0%*%b)
        W0 <- diag(sqrt(1/(mu0*(1-mu0))))%*%invR0%*%diag(sqrt(1/(mu0*(1-mu0))))
        mu1 <- c(X1%*%b)
        W1 <- diag(sqrt(1/(mu1*(1-mu1))))%*%invR1%*%diag(sqrt(1/(mu1*(1-mu1))))
      } else if (link=="log"){    # binomial with log link
        gmu0 <- c(X0%*%b)
        mu0 <- exp(gmu0)
        W0 <- diag(sqrt(mu0/(1-mu0)))%*%invR0%*%diag(sqrt(mu0/(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- exp(gmu1)
        W1 <- diag(sqrt(mu1/(1-mu1)))%*%invR1%*%diag(sqrt(mu1/(1-mu1)))
      }
      
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        Omega <- I*qc*(t(X0)%*%W0%*%X0) + I*(1-qc)*(t(X1)%*%W1%*%X1)
        vardelta <- solve(Omega)[(t+1),(t+1)]
        pred_power <- pt(qt(size/2, df=I-2)+abs(b[t+1])/sqrt(vardelta), df=I-2)
      }
    }
    
    else if (model=="LT-TI"){    # model  4
      X0 <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(0,t,1), matrix(0,t,1)), matrix(1,Kc,1))
      X1 <- kronecker(cbind(matrix(1,t,1), seq(1,t,1), matrix(1,t,1), seq(1,t,1)), matrix(1,Kt,1))
      
      if (link=="logit"){    # binomial with logit link
        gmu0 <- c(X0%*%b)
        mu0 <- plogis(gmu0)
        W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- plogis(gmu1)
        W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
      } else if (link=="identity"){    # binomial with identity link
        mu0 <- c(X0%*%b)
        W0 <- diag(sqrt(1/(mu0*(1-mu0))))%*%invR0%*%diag(sqrt(1/(mu0*(1-mu0))))
        mu1 <- c(X1%*%b)
        W1 <- diag(sqrt(1/(mu1*(1-mu1))))%*%invR1%*%diag(sqrt(1/(mu1*(1-mu1))))
      } else if (link=="log"){    # binomial with log link
        gmu0 <- c(X0%*%b)
        mu0 <- exp(gmu0)
        W0 <- diag(sqrt(mu0/(1-mu0)))%*%invR0%*%diag(sqrt(mu0/(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- exp(gmu1)
        W1 <- diag(sqrt(mu1/(1-mu1)))%*%invR1%*%diag(sqrt(mu1/(1-mu1)))
      }
      
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        Omega <- I*qc*(t(X0)%*%W0%*%X0) + I*(1-qc)*(t(X1)%*%W1%*%X1)
        vardelta <- solve(Omega)[3:4, 3:4]
        lambda <- (b[3:4])%*%(solve(vardelta))%*%(b[3:4])
        fc <- qf(1-size, df1=2, df2=I-3)
        pred_power <- 1-pf(fc, df1=2, df2=I-3, ncp=lambda)
      }
    }
    
    else if (model=="CT-TI"){    # model  5
      X0 <- kronecker(cbind(diag(t), matrix(0,t,t)), matrix(1,Kc,1))
      X1 <- kronecker(cbind(diag(t), diag(t)), matrix(1,Kt,1))
      
      if (link=="logit"){    # binomial with logit link
        gmu0 <- c(X0%*%b)
        mu0 <- plogis(gmu0)
        W0 <- diag(sqrt(mu0*(1-mu0)))%*%invR0%*%diag(sqrt(mu0*(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- plogis(gmu1)
        W1 <- diag(sqrt(mu1*(1-mu1)))%*%invR1%*%diag(sqrt(mu1*(1-mu1)))
      } else if (link=="identity"){    # binomial with identity link
        mu0 <- c(X0%*%b)
        W0 <- diag(sqrt(1/(mu0*(1-mu0))))%*%invR0%*%diag(sqrt(1/(mu0*(1-mu0))))
        mu1 <- c(X1%*%b)
        W1 <- diag(sqrt(1/(mu1*(1-mu1))))%*%invR1%*%diag(sqrt(1/(mu1*(1-mu1))))
      } else if (link=="log"){    # binomial with log link
        gmu0 <- c(X0%*%b)
        mu0 <- exp(gmu0)
        W0 <- diag(sqrt(mu0/(1-mu0)))%*%invR0%*%diag(sqrt(mu0/(1-mu0)))
        gmu1 <- c(X1%*%b)
        mu1 <- exp(gmu1)
        W1 <- diag(sqrt(mu1/(1-mu1)))%*%invR1%*%diag(sqrt(mu1/(1-mu1)))
      }
      
      I <- pred_power <- 0
      while (pred_power<(1-gamma)){
        I <- I + Im
        Omega <- I*qc*(t(X0)%*%W0%*%X0) + I*(1-qc)*(t(X1)%*%W1%*%X1)
        vardelta <- solve(Omega)[(t+1):(2*t), (t+1):(2*t)]
        lambda <- (b[(t+1):(2*t)])%*%(solve(vardelta))%*%(b[(t+1):(2*t)])
        fc <- qf(1-size, df1=t, df2=I-t-1)
        pred_power <- 1-pf(fc, df1=t, df2=I-t-1, ncp=lambda)
      }
    }
  }
  
  # required number of individuals
  N <- I*(qc*Kc+(1-qc)*Kt)

  # final results
  return(data.frame(N=N, I=I, qc=qc, Kt=Kt, Kc=Kc, t=t, pred_power=pred_power))
  
}


# ---------------------------------------------------------------------------------

# examples for a continuous variable
Kt <- 10; Kc <- 1; t <- 3; phit <- 1; phic <- 1; qc <- 10/11
alpha <- c(0.03, 0.015, 0.2); rho <- c(0, 0, 0.2)

# NT-CTE model
b <- 0.35
lirgt_ss(Kt, Kc, t, family="gaussian", link=NULL, model="NT-CTE", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# LT-CTE model
b <- 0.35
lirgt_ss(Kt, Kc, t, family="gaussian", link=NULL, model="LT-CTE", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# CT-CTE model
b <- 0.35
lirgt_ss(Kt, Kc, t, family="gaussian", link=NULL, model="CT-CTE", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# LT-TI model
b <- c(0.18, 0.1)
lirgt_ss(Kt, Kc, t, family="gaussian", link=NULL, model="LT-TI", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# CT-TI model
b <- cumsum(0.3*c(0.6^1, 0.6^2, 0.6^3))
lirgt_ss(Kt, Kc, t, family="gaussian", link=NULL, model="CT-TI", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)


# examples for a binary variable
Kt <- 10; Kc <- 1; t <- 3; phit <- 1; phic <- 1; qc <- 10/11
alpha <- c(0.03, 0.015, 0.2); rho <- c(0, 0, 0.2)

# NT-CTE model
b <- c(0, log(2))
lirgt_ss(Kt, Kc, t, family="binomial", link=NULL, model="NT-CTE", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# LT-CTE model
b <- c(log(0.7), log(1.2), log(1.9))
lirgt_ss(Kt, Kc, t, family="binomial", link=NULL, model="LT-CTE", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# CT-CTE model
b <- c(cumsum(0.1*c(0.5^1, 0.5^2, 0.5^3)), log(2))
lirgt_ss(Kt, Kc, t, family="binomial", link=NULL, model="CT-CTE", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# LT-TI model
b <- c(log(0.2), log(1.005), log(1.1), log(1.4))
lirgt_ss(Kt, Kc, t, family="binomial", link=NULL, model="LT-TI", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)

# CT-TI model
b <- c(cumsum(0.1*c(0.5^1, 0.5^2, 0.5^3)), cumsum(0.65*c(0.6^1, 0.6^2, 0.6^3)))
lirgt_ss(Kt, Kc, t, family="binomial", link=NULL, model="CT-TI", b=b, phit=1, phic=1, alpha, rho, qc=qc, size=0.05, gamma=0.15)


