##################################################
# Collection of results from raw data
# Size and Power
#
# NOTE: This file is not needed if using the
# intermediate results (i.e. the XLSX files)
# currently available in the folder "conResults".
##################################################
library(openxlsx)

conRES<-function(model,row,case,ind){
  path<-"conResults/conRData/"
  if(case==1){
    name<-"Size"
  } else{
    name<-"Power"
  }
  
  if(ind==1){
    namei<-"ind"
  } else{
    namei<-"MAEE"
  }
  
  load(paste0(path,model,"model/",row,"row/",name,"_beta_",namei,".RData"))
  
  if (model==4 | model==5) {
    # non-convergence
    err<-(results[,11]==1)
    nerror<-sum(err)
    errSIM<-results[err,1]
    
    # the correct results
    I<-as.numeric(results[1,9])
    p<-as.numeric(results[1,10])
    
    # remove non-convergence
    res1<-as.data.frame(results[!err,c(1:7)])
  } else {
    # non-convergence
    err<-(results[,10]==1)
    nerror<-sum(err)
    errSIM<-results[err,1]
    
    # the correct results
    I<-as.numeric(results[1,8])
    if (model==3){
      p<-2
    } else {
      p<-as.numeric(results[1,9])
    }
    
    # remove non-convergence
    if (ind==1){
      res1<-as.data.frame(results[!err,c(1:4,11:12)])
    } else{
      res1<-as.data.frame(results[!err,c(1:6)])
    }
    
    # calculate the average of KC and MD
    res1$AVG <- rowMeans(res1[,5:6], na.rm=TRUE)
  }
  
  # calculate power/size
  if (model==4){
    pval1<-colMeans(apply(res1[,3:7]/2,2,function(x){(1-pf(x,df1=2,df2=I-3))})<0.05)
  } else if (model==5){
    t<-as.numeric(results[1,2])
    pval1<-colMeans(apply(res1[,3:7]/t,2,function(x){(1-pf(x,df1=t,df2=I-t-1))})<0.05)
  } else{
    mcsd1<-sd(res1[,2])
    essd1<-colMeans(res1[,3:7])
    pval1<-colMeans(apply(abs(res1[,2])/res1[,3:7],2,function(x){2*(1-pt(x,df=I-p))})<0.05)
  }
  
  return(data.frame(MB=pval1[[1]], ROB=pval1[[2]], KC=pval1[[3]], MD=pval1[[4]], AVG=pval1[[5]], nerror=nerror))
}

# Collect data
source("conScenarios_Analytical.R")
aPower<-read.xlsx("conResults/pred_power_con.xlsx", colNames=TRUE)

for (ind in 1:2){
  sim_size<-NULL
  sim_power<-NULL
  for (model in c(1:5)){
    for (row in 1:6){
      # Size
      sim_size<-rbind(sim_size,conRES(model=model,row=row,case=1,ind=ind))
      
      # Power
      sim_power<-rbind(sim_power,conRES(model=model,row=row,case=2,ind=ind))
      
      if(ind==1){
        namei<-"ind"
      } else{
        namei<-"MAEE"
      }
    }
  }
  final_con<-cbind(aPower, sim_power, sim_size)
  write.xlsx(final_con, file = paste("conResults/final_con_",namei,".xlsx", sep=""), row.names = FALSE)
}


