######################################################################################
### Purpose: Wrapper of the sampling function to sample from posterier distribution 
###          under Sarmanov or indipendent beta prior distributions
### Input:   data(y1,n1,y2,n2), hyperparemeters(a1,a2,b1,b2,rho), comparative measure,
###          model, number of smaples (nsam)
### Output:  samples
### Note:    the implement function: "inde.sampling.posterior" and "sar.sampling.posterior"
### Author:  Sheng Luo and Xiao Su
### Data:    7/13/2012
######################################################################################
sampling <- function(a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,y1=0,n1=0,y2=0,n2=0,
                     model="Sarmanov",measure=measure,nsam=10000) {
  if (model=="Independent")
    samples <- inde.sampling.posterior(a1,b1,a2,b2,n1,y1,n2,y2,measure,nsam)
  else {
    if (model=="Sarmanov") {
      cc <- sqrt(a1*a2*b1*b2)/sqrt((a1+b1+1)*(a2+b2+1))
      upper.bound <- cc/max(a1*b2, a2*b1)
      lower.bound <- -cc/max(a1*a2, b1*b2)
      rho.range <- c(lower.bound,upper.bound)
      names(rho.range) <- c("lower.bound","upper.bound")
      if (rho > upper.bound | rho < lower.bound)
        stop(paste("rho is out of bound: ",
                   lower.bound, upper.bound))
      samples <- sar.sampling.posterior(a1,a2,b1,b2,n1,y1,n2,y2,rho,measure,nsam)
    }
  }
  return(samples)
}

##########################################################################################
### Purpose: sample from posterier distribution under independent beta prior distribution
### Input:   data(y1,n1,y2,n2), hyperparemeters(a1,a2,b1,b2,rho), comparative measure,
###          model, number of smaples (nsam)
### Output:  samples
### Note:    This function is called by the wrapper function "sampling"
### Author:  Sheng Luo and Xiao Su
### Data:    7/13/2012
##########################################################################################
inde.sampling.posterior <- function(a1,b1,a2,b2,n1,y1,n2,y2,measure,nsam) {
  alpha1 <- a1+y1
  beta1 <- b1+n1-y1
  alpha2 <- a2+y2
  beta2 <- b2+n2-y2
  p1 <- rbeta(nsam, shape1=alpha1, shape2=beta1)
  p2 <- rbeta(nsam, shape1=alpha2, shape2=beta2)

  if (measure=="OR") return(p2/(1-p2)*(1-p1)/p1)
  if (measure=="RR") return(p2/p1)
  if (measure=="RD") return(p2-p1)
}

######################################################################################
### Purpose: sample from posterier distribution under Sarmanov beta prior distribution
### Input:   data(y1,n1,y2,n2), hyperparemeters(a1,a2,b1,b2,rho), comparative measure,
###          model, number of smaples (nsam)
### Output:  samples
### Note:    This function is called by the wrapper function "sampling"
### Author:  Sheng Luo and Xiao Su
### Data:    7/13/2012
######################################################################################
sar.sampling.posterior<-function(a1,a2,b1,b2,n1,y1,n2,y2,rho,measure,nsam) {
  code <- file.path(.find.package("mmeta"), "winbugcode", "SarmanovPosterior.txt")
  mu1 <- a1/(a1+b1); mu1.1 <- 1-mu1
  mu2 <- a2/(a2+b2); mu2.1 <- 1-mu2
  omega <- rho/sqrt(mu1*mu1.1/(a1+b1+1))/sqrt(mu2*mu2.1/(a2+b2+1))
  data <- list(rho=rho, alpha1=a1+y1, beta1=n1-y1+b1, alpha2=a2+y2, beta2=n2-y2+b2,
               a1=a1, b1=b1, a2=a2, b2=b2)
  bugsData(data, fileName="DataR.txt", digits=6)
  init <- list(p1=0.5, p2=0.5)
  bugsData(init, fileName="InitR.txt", digits=6)
  modelCheck(code)     ###hwo to set relative paht???
  modelData("DataR.txt")
  modelCompile(numChains=1)
  modelInits("InitR.txt")
  modelGenInits()
  modelUpdate(nsam)
  samplesSet(c("OR","RR","RD"))
  samplesSetThin(2)
  modelUpdate(2*nsam)                                                             

  if (measure=="OR") return(samplesSample(node=c("OR")))
  if(measure=="RR")  return(samplesSample(node=c("RR")))
  if(measure=="RD")  return(samplesSample(node=c("RD")))

  unlink("InitR.txt",recursive=F)                                                              
  unlink("DataR.txt",recursive=F)
}
