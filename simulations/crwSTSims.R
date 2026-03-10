if(!requireNamespace("momentuHMM",quietly=TRUE) || packageVersion("momentuHMM")<"2.0.0"){
  remotes::install_github("bmcclintock/momentuHMM@develop",dependencies = TRUE) # requires momentuHMM version >= 2.0.0
}
library(momentuHMM)
library(raster)
library(fields)
#remotes::install_github("papayoun/Rhabit") 
library(Rhabit)

source("simulations/plotDens.R")
source("simulations/simCov.R")

nSims <- 100
n <- 1000 # obsPerAnimal
nbAnimals <- 10

scenarios <- c("low","moderate","high")

beta_1 <- c(10,1.4,1.75)
beta_2 <- c(0.275,.3,0.5)
sigma_1 <- c(2.7,0.6,5)
sigma_2 <- c(2.7,1.85,1)

crw1_1 <- c(0.1,0.5,0.25)
crw2_1 <- c(0,0,0) # lag 2 correlation
crw1_2 <- c(0.85,0.85,0.75)
crw2_2 <- c(0,0,0) # lag 2 correlation
crwsigma_1 <- c(0.25,0.25,0.75)
crwsigma_2 <- c(1.5,1,0.25)

mu_1 <- c(0.35,0.5,0.75)
mu_2 <- c(3,2,0.8)
sd_1 <- c(0.2,0.3,0.5)
sd_2 <- c(1.6,1,0.2)
rho_1 <- c(0.2,1,0.5)
rho_2 <- c(1.6,3,7.5)

formula <- ~cov
beta0 <- matrix(c(-1.5,-1.5,-.5,.5),2,2,byrow=TRUE)
delta0 <- c(0.5,0.5)

crwST <- ST <- pos <- crwpos <- datcrwST <- list()

# define extent for spatial covariate
x_seq <- seq(-350, 350, length.out = 70)
y_seq <- seq(-350, 350, length.out = 70)

posDM <- list(mu=list(mean.x=~1,
                      mean.y=~1,
                      sd.x=~1,
                      sd.y=~1,
                      corr.xy=~1))

poscrwDM <- list(mu=list(mean.x=~crw(mu.x_tm1),
                         mean.y=~crw(mu.y_tm1),
                         sd.x=~1,
                         sd.y=~1,
                         corr.xy=~1))

extractRes <- function(mod){
  list(mod=mod$mod[c("out","minimum","estimate","code","wpar","hessian","Sigma")],
       CIbeta=mod$CIbeta,
       CIreal=mod$CIreal,
       viterbi=viterbi(mod),
       stateProbs=stateProbs(mod),
       AIC=AIC(mod),
       states=mod$data$states,
       pseudoRes=pseudoRes(mod),
       plot=suppressMessages(plot(mod,plotCI=TRUE,plotStationary=TRUE,ask=FALSE,return=TRUE)$estimates))
}

crwSTplotDens <- list()

set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")

for(sc in 1:length(scenarios)){
  crwSTplotDens[[scenarios[sc]]] <- list()
  for(isim in 1:nSims){
    
    # generate spatial covariate
    spatialCov <- simCov(x_seq,y_seq)
    
    values(spatialCov$cov) <- scale(values(spatialCov$cov))
    
    datcrwST[[isim]] <- suppressMessages(simData(nbStates=2,nbAnimals=nbAnimals,obsPerAnimal=n,dist=list(step="crwrice",angle="crwvm"),formula=formula,Par=list(step=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc]),angle=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc])),beta=beta0,delta=delta0,spatialCovs=spatialCov,states=TRUE,TMB=TRUE,ncores=nbAnimals,retrySims=100))
    
    datcrwST[[isim]]$mu.x <- datcrwST[[isim]]$x
    datcrwST[[isim]]$mu.y <- datcrwST[[isim]]$y
    
    crwSTplotDens[[scenarios[sc]]][[isim]] <- plotDens(datcrwST[[isim]],scenarios[sc])
    
    plotSpatialCov(datcrwST[[isim]],spatialCov$cov,states=datcrwST[[isim]]$states)
    
    crwSTmod <- tryCatch(suppressMessages(fitHMM(datcrwST[[isim]],nbStates=2,dist=list(step="crwrice",angle="crwvm"),formula=formula,Par0=list(step=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc]),angle=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc])),
                            beta0=beta0,delta0=delta0,optMethod="TMB",control=list(silent=TRUE))),error=function(e) e)
    while(inherits(crwSTmod,"error") || inherits(crwSTmod$mod$Sigma,"error")){
      if(inherits(crwSTmod$mod$Sigma,"error")) parcrwST <- getPar(crwSTmod)
      else parcrwST <- list(Par=list(step=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc]),angle=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc])),beta=beta0,delta=delta0)
      crwSTmod <- tryCatch(suppressMessages(fitHMM(datcrwST[[isim]],nbStates=2,dist=list(step="crwrice",angle="crwvm"),formula=formula,Par0=parcrwST$Par,
                                                   beta0=parcrwST$beta,delta0=parcrwST$delta,retryFits=5,retrySD=list(step=0.1,beta=0,delta=0))),error=function(e) e)
    }
    
    crwST[[isim]] <- extractRes(crwSTmod)
  
    
    STmod <- suppressMessages(fitHMM(datcrwST[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),formula=formula,Par0=list(step=c(mu_1[sc],mu_2[sc],sd_1[sc],sd_2[sc]),angle=c(rho_1[sc],rho_2[sc])),
                                     beta0=beta0,delta0=delta0,optMethod = "TMB",control=list(silent=TRUE)))
    
    ST[[isim]] <- extractRes(STmod)
    
    
    posmod <- suppressMessages(fitHMM(datcrwST[[isim]],nbStates=2,dist=list(mu="rw_mvnorm2"),mvnCoords = "mu",DM=posDM,formula=formula,Par0=list(mu=c(1,1,1,1,log(crwsigma_1[sc]),log(crwsigma_2[sc]),log(crwsigma_1[sc]),log(crwsigma_2[sc]),0,0)),
                                      beta0=beta0,delta0=delta0,optMethod="TMB",fixPar=list(mu=c(NA,NA,NA,NA,1,2,1,2,NA,NA)),control=list(silent=TRUE)))
    
    pos[[isim]] <- extractRes(posmod)
    
    
    parcrw <- getPar0(posmod,DM=poscrwDM)
    crwposmod <- suppressMessages(fitHMM(datcrwST[[isim]],nbStates=2,dist=list(mu="rw_mvnorm2"),mvnCoords = "mu",DM=poscrwDM,formula=formula,Par0=parcrw$Par,
                                         beta0=parcrw$beta,delta0=parcrw$delta,optMethod="TMB",fixPar=list(mu=c(NA,1,NA,2,NA,1,NA,2,3,4,3,4,NA,NA)),control=list(silent=TRUE)))
    
    crwpos[[isim]] <- extractRes(crwposmod)
    
    cat(isim,mean(crwST[[isim]]$viterbi==crwST[[isim]]$states),mean(ST[[isim]]$viterbi==ST[[isim]]$states),mean(pos[[isim]]$viterbi==pos[[isim]]$states),mean(crwpos[[isim]]$viterbi==crwpos[[isim]]$states),"\n")
    rm(spatialCov,crwSTmod,STmod,posmod,crwposmod)
  }
  save(crwST,ST,pos,crwpos,file=paste0("simulations/crwSTSims_",scenarios[sc],".RData"))
}
save(crwSTplotDens,file="simulations/crwSTplotDens.RData")

crwSTplot <- list()
crwSTplot$viterbi <- crwSTplot$beta12 <- crwSTplot$beta21 <- list()
for(sc in 1:length(scenarios)){
  load(paste0("simulations/crwSTSims_",scenarios[sc],".RData"))
  sdf <- data.frame(model=factor(rep(c("crwST","ST",#"pos",
                                       "crwpos"),each=nSims)))
  sdf$viterbi <- c(unlist(lapply(crwST,function(x) mean(x$viterbi==x$states))),
                   unlist(lapply(ST,function(x) mean(x$viterbi==x$states))),
                   #unlist(lapply(pos,function(x) mean(x$viterbi==x$states))),
                   unlist(lapply(crwpos,function(x) mean(x$viterbi==x$states))))
  sdf$beta12 <- c(unlist(lapply(crwST,function(x) x$CIbeta$beta$est[2])),
                  unlist(lapply(ST,function(x) x$CIbeta$beta$est[2])),
                  #unlist(lapply(pos,function(x) x$CIbeta$beta$est[2])),
                  unlist(lapply(crwpos,function(x) x$CIbeta$beta$est[2])))
  sdf$beta21 <- c(unlist(lapply(crwST,function(x) x$CIbeta$beta$est[4])),
                  unlist(lapply(ST,function(x) x$CIbeta$beta$est[4])),
                  #unlist(lapply(pos,function(x) x$CIbeta$beta$est[4])),
                  unlist(lapply(crwpos,function(x) x$CIbeta$beta$est[4])))
  
  #vioplot::vioplot(viterbi~model,data=sdf)
  #abline(h=0.95,col=2,lty=2)
  crwSTplot$viterbi[[scenarios[sc]]] <- ggplot(sdf, aes(x=model, y=viterbi)) + ylab("Proportion correct state assignments")+
    geom_violin(trim=TRUE,fill="grey")+ geom_boxplot(width=0.05)+geom_hline(yintercept=0.95,linetype=2,col=2) + scale_y_continuous(limits = c(0.5,1)) 
  #vioplot::vioplot(beta12~model,data=sdf)
  #abline(h=beta0[2],col=2,lty=2)
  crwSTplot$beta12[[scenarios[sc]]] <- ggplot(sdf, aes(x=model, y=beta12)) +ylab(expression(beta[12]))+
    geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[2],linetype=2,col=2)+ scale_y_continuous(limits = c(-1.3,0.3)) 
  #vioplot::vioplot(beta21~model,data=sdf)
  #abline(h=beta0[4],col=2,lty=2)
  crwSTplot$beta21[[scenarios[sc]]] <- ggplot(sdf, aes(x=model, y=beta21)) +ylab(expression(beta[21]))+
    geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[4],linetype=2,col=2)+ scale_y_continuous(limits = c(-0.25,1))
}
save(crwSTplot,file="simulations/crwSTplots.RData")
