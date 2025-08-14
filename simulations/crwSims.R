library(momentuHMM)
library(raster)
#remotes::install_version("RandomFieldsUtils", version= "1.2.5", repos = "http://cran.us.r-project.org") # most recent archived version; required by RandomFields
#remotes::install_version("RandomFields", version="3.3.14", repos = "http://cran.us.r-project.org") # most recent archived version
library(RandomFields)
RFoptions(spConform = FALSE)
source("simulations/plotDens.R")

nSims <- 100
n <- 1000 # obsPerAnimal
nbAnimals <- 10

scenarios <- c("low","moderate","high")

# Cartesian correlated random walk parameters
crw1_1 <- c(0.1,0.5,0.25) # state 1 lag 1 correlation (low, moderate, high)
crw2_1 <- c(0,0,0) # state 2 lag 2 correlation
crw1_2 <- c(0.85,0.85,0.75) # state 1 log 1 correlation
crw2_2 <- c(0,0,0) # state 2 lag 2 correlation
crwsigma_1 <- c(0.25,0.25,0.75) # state 1 sqrt(variance)
crwsigma_2 <- c(1.5,1,0.25) # state 2 sqrt(variance)

# correlated step and turn initial values
beta_1 <- c(10,1.4,1.75) # state 1 autocorrelation (low, moderate, high)
beta_2 <- c(0.275,.3,0.5) # state 2 autocorrelation
sigma_1 <- c(2.7,0.6,2) # state 1 sqrt(variance)
sigma_2 <- c(2.7,1.85,0.5) # state 2 sqrt(variance)

# conventional step and turn initial values
mu_1 <- c(0.35,0.5,0.75) # state 1 step length mean (low, moderate, high)
mu_2 <- c(3,2,0.8) # state 2 step length mean
sd_1 <- c(0.2,0.3,0.5) # state 1 step length sd
sd_2 <- c(1.6,1,0.2) # state 2 step length sd
rho_1 <- c(0.2,1,0.5) # state 1 turn angle concentration parameter
rho_2 <- c(1.6,3,7.5) # state 2 turn angle concentration parameter

formula <- ~cov # state transition  probability formula
beta0 <- matrix(c(-1.5,-1.5,-.5,.5),2,2,byrow=TRUE) # logit-scale state transition probability parameters
delta0 <- c(0.5,0.5) # initial distribution

crwST <- ST <- pos <- crwpos <- datcrwpos <- list()

# define extent for spatial covariate
x_seq <- seq(-350, 350, length.out = 70)
y_seq <- seq(-350, 350, length.out = 70)

simDM <- list(mu=list(mean.x=~crw(mu.x_tm1)+crw(mu.x_tm1,lag=2),
                      mean.y=~crw(mu.y_tm1)+crw(mu.y_tm1,lag=2),
                      sd.x=~1,
                      sd.y=~1,
                      corr.xy=~1))

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

crwplotDens <- list()

set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")

for(sc in 1:length(scenarios)){
  crwplotDens[[scenarios[sc]]] <- list()
  for(isim in 1:nSims){
    
    # generate spatial covariate
    spatialCov <- list(cov=raster::raster(RandomFields::RFsimulate(
      model =  RandomFields::RMmatern(var = 0.1, scale = 50, nu = 0.6), 
      x = x_seq, 
      y = y_seq, 
      grid = TRUE
    ),xmn=min(x_seq),xmx=max(x_seq),ymn=min(y_seq),ymx=max(y_seq)))
    
    # scale covariate
    values(spatialCov$cov) <- scale(values(spatialCov$cov))
    
    datcrwpos[[isim]] <- tryCatch(suppressMessages(simData(nbStates=2,nbAnimals=nbAnimals,obsPerAnimal=n,dist=list(mu="rw_mvnorm2"),spatialCovs=spatialCov,DM=simDM,formula=formula,Par=list(mu=c(1,crw1_1[sc],crw2_1[sc],1,crw1_2[sc],crw2_2[sc],1,crw1_1[sc],crw2_1[sc],1,crw1_2[sc],crw2_2[sc],log(crwsigma_1[sc]),log(crwsigma_2[sc]),log(crwsigma_1[sc]),log(crwsigma_2[sc]),0,0)),beta=beta0, delta=delta0, TMB=TRUE, mvnCoords = "mu",states=TRUE,ncores=nbAnimals,retrySims=100)),error=function(e) e)
    while(inherits(datcrwpos[[isim]],"error")){
      datcrwpos[[isim]] <- tryCatch(suppressMessages(simData(nbStates=2,nbAnimals=nbAnimals,obsPerAnimal=n,dist=list(mu="rw_mvnorm2"),spatialCovs=spatialCov,DM=simDM,formula=formula,Par=list(mu=c(1,crw1_1[sc],crw2_1[sc],1,crw1_2[sc],crw2_2[sc],1,crw1_1[sc],crw2_1[sc],1,crw1_2[sc],crw2_2[sc],log(crwsigma_1[sc]),log(crwsigma_2[sc]),log(crwsigma_1[sc]),log(crwsigma_2[sc]),0,0)),beta=beta0, delta=delta0, TMB=TRUE, mvnCoords = "mu",states=TRUE,ncores=nbAnimals,retrySims=100)),error=function(e) e)
    }
    datcrwpos[[isim]] <- suppressMessages(prepData(datcrwpos[[isim]],coordNames=c("mu.x","mu.y"),altCoordNames = "mu"))
    
    crwplotDens[[scenarios[sc]]][[isim]] <- plotDens(datcrwpos[[isim]],scenarios[sc])
    
    plotSpatialCov(datcrwpos[[isim]],spatialCov$cov,states=datcrwpos[[isim]]$states)
    
    crwSTmod <- tryCatch(suppressMessages(fitHMM(datcrwpos[[isim]],nbStates=2,dist=list(step="crwrice",angle="crwvm"),formula=formula,Par0=list(step=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc]),angle=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc])),
                            beta0=beta0,delta0=delta0,optMethod="TMB",control=list(silent=TRUE))),error=function(e) e)
    while(inherits(crwSTmod,"error") || inherits(crwSTmod$mod$Sigma,"error")){
      if(inherits(crwSTmod$mod$Sigma,"error")) parcrwST <- getPar(crwSTmod)
      else parcrwST <- list(Par=list(step=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc]),angle=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc])),beta=beta0,delta=delta0)
      crwSTmod <- tryCatch(suppressMessages(fitHMM(datcrwpos[[isim]],nbStates=2,dist=list(step="crwrice",angle="crwvm"),formula=formula,Par0=parcrwST$Par,
                                                   beta0=parcrwST$beta,delta0=parcrwST$delta,retryFits=5,retrySD=list(step=0.1,beta=0,delta=0))),error=function(e) e)
    }
    crwST[[isim]] <- extractRes(crwSTmod)
  
    
    STmod <- suppressMessages(fitHMM(datcrwpos[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),formula=formula,Par0=list(step=c(mu_1[sc],mu_2[sc],sd_1[sc],sd_2[sc]),angle=c(rho_1[sc],rho_2[sc])),
                         beta0=beta0,delta0=delta0,optMethod = "TMB",control=list(silent=TRUE)))
    
    ST[[isim]] <- extractRes(STmod)
    
    
    posmod <- suppressMessages(fitHMM(datcrwpos[[isim]],nbStates=2,dist=list(mu="rw_mvnorm2"),mvnCoords = "mu",DM=posDM,formula=formula,Par0=list(mu=c(1,1,1,1,log(crwsigma_1[sc]),log(crwsigma_2[sc]),log(crwsigma_1[sc]),log(crwsigma_2[sc]),0,0)),
                          beta0=beta0,delta0=delta0,optMethod="TMB",fixPar=list(mu=c(NA,NA,NA,NA,1,2,1,2,NA,NA)),control=list(silent=TRUE)))
    
    pos[[isim]] <- extractRes(posmod)
    
    
    parcrw <- getPar0(posmod,DM=poscrwDM)
    crwposmod <- suppressMessages(fitHMM(datcrwpos[[isim]],nbStates=2,dist=list(mu="rw_mvnorm2"),mvnCoords = "mu",DM=poscrwDM,formula=formula,Par0=parcrw$Par,
                             beta0=parcrw$beta,delta0=parcrw$delta,optMethod="TMB",fixPar=list(mu=c(NA,1,NA,2,NA,1,NA,2,3,4,3,4,NA,NA)),control=list(silent=TRUE)))
  
    crwpos[[isim]] <- extractRes(crwposmod)
    
    cat(isim,mean(crwST[[isim]]$viterbi==crwST[[isim]]$states),mean(ST[[isim]]$viterbi==ST[[isim]]$states),mean(pos[[isim]]$viterbi==pos[[isim]]$states),mean(crwpos[[isim]]$viterbi==crwpos[[isim]]$states),"\n")
    rm(spatialCov,crwSTmod,STmod,posmod,crwposmod)
  }
  save(crwST,ST,pos,crwpos,file=paste0("simulations/crwSims_",scenarios[sc],".RData"))
}
save(crwplotDens,file="simulations/crwplotDens.RData")

crwplot <- list()
crwplot$viterbi <- crwplot$beta12 <- crwplot$beta21 <- list()
for(sc in 1:length(scenarios)){
  load(paste0("simulations/crwSims_",scenarios[sc],".RData"))  
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
  crwplot$viterbi[[scenarios[sc]]] <- ggplot(sdf, aes(x=model, y=viterbi)) + ylab("Proportion correct state assignments")+
    geom_violin(trim=TRUE,fill="grey")+ geom_boxplot(width=0.05)+geom_hline(yintercept=0.95,linetype=2,col=2) + scale_y_continuous(limits = c(0.5,1)) 
  #vioplot::vioplot(beta12~model,data=sdf)
  #abline(h=beta0[2],col=2,lty=2)
  crwplot$beta12[[scenarios[sc]]] <- ggplot(sdf, aes(x=model, y=beta12)) +ylab(expression(beta[12]))+
         geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[2],linetype=2,col=2)+ scale_y_continuous(limits = c(-1.1,0.3)) 
  #vioplot::vioplot(beta21~model,data=sdf)
  #abline(h=beta0[4],col=2,lty=2)
  crwplot$beta21[[scenarios[sc]]] <- ggplot(sdf, aes(x=model, y=beta21)) +ylab(expression(beta[21]))+
    geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[4],linetype=2,col=2)+ scale_y_continuous(limits = c(-0.1,0.9)) 
}
save(crwplot,file="simulations/crwplots.RData")
