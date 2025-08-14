library(momentuHMM)
library(raster)
library(geoR)

nSims <- 100
n <- 1000 # obsPerAnimals
nbAnimals <- 10

scenarios <- c("low","moderate","high")

beta_1 <- c(0.05,0.5,0.4)
beta_2 <- c(0.8,0.8,0.6)
sigma_1 <- c(0.1,1/3,0.5)
sigma_2 <- c(1,2/3,0.75)

crw1_1 <- c(0.1,0.5,0.25)
crw2_1 <- c(0,0,0) # lag 2 correlation
crw1_2 <- c(0.85,0.85,0.75)
crw2_2 <- c(0,0,0) # lag 2 correlation
crwsigma_1 <- c(0.25,0.25,.25)
crwsigma_2 <- c(1.5,1,0.75)


formula <- ~cov
beta0 <- matrix(c(-3.5,-3.5,-.5,.5),2,2,byrow=TRUE)
delta0 <- c(0.5,0.5)

ctcrw <- pos <- crwpos <- datctcrw <- list()

posDM <- list(mu=list(mean.x=~1,
                      mean.y=~1,
                      sd.x=~1,
                      sd.y=~1,
                      corr.xy=~1))

poscrwDM <- list(mu=list(mean.x=~crw(mu.x_tm1,dt=dt),
                         mean.y=~crw(mu.y_tm1,dt=dt),
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

set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")

for(sc in 1:length(scenarios)){
  for(isim in 1:nSims){
    
    spatialCov <- list(cov=raster::raster(geoR::grf(70*70, xlims =c(-350,350), ylims=c(-350,350),grid = "reg", kappa=0.6,cov.pars = c(0.1, 50))))
    
    values(spatialCov$cov) <- scale(values(spatialCov$cov))
    
    datctcrw[[isim]] <- suppressMessages(simCTHMM(nbStates=2,nbAnimals=nbAnimals,obsPerAnimal=n,dist=list(mu="ctcrw"),mvnCoords="mu",formula=formula,Par=list(mu=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc])),beta=beta0,delta=delta0,spatialCovs=spatialCov,states=TRUE,TMB=TRUE,ncores=nbAnimals,retrySims=100))
    
    plotSpatialCov(datctcrw[[isim]],spatialCov$cov,states=datctcrw[[isim]]$states)
    
    ctcrwmod <- suppressMessages(fitCTHMM(datctcrw[[isim]],nbStates=2,dist=list(mu="ctcrw"),mvnCoords="mu",formula=formula,Par0=list(mu=c(beta_1[sc],beta_2[sc],sigma_1[sc],sigma_2[sc])),beta0=beta0,delta0=delta0,optMethod="TMB",control=list(silent=TRUE)))
    
    ctcrw[[isim]] <- extractRes(ctcrwmod)
    
    posmod <- suppressMessages(fitCTHMM(datctcrw[[isim]],nbStates=2,dist=list(mu="rw_mvnorm2"),mvnCoords = "mu",DM=posDM,formula=formula,Par0=list(mu=c(1,1,1,1,log(crwsigma_1[sc]),log(crwsigma_2[sc]),log(crwsigma_1[sc]),log(crwsigma_2[sc]),0,0)),
                                      beta0=beta0,delta0=delta0,optMethod="TMB",fixPar=list(mu=c(NA,NA,NA,NA,1,2,1,2,NA,NA)),control=list(silent=TRUE)))
    
    pos[[isim]] <- extractRes(posmod)
    
    
    parcrw <- getPar0(posmod,DM=poscrwDM)
    crwposmod <- suppressMessages(fitCTHMM(datctcrw[[isim]],nbStates=2,dist=list(mu="rw_mvnorm2"),mvnCoords = "mu",DM=poscrwDM,formula=formula,Par0=parcrw$Par,
                                         beta0=parcrw$beta,delta0=parcrw$delta,optMethod="TMB",fixPar=list(mu=c(NA,1,NA,2,NA,1,NA,2,3,4,3,4,NA,NA)),control=list(silent=TRUE)))
    
    crwpos[[isim]] <- extractRes(crwposmod)
    
    cat(isim,mean(ctcrw[[isim]]$viterbi==ctcrw[[isim]]$states),mean(pos[[isim]]$viterbi==pos[[isim]]$states),mean(crwpos[[isim]]$viterbi==crwpos[[isim]]$states),"\n")
    rm(spatialCov,ctcrwmod,posmod,crwposmod)
  }
  save(ctcrw,pos,crwpos,file=paste0("simulations/ctcrwSims_",scenarios[sc],".RData"))
}

sdf <- data.frame(model=factor(rep(c("ctcrw","pos","crwpos"),each=nSims)))
sdf$viterbi <- c(unlist(lapply(ctcrw,function(x) mean(x$viterbi==x$states))),
                 unlist(lapply(pos,function(x) mean(x$viterbi==x$states))),
                 unlist(lapply(crwpos,function(x) mean(x$viterbi==x$states))))
sdf$beta12 <- c(unlist(lapply(ctcrw,function(x) x$CIbeta$beta$est[2])),
                unlist(lapply(pos,function(x) x$CIbeta$beta$est[2])),
                unlist(lapply(crwpos,function(x) x$CIbeta$beta$est[2])))
sdf$beta21 <- c(unlist(lapply(ctcrw,function(x) x$CIbeta$beta$est[4])),
                unlist(lapply(pos,function(x) x$CIbeta$beta$est[4])),
                unlist(lapply(crwpos,function(x) x$CIbeta$beta$est[4])))

#vioplot::vioplot(viterbi~model,data=sdf)
#abline(h=0.95,col=2,lty=2)
ggplot(sdf, aes(x=model, y=viterbi)) + 
  geom_violin(trim=TRUE,fill="grey")+ geom_boxplot(width=0.05)+geom_hline(yintercept=0.95,linetype=2,col=2)
#vioplot::vioplot(beta12~model,data=sdf)
#abline(h=beta0[2],col=2,lty=2)
ggplot(sdf, aes(x=model, y=beta12)) +
  geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[2],linetype=2,col=2)
#vioplot::vioplot(beta21~model,data=sdf)
#abline(h=beta0[4],col=2,lty=2)
ggplot(sdf, aes(x=model, y=beta21)) +
  geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[4],linetype=2,col=2)

