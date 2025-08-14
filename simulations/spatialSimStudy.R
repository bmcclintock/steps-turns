options("rgdal_show_exportToProj4_warnings"="none") # suppress annoying warnings
library(tidyverse)
library(raster)
library(rasterVis)
library(viridis)
library(ggplot2)
library(doFuture)
library(doRNG)
#remotes::install_version("RandomFieldsUtils", version= "1.2.5", repos = "http://cran.us.r-project.org") # most recent archived version; required by RandomFields
#remotes::install_version("RandomFields", version="3.3.14", repos = "http://cran.us.r-project.org") # most recent archived version; required by Rhabit
#remotes::install_version("geoR", version = "1.8-1", repos = "http://cran.us.r-project.org") # last version before RandomFields was removed from dependencies
#remotes::install_github("papayoun/Rhabit@31ddf44",dependencies = TRUE) # last commit before RandomFields was removed from dependencies
library(Rhabit)
#if(!requireNamespace("momentuHMM",quietly=TRUE) || packageVersion("momentuHMM")<2){
#  remotes::install_github("bmcclintock/momentuHMM@develop",dependencies = TRUE) # requires momentuHMM version >= 2.0.0
#}
library(momentuHMM)

source("simulations/plotDens.R")

nsims <- 100
nbAnimals <- 10
obsPerAnimal <- 1000
ncores <- 10

lim <- c(-1, 1, -1, 1)*100
cropExtent <- extent(lim)
resol <- 1
ncov <- 3

sc <- "low"
beta1 <- c(0*resol,0*resol,0*resol,-0.2*resol) # state1 resource selection coefficients for the spatial covariates
crw1 <- 0
sd_1 <- sqrt(0.5)
beta2 <- c(-2*resol,3*resol,2.5*resol,-0.2*resol) # state2 resource selection coefficients for the spatial covariates
crw2 <- 0.15
sd_2 <- sqrt(5)

beta0 <- matrix(c(-2.5,-2.5,
                  -0.5,0.5,
                  0.5,-0.5,
                  0.5,-0.5,
                  0,0),ncov+2,2,byrow=TRUE)

formula <- ~ cov1 + cov2 + cov3 + d2c
fixPar2 <- list(mu=c(NA,1,2,3,4,5,NA,6,7,8,9,10,NA,1,2,3,4,5,NA,6,7,8,9,10,11,11,12,11,11,12,NA,NA),
                delta=c(NA,NA))

spatialStudyplotDens <- list()
spatialStudyplotDens[[sc]] <- list()

set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")

results <- list()
covlist <- list()
origData <- list()
for(isim in 1:nsims){
  
  #######################
  ## Define covariates ##
  #######################
  # Generate ncov spatial covariates
  covlist[[isim]] <- list()
  for(i in 1:ncov) {
    covlist[[isim]][[i]] <- Rhabit::simSpatialCov(lim = lim, nu = 0.6, rho = 50, sigma2 = 0.1, 
                                          resol = resol, raster_like = TRUE)
  }
  
  # Include squared distance to origin as covariate
  xgrid <- seq(lim[1], lim[2], by=resol)
  ygrid <- seq(lim[3], lim[4], by=resol)
  xygrid <- expand.grid(xgrid,ygrid)
  dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/100
  covlist[[isim]][[4]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))
  
  # Compute utilization distribution for states 1 and 2
  UD1 <- Rhabit::getUD(covariates=covlist[[isim]], beta=beta1,log=TRUE)
  UD2 <- Rhabit::getUD(covariates=covlist[[isim]], beta=beta2,log=TRUE)
  
  # Plot covariates
  ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                   legend.title = element_text(size=12), legend.text = element_text(size=12))
  c1plot <- Rhabit::plotRaster(rhabitToRaster(covlist[[isim]][[1]]), scale.name = expression(c[1])) + ggtheme
  c2plot <- Rhabit::plotRaster(rhabitToRaster(covlist[[isim]][[2]]), scale.name = expression(c[2])) + ggtheme
  UD1plot <- Rhabit::plotRaster(rhabitToRaster(UD1), scale.name = expression(pi)) + ggtheme
  UD2plot <- Rhabit::plotRaster(rhabitToRaster(UD2), scale.name = expression(pi)) + ggtheme
  
  names(covlist[[isim]]) <- c("cov1","cov2","cov3","d2c")
  spatialCovs <- lapply(lapply(covlist[[isim]],rhabitToRaster),function(x) {proj4string(x) <- CRS("+init=epsg:3416");return(x)})
  for(j in c("cov1","cov2","cov3")){
    spatialCovs[[j]] <- scale(spatialCovs[[j]])
  }
  names(spatialCovs) <- c("cov1","cov2","cov3","d2c")
  for(i in names(spatialCovs)){
    names(spatialCovs[[i]]) <- i
  }
  
  dist <- list(mu="rw_mvnorm2")   # bivariate normal random walk
  
  # specify 2-state langevin pseudo-design matrix
  DM2 <- list(mu=list(mean.x=~0+mu.x_tm1+crw(mu.x_tm1)+langevin(cov1.x)+langevin(cov2.x)+langevin(cov3.x)+langevin(d2c.x),
                      mean.y=~0+mu.y_tm1+crw(mu.y_tm1)+langevin(cov1.y)+langevin(cov2.y)+langevin(cov3.y)+langevin(d2c.y),
                      sd.x=~state2(intercept),
                      sd.y=~state2(intercept),
                      corr.xy=~1))
  
  # simulate "high resolution" tracks; this can take a while...
  origData[[isim]] <- tryCatch(stop(),error=function(e) e)
  while(inherits(origData[[isim]],"error")){
    initState <- sample.int(2,nbAnimals,replace=TRUE)
    initPos <- matrix(c(sample(ncell(rhabitToRaster(UD1)),nbAnimals,replace=FALSE,prob=exp(getValues(rhabitToRaster(UD1)))/sum(exp(getValues(rhabitToRaster(UD1))))),
                        sample(ncell(rhabitToRaster(UD2)),nbAnimals,replace=FALSE,prob=exp(getValues(rhabitToRaster(UD2)))/sum(exp(getValues(rhabitToRaster(UD2)))))),2,nbAnimals,byrow=TRUE)
    initialPosition <- mapply(function(x) xyFromCell(rhabitToRaster(UD1),initPos[initState[x],x]),1:nbAnimals,SIMPLIFY = FALSE)
    origData[[isim]] <- tryCatch(simData(nbAnimals=nbAnimals, obsPerAnimal=obsPerAnimal, ncores=ncores,
                                  nbStates=2,
                                  initialPosition=initialPosition,
                                  dist=dist["mu"],
                                  DM=DM2,
                                  Par=list(mu=c(1,crw1,beta1,1,crw2,beta2,1,crw1,beta1,1,crw2,beta2,log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),0,0)),#,aux=c(10,2)),
                                  beta=beta0,
                                  formula=formula,
                                  formulaDelta=~0+ID,
                                  delta = matrix(ifelse(initState==1,-1.e+100,1.e+100),nbAnimals,1),
                                  spatialCovs = spatialCovs, 
                                  covs=data.frame(intercept=1),
                                  gradient=TRUE,
                                  mvnCoords = "mu",
                                  states=TRUE,
                                  TMB=TRUE),error=function(e) e) 
  }
  
  plotSpatialCov(origData[[isim]],spatialCov=rhabitToRaster(UD1),states=origData[[isim]]$states)
  
  Par0 <- list(mu=c(1,crw1,beta1,1,crw2,beta2,1,crw1,beta1,1,crw2,beta2,log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),0,0))
  
  fit <- suppressMessages(fitHMM(origData[[isim]],nbStates=2,
                                       dist=dist["mu"],
                                       DM=DM2,
                                       Par0=Par0,
                                       beta0=beta0,
                                       delta0=c(0.5,0.5),
                                       formula=formula,
                                       fixPar=fixPar2,
                                       mvnCoords = "mu",
                                       optMethod="TMB",
                                       control=list(silent=TRUE,maxit=10000,eval.max=10000)))
      
  results[[isim]] <- list()
  results[[isim]][[1]] <- list()
  results[[isim]][[1]]$CIbeta <- fit$CIbeta
  results[[isim]][[1]]$CIreal <- fit$CIreal
  results[[isim]][[1]]$states <- fit$data$states
  results[[isim]][[1]]$viterbi <- viterbi(fit)
  results[[isim]][[1]]$stateProbs <- stateProbs(fit)
  results[[isim]][[1]]$mod <- fit$mod$out
  
  for(j in names(spatialCovs)){
    origData[[isim]][[paste0(j,".w")]] <- sqrt(origData[[isim]][[paste0(j,".x")]]^2 + origData[[isim]][[paste0(j,".y")]]^2)
    origData[[isim]][[paste0(j,".d")]] <- origData[[isim]][[paste0(j,".r")]] <- atan2(origData[[isim]][[paste0(j,".y")]],origData[[isim]][[paste0(j,".x")]])
    class(origData[[isim]][[paste0(j,".d")]]) <- class(origData[[isim]][[paste0(j,".r")]]) <- "angle"
  }
  
  origData[[isim]] <- suppressMessages(prepData(origData[[isim]],coordNames=c("mu.x","mu.y"),angleCovs=paste0(names(spatialCovs),".d")))
  origData[[isim]]$bearing <- c(atan2(origData[[isim]]$y[2:nrow(origData[[isim]])]-origData[[isim]]$y[1:(nrow(origData[[isim]])-1)],origData[[isim]]$x[2:nrow(origData[[isim]])]-origData[[isim]]$x[1:(nrow(origData[[isim]])-1)]),NA)
  class(origData[[isim]]$bearing) <- "angle"
  
  for(j in names(spatialCovs)){
    origData[[isim]][[paste0(j,".angle")]] <- cos(origData[[isim]]$bearing-origData[[isim]][[paste0(j,".r")]])
    origData[[isim]][[paste0(j,".angle")]][is.na(origData[[isim]][[paste0(j,".angle")]])] <- 0
    class(origData[[isim]][[paste0(j,".angle")]]) <- "angle"
  }
  
  spatialStudyplotDens[[sc]][[isim]] <- plotDens(origData[[isim]],sc)
  
  ST <- suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),Par0=list(step=c(1.7,3,1,1.3),angle=c(0.01,0.3)),formula=formula,beta0=beta0,delta0=c(0.5,0.5),optMethod="TMB",control=list(silent=TRUE,maxit=10000,eval.max=10000),fixPar=list(delta=c(NA,NA))))
  results[[isim]][[2]] <- list()
  results[[isim]][[2]]$CIbeta <- ST$CIbeta
  results[[isim]][[2]]$CIreal <- ST$CIreal
  results[[isim]][[2]]$states <- ST$data$states
  results[[isim]][[2]]$viterbi <- viterbi(ST)
  results[[isim]][[2]]$stateProbs <- stateProbs(ST)
  results[[isim]][[2]]$AIC <- AIC(ST)
  results[[isim]][[2]]$mod <- ST$mod$out
  
  stpar <- getPar0(ST,DM=list(step=list(mean=~cov1+cov2+cov3+d2c,sd=~1),angle=list(concentration=~cov1+cov2+cov3+d2c)))
  STcov <- tryCatch(suppressWarnings(suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),DM=list(step=list(mean=formula,sd=~1),angle=list(concentration=formula)),
                                                             Par0=stpar$Par,formula=formula,beta0=stpar$beta,delta0=stpar$delta,optMethod="TMB",control=list(silent=TRUE,maxit=100000,eval.max=100000),fixPar=list(delta=c(NA,NA))))),error=function(e) e)
  if(inherits(STcov,"error") || STcov$mod$code>0){
    if(!inherits(STcov,"error")) stpar <- getPar(STcov)
    tmp <- tryCatch(suppressWarnings(suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),DM=list(step=list(mean=formula,sd=~1),angle=list(concentration=formula)),
                                                      Par0=stpar$Par,formula=formula,beta0=stpar$beta,delta0=stpar$delta,optMethod="TMB",control=list(silent=TRUE,maxit=100000,method="Nelder-Mead"),fixPar=list(delta=c(NA,NA))))),error=function(e) e)
    if(inherits(tmp,"error")){
      tmp <- tryCatch(suppressWarnings(suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),DM=list(step=list(mean=formula,sd=~1),angle=list(concentration=formula)),
                                                                 Par0=stpar$Par,formula=formula,beta0=stpar$beta,delta0=stpar$delta,nlmPar=list(iterlim=100000),fixPar=list(delta=stfullpar$delta)))),error=function(e) e)
      if(!inherits(tmp,"error")) STcov <- tmp
    } else STcov <- tmp
    stpar <- getPar(STcov)
    tmp <- tryCatch(suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),DM=list(step=list(mean=formula,sd=~1),angle=list(concentration=formula)),
                                     Par0=stpar$Par,formula=formula,beta0=stpar$beta,delta0=stpar$delta,optMethod="TMB",control=list(silent=TRUE,maxit=100000,eval.max=100000),fixPar=list(delta=c(NA,NA)))),error=function(e) e)
    if(!inherits(tmp,"error")) STcov <- tmp
  }
  
  results[[isim]][[3]] <- list()
  results[[isim]][[3]]$CIbeta <- STcov$CIbeta
  results[[isim]][[3]]$CIreal <- STcov$CIreal
  results[[isim]][[3]]$states <- STcov$data$states
  results[[isim]][[3]]$viterbi <- viterbi(STcov)
  results[[isim]][[3]]$stateProbs <- stateProbs(STcov)
  results[[isim]][[3]]$AIC <- AIC(STcov)
  if(STcov$conditions$optMethod=="TMB") results[[isim]][[3]]$mod <- STcov$mod$out
  else results[[isim]][[3]]$mod <- STcov$mod
  
  DMSTfull<-list(step=list(mean=~cov1.w:cov1.angle+cov2.w:cov2.angle+cov3.w:cov3.angle+d2c.w:d2c.angle,sd=~1),
                 angle=list(mean=~cov1.d+cov2.d+cov3.d+d2c.d,concentration=~1))
  stfullpar <- getPar0(ST,DM=DMSTfull,estAngleMean = list(angle=TRUE),circularAngleMean = list(angle=TRUE))
  stfullpar$Par$step[c(2:5,7:10)] <- c(beta1,beta2) * 0.25
  stfullpar$Par$angle[1:10] <- c(c(beta1,beta2) * .8,-0.4,0.15)
  STfull <- tryCatch(suppressWarnings(suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),DM=DMSTfull,estAngleMean = list(angle=TRUE),circularAngleMean = list(angle=TRUE),
                                                              Par0=stfullpar$Par,formula=formula,beta0=stfullpar$beta,delta0=stfullpar$delta,nlmPar=list(iterlim=100000),fixPar=list(delta=stfullpar$delta)))),error=function(e) e)
  
  if(inherits(STfull,"error") || STfull$mod$code>0){
    if(!inherits(STfull,"error")) stfullpar <- getPar(STfull)
    tmp <- tryCatch(suppressWarnings(suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),DM=DMSTfull,estAngleMean = list(angle=TRUE),circularAngleMean = list(angle=TRUE),
                                                             Par0=stfullpar$Par,formula=formula,beta0=stfullpar$beta,delta0=stfullpar$delta,optMethod="Nelder-Mead",control=list(maxit=100000,hessian=FALSE),fixPar=list(delta=stfullpar$delta)))),error=function(e) e)
    if(!inherits(tmp,"error")) STfull <- tmp
    stfullpar <- getPar(STfull)
    tmp <- tryCatch(suppressMessages(fitHMM(origData[[isim]],nbStates=2,dist=list(step="gamma",angle="vm"),DM=DMSTfull,estAngleMean = list(angle=TRUE),circularAngleMean = list(angle=TRUE),
                                            Par0=stfullpar$Par,formula=formula,beta0=stfullpar$beta,delta0=stfullpar$delta,nlmPar=list(iterlim=100000),fixPar=list(delta=stfullpar$delta))),error=function(e) e)
    if(!inherits(tmp,"error")) STfull <- tmp
  }
  
  results[[isim]][[4]] <- list()
  results[[isim]][[4]]$CIbeta <- STfull$CIbeta
  results[[isim]][[4]]$CIreal <- STfull$CIreal
  results[[isim]][[4]]$states <- STfull$data$states
  results[[isim]][[4]]$viterbi <- viterbi(STfull)
  results[[isim]][[4]]$stateProbs <- stateProbs(STfull)
  results[[isim]][[4]]$AIC <- AIC(STfull)
  results[[isim]][[4]]$mod <- STfull$mod
  
  cat(isim,mean(results[[isim]][[1]]$viterbi==results[[isim]][[1]]$states),mean(results[[isim]][[2]]$viterbi==results[[isim]][[2]]$states),mean(results[[isim]][[3]]$viterbi==results[[isim]][[3]]$states),mean(results[[isim]][[4]]$viterbi==results[[isim]][[4]]$states),"\n")
  rm(fit,ST,STcov,STfull)
}
save(beta0,beta1,beta2,crw1,crw2,results,sd_1,sd_2,file=paste0("simulations/spatialSimStudy_",sc,".RData"))
save(spatialStudyplotDens,file="simulations/spatialStudyplotDens.RData")


spatialplot <- list()
spatialplot$viterbi <- spatialplot$beta12 <- spatialplot$beta21 <- list()
for(sc in c("low","moderate","high")){
  load(paste0("simulations/spatialSimStudy_",sc,".RData"))
  sdf <- data.frame(model=factor(rep(c("lang","ST","STcov","STfull"),each=nsims)))
  sdf$viterbi <- c(unlist(lapply(results,function(x) mean(x[[1]]$viterbi==x[[1]]$states))),
                   unlist(lapply(results,function(x) mean(x[[2]]$viterbi==x[[2]]$states))),
                   unlist(lapply(results,function(x) mean(x[[3]]$viterbi==x[[3]]$states))),
                   unlist(lapply(results,function(x) mean(x[[4]]$viterbi==x[[4]]$states))))
  #vioplot::vioplot(viterbi~model,data=sdf)
  #abline(h=0.95,col=2,lty=2)
  spatialplot$viterbi[[sc]] <- ggplot(sdf, aes(x=model, y=viterbi)) + ylab("Proportion correct state assignments")+
       geom_violin(trim=TRUE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=0.95,linetype=2,col=2)+ scale_y_continuous(limits = c(0.4,1)) 
  covNames <- c("cov1","cov2","cov3","d2c")
  spatialplot$beta12[[sc]] <- spatialplot$beta21[[sc]] <- list()
  for(j in 1:length(covNames)){
    cov <- covNames[j]
    sdf[[paste0("beta12",cov)]] <- c(unlist(lapply(results,function(x) x[[1]]$CIbeta$beta$est[1+j])),
                                     unlist(lapply(results,function(x) x[[2]]$CIbeta$beta$est[1+j])),
                                     unlist(lapply(results,function(x) x[[3]]$CIbeta$beta$est[1+j])),
                                     unlist(lapply(results,function(x) x[[4]]$CIbeta$beta$est[1+j])))
    sdf[[paste0("beta21",cov)]] <- c(unlist(lapply(results,function(x) x[[1]]$CIbeta$beta$est[6+j])),
                                     unlist(lapply(results,function(x) x[[2]]$CIbeta$beta$est[6+j])),
                                     unlist(lapply(results,function(x) x[[3]]$CIbeta$beta$est[6+j])),
                                     unlist(lapply(results,function(x) x[[4]]$CIbeta$beta$est[6+j])))
    #vioplot::vioplot(as.formula(paste0("beta12",cov,"~model")),data=sdf)
    #abline(h=beta0[1+j],col=2,lty=2)
    spatialplot$beta12[[sc]][[cov]] <- ggplot(sdf, aes(x=model, y=!!sym(paste0("beta12",cov)))) + ylab(expression(beta[12]))+
      geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[1+j],linetype=2,col=2)+ scale_y_continuous(limits = c(-1.75,1.75)) 
    #vioplot::vioplot(as.formula(paste0("beta21",cov,"~model")),data=sdf)
    #abline(h=beta0[6+j],col=2,lty=2)
    spatialplot$beta21[[sc]][[cov]] <- ggplot(sdf, aes(x=model, y=!!sym(paste0("beta21",cov)))) + ylab(expression(beta[21]))+
      geom_violin(trim=FALSE,fill="grey")+ geom_boxplot(width=0.1)+geom_hline(yintercept=beta0[6+j],linetype=2,col=2)+ scale_y_continuous(limits = c(-1.75,1.75)) 
  }
}
save(spatialplot,file="simulations/spatialSimStudyplots.RData")
