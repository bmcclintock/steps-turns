if(!requireNamespace("momentuHMM",quietly=TRUE) || packageVersion("momentuHMM")<"2.0.0"){
  remotes::install_github("bmcclintock/momentuHMM@develop",dependencies = TRUE) # requires momentuHMM version >= 2.0.0
}
library(momentuHMM)
library(sp)
library(raster)
library(cylcop)
library(CircStats)
library(marmap) 
library(maps) 
library(mapdata) 
library(giscoR)
library(sf)
library(tidyverse)
library(ggplot2)
library(qqplotr)
library(patchwork)
library(caret)
library(stringr)

source("examples/ses/sesFunctions.R")

# load ses data from github
tracks <- read.csv(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/sesData.csv"))
subSeals <- which(tracks$ID %in% c(2,3,4,7,11,14))
tracks <- tracks[subSeals,c(1,2,3)]
colnames(tracks) <- c("ID","x","y")

centerLL <- c(70,-49) # center long/lat coordinates

# add activity center to tracks for projection to epsg:3032
tracks[nrow(tracks)+1,]<-c(NA,centerLL)
coordinates(tracks)<-c("x","y")
proj4string(tracks)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
myCRS <- CRS("+init=epsg:3032 +units=km")

projCoords<-spTransform(tracks,CRS=myCRS)
projCoords<-data.frame(projCoords)
center <- matrix(as.numeric(projCoords[nrow(projCoords),c("coords.x1","coords.x2")]),nrow=1,dimnames=list("colony"))

# remove activity center from tracks
tracks <- tracks[-nrow(tracks),]
projCoords <- projCoords[-nrow(projCoords),]

meanx <- mean(projCoords$coords.x1)
meany <- mean(projCoords$coords.x2)
center <- (center - c(meanx,meany))/1000
projCoords$coords.x1 <- (projCoords$coords.x1 - meanx)/1000
projCoords$coords.x2 <- (projCoords$coords.x2 - meany)/1000

load("examples/ses/spatialCovs.RData")
bathdat <- getNOAA.bathy(-10, 155, -72, -35)

# time spent since left colony
time <- NULL
for(id in unique(projCoords$ID)) {
  nbSubObs <- length(which(projCoords$ID==id))
  
  # approximately in months for interval = 9.6h
  time <- c(time, (1:nbSubObs)/75)
}

projCoords$time <- time

# calculate distance/angle to center
data <- prepData(data=projCoords[,c(1:3,5)], coordNames = c("coords.x1","coords.x2"),covNames=c("time","ptime"),centers=center,spatialCovs=spatialCovs,gradient=TRUE,altCoordNames = "mu")

# plot correlation between steps and turns
allCor <- cylcop::plot_joint_box(x=data$step*1000,theta=data$angle)+labs(title="All individuals")
indCor <- list()
count <- 1
for(i in unique(data$ID)){
  indCor[[i]] <- plot_joint_box(x=data$step[which(data$ID==i)]*1000,theta=data$angle[which(data$ID==i)])+labs(title=paste("Individual",count))
  count <- count + 1
}
bottom_plots <- wrap_plots(indCor, ncol = 2)
pdf("examples/ses/sesStepTurnCor.pdf", height = 11, width = 8)
allCor / bottom_plots + plot_layout(heights = c(2.5, 5))
dev.off()

# include time since departure on 2->4 (column 6) and distance to colony and coast on 1->2 (column 1) in t.p.m. formula
formula <- ~betaCol1(colony.dist) + betaCol6(time)

# constrain transition probabilities
fixPar <- list(beta=matrix(c(NA,-1.e+10,-1.e+10,-1.e+10,NA,NA,-1.e+10,NA,-1.e+10,-1.e+10,-1.e+10,-1.e+10,
                             NA,   0,   0,   0, 0, 0,   0, 0,   0,   0,   0,   0,
                             0,   0,   0,   0, 0, NA,   0, 0,   0,   0,   0,   0),
                           nrow=3,byrow=TRUE),delta=matrix(-1.e+100,1,3))

stateNames <- c("outbound","state 2","state 3","inbound")
latexNames <- ifelse(stateNames %in% c("outbound", "inbound"), 
                     paste0("``", stateNames, "''"), 
                     stateNames)

initPos <- mapply(function(x) c(data[which(data$ID==x)[1],c("mu.x")],data[which(data$ID==x)[1],c("mu.y")]),unique(data$ID),SIMPLIFY = FALSE)

centroid <- xyFromCell(spatialCovs$d2coast,cellFromXY(spatialCovs$d2coast,center))

#############################################################################################
# Fit covariate-free models 
#############################################################################################

###############################
# standard step and turn model
###############################
m1 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), formula = formula, formulaDelta=~1, Par0=list(step=c(0.041,0.018,0.008,0.04,0.013,0.008,0.005,0.012),angle=c(10.262,3.448,1.e-3,9.807)),retryFits=0,retrySD=list(step=1,beta=1,delta=0),
             fixPar=fixPar, stateNames = stateNames)

plotPR(m1)

prm1 <- prPlot(m1,what=c("step","angle"))

plotTracks(m1$data %>% rename(x=mu.x,y=mu.y) %>% mutate(states=viterbi(m1)))

m1dens <- plotDens(m1)

# note simulation is blind to spatial process
set.seed(97,kind="Mersenne-Twister",normal.kind = "Inversion")
simm1 <- simMod(model=m1,initialPosition=initPos,centers=center,matchModelObs=TRUE)

pdf(file="examples/ses/m1_cor.pdf", width = 16, height = 22)
plotCor(data,m1,simm1,patchwork=TRUE)
dev.off()

###############################
# correlated step and turn model
###############################

crwm1 <- fitHMM(data=data, nbStates=4, dist=list(step="crwrice",angle="crwvm"), formula=formula, formulaDelta=~1, Par0=list(step=c(0.24,0.41,3.66,0.17,c(23.34,15.77,27.1,23.66)/1000),angle=c(0.24,0.41,3.66,0.17,c(23.34,15.77,27.1,23.66)/1000)),retryFits=0,retrySD=list(step=1,beta=1,delta=0),
                fixPar=fixPar, stateNames = stateNames)

plotPR(crwm1)

prcrwm1 <- prPlot(crwm1,what=c("step","angle"))

plotTracks(crwm1$data %>% rename(x=mu.x,y=mu.y) %>% mutate(states=viterbi(crwm1)))

crwm1dens <- plotDens(crwm1)

# note simulation is blind to spatial process
set.seed(97,kind="Mersenne-Twister",normal.kind = "Inversion")
simcrwm1 <- simMod(model=crwm1,initialPosition=initPos,centers=center,matchModelObs = TRUE)

pdf(file="examples/ses/crwm1_cor.pdf", width = 16, height = 22)
plotCor(data,crwm1,simcrwm1,patchwork=TRUE)
dev.off()

## AIC heavily favors correlated step and turn model
AICweights(crwm1,m1)

# agreement in decoded states
mean(viterbi(m1)==viterbi(crwm1))
cm <- caret::confusionMatrix(data=factor(viterbi(m1)),reference=factor(viterbi(crwm1)))
#cm_to_latex(
#  cm = cm,
#  caption = "Contingency table comparing the Viterbi-decoded state assignments for the 4-state standard step and turn model (``ST'') and the correlated step and turn model (``crwST'') in the southern elephant seal example.",
#  label = "tab:STstates",
#  pred_name = "ST",
#  ref_name = "crwST",
#  row_names = latexNames,
#  col_names = latexNames
#)

pdf("examples/ses/STdens.pdf",width=11,height=11)
m1dens$stepPlot + m1dens$anglePlot + crwm1dens$stepPlot + crwm1dens$anglePlot + plot_layout(2,2,byrow=TRUE,guides="collect",axes="collect")
dev.off()

pdf("examples/ses/STsteppr.pdf",width=11,height=11)
prm1$step$pr + prm1$step$qq + prm1$step$acf + prcrwm1$step$pr + prcrwm1$step$qq + prcrwm1$step$acf + plot_layout(axes="collect")
dev.off()

pdf("examples/ses/STanglepr.pdf",width=11,height=11)
prm1$angle$pr + prm1$angle$qq + prm1$angle$acf + prcrwm1$angle$pr + prcrwm1$angle$qq + prcrwm1$angle$acf + plot_layout(axes="collect")
dev.off()

###############################
# position-based correlated random walk model
###############################

DMp <- list(mu=list(mean.x=~crw(mu.x_tm1),
                    mean.y=~crw(mu.y_tm1),
                    sd.x=~ 1,
                    sd.y=~ 1,
                    corr.xy=~1))

fixParmp <- list(mu=c(NA,1,NA,2,NA,3,NA,4,NA,5,NA,6,NA,7,NA,8,9,10,11,12,9,10,11,12,rep(NA,4)),
                 beta=matrix(c(1,2,NA,rep(NA,9),3,NA,NA,4,NA,5,rep(NA,3),6,NA,NA,rep(NA,12)),nrow=3),delta=matrix(NA,1,3))

mp <- fitHMM(data=data, nbStates=4, dist=list(mu="rw_mvnorm2"), formula=formula, DM=DMp, Par0=list(mu=c(1,0.837,1,0.682,1,-0.354,1,0.885,1,0.903,1,0.7,1,-0.249,1,0.932,-4.376,-4.828,-5.243,-4.41,-4.376,-4.828,-5.243,-4.41,0,0,0,0)),beta0=getPar(m1)$beta, optMethod = "TMB",control=list(silent=TRUE,iter.max=10000,eval.max=10000,maxit=10000),retryFits=0,retrySD=list(step=0.1,beta=1,delta=0),
             fixPar=fixParmp, stateNames = stateNames, formulaDelta=~1, delta0=matrix(-1.e+100,1,3),mvnCoords = "mu")

plotPR(mp)

prmp <- prPlot(mp,what="mu")

plotTracks(mp$data %>% rename(x=mu.x_tm1,y=mu.y_tm1) %>% mutate(states=viterbi(mp)))

mpdens <- plotDens(mp)

# note simulation is blind to spatial process
set.seed(97,kind="Mersenne-Twister",normal.kind = "Inversion")
simmp <- simMod(model=mp,initialPosition=initPos,centers=center,matchModelObs = TRUE)

pdf(file="examples/ses/mp_cor.pdf", width = 16, height = 22)
plotCor(data,mp,simmp,patchwork=TRUE)
dev.off()

# agreement in decoded states much higher between crwm1 and mp
mean(viterbi(m1)==viterbi(mp))
mean(viterbi(crwm1)==viterbi(mp))

#############################################################################################
# biased correlated random walk model including colony covariates on step length and turning angle concentration parameters
#############################################################################################
distFormula <- ~ state1(colony.dist) + state4(colony.dist)
angleFormula <- ~ state1(colony.angle) + state4(colony.angle)

stepDM <- list(mean=distFormula, sd=distFormula)
angleDM <- list(mean=angleFormula, concentration=distFormula)

m2 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), Par0=list(step=c(-3.1,-0.07,-4,-4.86,-3.1,-0.07,-4.25,-0.07,-4.81,-5.3,-4.02,-0.34), angle=c(-1.1,0.45,2.18,0.2,1.24,-17.73,2.3,0.1)),
             beta0=getPar(m1)$beta, fixPar=fixPar, formula=formula, formulaDelta=~1, DM = list(step=stepDM, angle=angleDM), 
             estAngleMean=list(angle=TRUE), circularAngleMean=list(angle=TRUE), stateNames = stateNames)

AICweights(m1,crwm1,m2)

plotPR(m2)

prm2 <- prPlot(m2,what=c("step","angle"))

pdf("examples/ses/m2.pdf",width=12,height=6)
plotTracks(m2$data %>% rename(x=mu.x,y=mu.y) %>% mutate(states=viterbi(m2)))
dev.off()

m2dens <- plotDens(m2)

set.seed(97,kind="Mersenne-Twister",normal.kind = "Inversion")
simm2 <- simMod(model=m2,initialPosition=initPos,centers=center,spCovs=spatialCovs)

pdf(file=paste0("examples/ses/simm2.pdf"),width=12,height=6)
plotTracks(simm2)
dev.off()

pdf(file="examples/ses/m2_cor.pdf", width = 16, height = 22)
plotCor(data,m2,simm2,patchwork=TRUE)
dev.off()

bathm2 <- list()
bathm2$state2 <- plotHist(m2,state=2,title="State 2")
bathm2$state3 <- plotHist(m2,state=3,title="State 3")

#############################################################################################
# position-based potential function model
#############################################################################################

DMpos <- list(mu=list(mean.x=~state1(crw(mu.x_tm1))+state1(d2coast.x+dfrcol.x+land.x+boundary.x)+state2(crw(mu.x_tm1))+state2(d2coast.x+land.x)+state3(crw(mu.x_tm1))+state3(d2coast.x+land.x)+state4(crw(mu.x_tm1))+state4(d2col.x+d2coast.x+land.x),
                      mean.y=~state1(crw(mu.y_tm1))+state1(d2coast.y+dfrcol.y+land.y+boundary.y)+state2(crw(mu.y_tm1))+state2(d2coast.y+land.y)+state3(crw(mu.y_tm1))+state3(d2coast.y+land.y)+state4(crw(mu.y_tm1))+state4(d2col.y+d2coast.y+land.y),
                      sd.x=~ bath + d2coast + land + state1(colony.dist) + state4(colony.dist),
                      sd.y=~ bath + d2coast + land + state1(colony.dist) + state4(colony.dist),
                      corr.xy=~1))

formulapos <- ~betaCol1(bath+colony.dist+I(colony.dist^2)+I(colony.dist^3)+d2coast+I(d2coast^2)+I(d2coast^3)) + betaCol5(bath) + betaCol6(time) + betaCol8(bath)

par <- list(Par=list(mu=c(1,0.42,4.75,0.25,4.61,0.44,1,0.67,7.07,5.52,1,-0.34,39.43,6.91,1,0.51,3.9,0.01,4.61,1,0.61,4.45,0,4.61,-0.84,1,0.66,19.68,5.52,1,-0.16,8.25,6.91,1,0.55,4.62,0.01,4.61,-4.11,0.05,-0.09,5.66,-0.08,-4.42,0.14,0.2,3.66,-5.36,-0.07,-0.49,0.56,-4.23,0.07,0.01,2.87,-0.09,-4.11,0.05,-0.09,5.66,-0.08,-4.42,0.14,0.2,3.66,-5.36,-0.07,-0.49,0.56,-4.23,0.07,0.01,2.87,-0.09,0,0,0,0)),
            beta=matrix(c(-27.88,-0.67,39.9,-14.21,1.65,-89.89,136.32,-62.53,0,-1e+10,0,0,0,0,0,0,0,0,-1e+10,0,0,0,0,0,0,0,0,-1e+10,0,0,0,0,0,0,0,0,-0.9,0.75,0,0,0,0,0,0,0,-16.64,0,0,0,0,0,0,0,1.97,-1e+10,0,0,0,0,0,0,0,0,-2.95,-0.14,0,0,0,0,0,0,0,-1e+10,0,0,0,0,0,0,0,0,-1e+10,0,0,0,0,0,0,0,0,-1e+10,0,0,0,0,0,0,0,0,-1e+10,0,0,0,0,0,0,0,0),nrow=9))
workBounds <- list(mu=matrix(c(-Inf,Inf),nrow=length(par$Par$mu),ncol=2,byrow = TRUE))
workBounds$mu[c(3,5,9,10,13,14,17,19,22,24,28,29,32,33,36,38),2] <- 0
workBounds$mu[c(6,25),1] <- 0

fixParpos <- list(mu=c(NA,1:3,NA,4,NA,5:6,NA,NA,7:8,NA,NA,9:11,NA,NA,11+1:3,NA,15,NA,11+5:6,NA,NA,11+7:8,NA,NA,11+9:11,NA,11+12:29,11+12:29,rep(NA,4)),
                  beta=matrix(c(1,rep(NA,3),9,11,NA,13,rep(NA,4),
                                2,rep(NA,3),10,NA,NA,14,rep(NA,4),
                                3,rep(NA,11),
                                4,rep(NA,11),
                                5,rep(NA,11),
                                6,rep(NA,11),
                                7,rep(NA,11),
                                8,rep(NA,11),
                                rep(NA,5),12,rep(NA,6)),nrow=9,byrow=TRUE),
                  delta=matrix(NA,1,3))

mpos <- fitHMM(data=data, nbStates=4, dist=list(mu="rw_mvnorm2"), formula=formulapos, DM=DMpos, Par0=par$Par, beta0=par$beta,optMethod = "TMB",control=list(silent=TRUE,iter.max=10000,eval.max=10000,maxit=10000),retryFits=0,retrySD=list(mu=1,beta=1,delta=0),
               fixPar=fixParpos, stateNames = stateNames, formulaDelta=~1, delta0=matrix(-1.e+100,1,3),mvnCoords = "mu",workBounds=workBounds)

AICweights(mpos,mp)

muEst <- t(do.call(rbind,mpos$CIbeta$mu))[!is.na(fixParpos$mu) & !duplicated(fixParpos$mu),]
colnames(muEst) <- c("est","se","lower","upper")
muEst

betaEst <- do.call(cbind,lapply(mpos$CIbeta$beta,function(x) x[!is.na(fixParpos$beta)]))
bInd <- arrayInd(which(!is.na(fixParpos$beta)),dim(mpos$CIbeta$beta$est),.dimnames=dimnames(mpos$CIbeta$beta$est),useNames = TRUE)
rownames(betaEst) <- paste0(colnames(mpos$CIbeta$beta$est)[bInd[,"col"]]," : ",rownames(mpos$CIbeta$beta$est)[bInd[,"row"]])
#create_latex_table(muEst, betaEst)

plotPR(mpos)

prmpos <- prPlot(mpos,what="mu")

plmpos <- plot(mpos,plotCI=TRUE,ask=FALSE,return=TRUE)

df <- data.frame(do.call(rbind,plmpos$estimates$beta$bath[c("1 -> 2","2 -> 3","3 -> 2")]),gamma=rep(c("1 -> 2","2 -> 3","3 -> 2"),each=101),cols=rep(c("#78c679","#F0E442","#E31A1C"),each=101))
color_mapping <- setNames(unique(df$cols), unique(df$gamma))
pdf("examples/ses/gammaBath.pdf",width=11,height=8)
ggplot(df, aes(x = bath, group = gamma)) +
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = gamma), alpha = 0.5) +
  geom_line(aes(y = est, color = gamma), linewidth = 1.5, alpha = 0.5) +
  scale_fill_manual(name = "State transition", values = color_mapping) +
  scale_color_manual(name = "State transition", values = color_mapping) +
  labs(
    x = "depth (km)",
    y = "Probability"
  ) + theme(text = element_text(size = 15))
dev.off()

pdf("examples/ses/mpos.pdf",width=12,height=6)
plotTracks(mpos$data %>% rename(x=mu.x_tm1,y=mu.y_tm1) %>% mutate(states=viterbi(mpos)))
dev.off()

mposdens <- plotDens(mpos)

## agreement in decoded states
mean(viterbi(mpos)==viterbi(m2))
cmb <- caret::confusionMatrix(data=factor(viterbi(m2)),reference=factor(viterbi(mpos)))
#cm_to_latex(
#  cm = cmb,
#  caption = "Contingency table comparing the Viterbi-decoded state assignments for the 4-state biased random walk step and turn model (``STfull'') and the position-based potential function model (``pot'') in the southern elephant seal example.",
#  label = "tab:bcrwstates",
#  pred_name = "STfull",
#  ref_name = "pot",
#  row_names = latexNames,
#  col_names = latexNames
#)

timeInStates(mpos)
## do a bunch of simulations (takes a while)
#simout <- list()
#for(seed in 1:100){
#  set.seed(seed,kind="Mersenne-Twister",normal.kind = "Inversion")
#  simout[[seed]] <- simMod(mpos,initialPosition=initPos,centers=center,spCovs=spatialCovs,gradient=TRUE,retrySims=100)
#}
#apply(do.call(rbind,lapply(simout,function(x) table(x$states)/nrow(x))),2,mean)

bathmpos <- list()
bathmpos$state2 <- plotHist(mpos,state=2)
bathmpos$state3 <- plotHist(mpos,state=3)

set.seed(97,kind="Mersenne-Twister",normal.kind = "Inversion")
simpos <- simMod(mpos,initialPosition=initPos,centers=center,spCovs=spatialCovs,gradient=TRUE,retrySims=100)

pdf(file=paste0("examples/ses/simmpos.pdf"),width=12,height=6)
plotTracks(simpos)
dev.off()

pdf(file="examples/ses/mpos_cor.pdf", width = 16, height = 22)
plotCor(data,mpos,simpos,patchwork=TRUE)
dev.off()

m1dens$stepPlot + m1dens$anglePlot + crwm1dens$stepPlot + crwm1dens$anglePlot + mpdens$stepPlot + mpdens$anglePlot + m2dens$stepPlot + m2dens$anglePlot + mposdens$stepPlot + mposdens$anglePlot + 
  plot_layout(2,5,byrow=TRUE,axes="collect",guides = "collect",tag_level="new") + plot_annotation(tag_levels = "a")

pdf("examples/ses/bcrwDens.pdf",height=11,width=11)
((m2dens$stepPlot + m2dens$anglePlot + mposdens$stepPlot + mposdens$anglePlot)  &  theme(plot.tag.position = "topleft")) +
  plot_layout(2,2,byrow=TRUE,axes="collect",guides = "collect")#,tag_level="new") + plot_annotation(tag_levels = "a") 
dev.off()

pdf("examples/ses/bcrwpr.pdf",width=11,height=11)
prm2$step$pr + prm2$step$qq + prm2$step$acf +
  prm2$angle$pr + prm2$angle$qq + prm2$angle$acf + 
  prmpos$mu$pr + prmpos$mu$qq + prmpos$mu$acf + plot_layout(axes="collect")
dev.off()

pdf("examples/ses/bathComp.pdf",width=11, height=11)
bathm2$state2 + bathm2$state3 + bathmpos$state2 + bathmpos$state3 + plot_layout(2,2,axes = "collect")
dev.off()

save.image("examples/ses/sesExample.RData")
