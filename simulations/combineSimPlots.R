library(ggplot2)
library(patchwork)

load("simulations/crwplots.RData")
load("simulations/crwSTplots.RData")
load("simulations/spatialSimStudyplots.RData")
load("simulations/STplots.RData")
load("simulations/spatialSimStudyplots_noLangevin.RData")
load("simulations/crwplotDens.RData")
load("simulations/crwSTplotDens.RData")
load("simulations/spatialStudyplotDens.RData")
load("simulations/STplotDens.RData")
load("simulations/spatialStudyplotDens_noLangevin.RData")

base_size <- 14

pdf("simulations/crwViterbiPlots.pdf",width=8,height=11)
crwplot$viterbi$low +theme_gray(base_size = base_size)+ labs(title="crwpos",subtitle="low") + crwplot$viterbi$moderate +theme_gray(base_size = base_size)+labs(subtitle="moderate")+ crwplot$viterbi$high +theme_gray(base_size = base_size)+labs(subtitle="high")+
  crwSTplot$viterbi$low +theme_gray(base_size = base_size)+ ggtitle("crwST")+crwSTplot$viterbi$moderate +theme_gray(base_size = base_size)+ crwSTplot$viterbi$high +theme_gray(base_size = base_size)+ plot_layout(2,3,axes="collect",byrow=FALSE)
dev.off()

pdf("simulations/crwBeta12Plots.pdf",width=8,height=11)
crwplot$beta12$low +theme_gray(base_size = base_size)+ labs(title="crwpos",subtitle="low")+ylab(expression(beta[paste(1,",",2,",",1)]))+ crwplot$beta12$moderate +theme_gray(base_size = base_size)+labs(subtitle="moderate")+ylab(expression(beta[paste(1,",",2,",",1)]))+ crwplot$beta12$high +theme_gray(base_size = base_size)+labs(subtitle="high")+ylab(expression(beta[paste(1,",",2,",",1)]))+
  crwSTplot$beta12$low +theme_gray(base_size = base_size)+ ggtitle("crwST")+ylab(expression(beta[paste(1,",",2,",",1)]))+crwSTplot$beta12$moderate +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",1)]))+ crwSTplot$beta12$high +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",1)]))+ plot_layout(2,3,axes="collect",byrow=FALSE)
dev.off()

pdf("simulations/crwBeta21Plots.pdf",width=8,height=11)
crwplot$beta21$low +theme_gray(base_size = base_size)+ labs(title="crwpos",subtitle="low")+ylab(expression(beta[paste(2,",",1,",",1)]))+ crwplot$beta21$moderate +theme_gray(base_size = base_size)+labs(subtitle="moderate")+ylab(expression(beta[paste(2,",",1,",",1)]))+ crwplot$beta21$high +theme_gray(base_size = base_size)+labs(subtitle="high")+ylab(expression(beta[paste(2,",",1,",",1)]))+
  crwSTplot$beta21$low +theme_gray(base_size = base_size)+ ggtitle("crwST")+ylab(expression(beta[paste(2,",",1,",",1)]))+crwSTplot$beta21$moderate +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",1)]))+ crwSTplot$beta21$high +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",1)]))+ plot_layout(2,3,axes="collect",byrow=FALSE)
dev.off()

pdf("simulations/spatViterbiPlots.pdf",width=8,height=11)
spatialplot$viterbi$low +theme_gray(base_size = base_size)+ ggtitle("low")+spatialplot$viterbi$moderate +theme_gray(base_size = base_size)+ ggtitle("moderate")+ spatialplot$viterbi$high +theme_gray(base_size = base_size)+ ggtitle("high")+ plot_layout(1,3,axes="collect")
dev.off()

pdf("simulations/spatBeta12Plots.pdf",width=8,height=11)
spatialplot$beta12$low$cov1 +theme_gray(base_size = base_size)+ labs(title="low")+ylab(expression(beta[paste(1,",",2,",",1)]))+spatialplot$beta12$moderate$cov1 +theme_gray(base_size = base_size)+ ggtitle("moderate")+ spatialplot$beta12$high$cov1 +theme_gray(base_size = base_size)+ ggtitle("high")+
  spatialplot$beta12$low$cov2 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",2)]))+spatialplot$beta12$moderate$cov2 +theme_gray(base_size = base_size)+ spatialplot$beta12$high$cov2 +theme_gray(base_size = base_size)+ 
  spatialplot$beta12$low$cov3 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",3)]))+spatialplot$beta12$moderate$cov3 +theme_gray(base_size = base_size)+ spatialplot$beta12$high$cov2 +theme_gray(base_size = base_size)+ 
  spatialplot$beta12$low$d2c +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",4)]))+spatialplot$beta12$moderate$d2c +theme_gray(base_size = base_size)+ spatialplot$beta12$high$d2c +theme_gray(base_size = base_size)+ plot_layout(3,4,axes="collect",byrow=TRUE)
dev.off()

pdf("simulations/spatBeta21Plots.pdf",width=8,height=11)
spatialplot$beta21$low$cov1 +theme_gray(base_size = base_size)+ labs(title="low")+ylab(expression(beta[paste(2,",",1,",",1)]))+spatialplot$beta21$moderate$cov1 +theme_gray(base_size = base_size)+ ggtitle("moderate")+ spatialplot$beta21$high$cov1 +theme_gray(base_size = base_size)+ ggtitle("high")+
  spatialplot$beta21$low$cov2 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",2)]))+spatialplot$beta21$moderate$cov2 +theme_gray(base_size = base_size)+ spatialplot$beta21$high$cov2 +theme_gray(base_size = base_size)+ 
  spatialplot$beta21$low$cov3 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",3)]))+spatialplot$beta21$moderate$cov3 +theme_gray(base_size = base_size)+ spatialplot$beta21$high$cov2 +theme_gray(base_size = base_size)+ 
  spatialplot$beta21$low$d2c +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",4)]))+spatialplot$beta21$moderate$d2c +theme_gray(base_size = base_size)+ spatialplot$beta21$high$d2c +theme_gray(base_size = base_size)+ plot_layout(3,4,axes="collect",byrow=TRUE)
dev.off()

pdf("simulations/STViterbiPlots.pdf",width=11,height=4)
STplot$viterbi$low +theme_gray(base_size = base_size)+ labs(subtitle="low") + STplot$viterbi$moderate +theme_gray(base_size = base_size)+labs(subtitle="moderate")+ STplot$viterbi$high +theme_gray(base_size = base_size)+labs(subtitle="high")+
 plot_layout(3,1,axes="collect")
dev.off()

pdf("simulations/STBetaPlots.pdf",width=11,height=8)
STplot$beta12$low +theme_gray(base_size = base_size)+ labs(subtitle="low") +STplot$beta12$moderate +theme_gray(base_size = base_size)+labs(subtitle="moderate") +STplot$beta12$high +theme_gray(base_size = base_size)+labs(subtitle="high") +
  STplot$beta21$low +theme_gray(base_size = base_size)+ STplot$beta21$moderate +theme_gray(base_size = base_size)+STplot$beta21$high +theme_gray(base_size = base_size)+plot_layout(3,2,axes="collect",byrow=TRUE)
dev.off()

pdf("simulations/spat_noLangevinViterbiPlots.pdf",width=8,height=11)
spatial_noLangevinplot$viterbi$low +theme_gray(base_size = base_size)+ ggtitle("low")+spatial_noLangevinplot$viterbi$moderate +theme_gray(base_size = base_size)+ ggtitle("moderate")+ spatial_noLangevinplot$viterbi$high +theme_gray(base_size = base_size)+ ggtitle("high")+ plot_layout(1,3,axes="collect")
dev.off()

pdf("simulations/spat_noLangevinBeta12Plots.pdf",width=8,height=11)
spatial_noLangevinplot$beta12$low$cov1 +theme_gray(base_size = base_size)+ labs(title="low")+ylab(expression(beta[paste(1,",",2,",",1)]))+spatial_noLangevinplot$beta12$moderate$cov1 +theme_gray(base_size = base_size)+ ggtitle("moderate")+ spatial_noLangevinplot$beta12$high$cov1 +theme_gray(base_size = base_size)+ ggtitle("high")+
  spatial_noLangevinplot$beta12$low$cov2 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",2)]))+spatial_noLangevinplot$beta12$moderate$cov2 +theme_gray(base_size = base_size)+ spatial_noLangevinplot$beta12$high$cov2 +theme_gray(base_size = base_size)+ 
  spatial_noLangevinplot$beta12$low$cov3 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",3)]))+spatial_noLangevinplot$beta12$moderate$cov3 +theme_gray(base_size = base_size)+ spatial_noLangevinplot$beta12$high$cov2 +theme_gray(base_size = base_size)+ 
  spatial_noLangevinplot$beta12$low$d2c +theme_gray(base_size = base_size)+ylab(expression(beta[paste(1,",",2,",",4)]))+spatial_noLangevinplot$beta12$moderate$d2c +theme_gray(base_size = base_size)+ spatial_noLangevinplot$beta12$high$d2c +theme_gray(base_size = base_size)+ plot_layout(3,4,axes="collect",byrow=TRUE)
dev.off()

pdf("simulations/spat_noLangevinBeta21Plots.pdf",width=8,height=11)
spatial_noLangevinplot$beta21$low$cov1 +theme_gray(base_size = base_size)+ labs(title="low")+ylab(expression(beta[paste(2,",",1,",",1)]))+spatial_noLangevinplot$beta21$moderate$cov1 +theme_gray(base_size = base_size)+ ggtitle("moderate")+ spatial_noLangevinplot$beta21$high$cov1 +theme_gray(base_size = base_size)+ ggtitle("high")+
  spatial_noLangevinplot$beta21$low$cov2 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",2)]))+spatial_noLangevinplot$beta21$moderate$cov2 +theme_gray(base_size = base_size)+ spatial_noLangevinplot$beta21$high$cov2 +theme_gray(base_size = base_size)+ 
  spatial_noLangevinplot$beta21$low$cov3 +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",3)]))+spatial_noLangevinplot$beta21$moderate$cov3 +theme_gray(base_size = base_size)+ spatial_noLangevinplot$beta21$high$cov2 +theme_gray(base_size = base_size)+ 
  spatial_noLangevinplot$beta21$low$d2c +theme_gray(base_size = base_size)+ylab(expression(beta[paste(2,",",1,",",4)]))+spatial_noLangevinplot$beta21$moderate$d2c +theme_gray(base_size = base_size)+ spatial_noLangevinplot$beta21$high$d2c +theme_gray(base_size = base_size)+ plot_layout(3,4,axes="collect",byrow=TRUE)
dev.off()

pdf("simulations/crwDensPlots.pdf",height=11,width=8)
crwplotDens$low[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("low") + crwplotDens$low[[1]]$anglePlot +theme_gray(base_size = base_size)+
  crwplotDens$moderate[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("moderate")+ crwplotDens$moderate[[1]]$anglePlot +theme_gray(base_size = base_size)+
  crwplotDens$high[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("high")+ crwplotDens$high[[1]]$anglePlot +theme_gray(base_size = base_size)+ plot_layout(2,3,axes="collect",byrow=TRUE,guides="collect")
dev.off()

pdf("simulations/crwSTDensPlots.pdf",height=11,width=8)
crwSTplotDens$low[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("low") + crwSTplotDens$low[[1]]$anglePlot +theme_gray(base_size = base_size)+
  crwSTplotDens$moderate[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("moderate")+ crwSTplotDens$moderate[[1]]$anglePlot +theme_gray(base_size = base_size)+
  crwSTplotDens$high[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("high")+ crwSTplotDens$high[[1]]$anglePlot +theme_gray(base_size = base_size)+ plot_layout(2,3,axes="collect",byrow=TRUE,guides="collect")
dev.off()

pdf("simulations/spatialStudyplotDens.pdf",height=11,width=8)
spatialStudyplotDens$low$stepPlot +theme_gray(base_size = base_size)+ xlim(0,15)+ggtitle("low") + spatialStudyplotDens$low$anglePlot +theme_gray(base_size = base_size)+
  spatialStudyplotDens$moderate$stepPlot +theme_gray(base_size = base_size)+ xlim(0,15)+ ggtitle("moderate")+ spatialStudyplotDens$moderate$anglePlot+theme_gray(base_size = base_size)+
  spatialStudyplotDens$high$stepPlot +theme_gray(base_size = base_size)+ xlim(0,15)+ ggtitle("high")+ spatialStudyplotDens$high$anglePlot +theme_gray(base_size = base_size)+ plot_layout(2,3,axes="collect",byrow=TRUE,guides="collect")
dev.off()

pdf("simulations/STDensPlots.pdf",height=11,width=8)
STplotDens$low[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("low") + STplotDens$low[[1]]$anglePlot +theme_gray(base_size = base_size)+
  STplotDens$moderate[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("moderate")+ STplotDens$moderate[[1]]$anglePlot +theme_gray(base_size = base_size)+
  STplotDens$high[[1]]$stepPlot +theme_gray(base_size = base_size)+ ggtitle("high")+ STplotDens$high[[1]]$anglePlot +theme_gray(base_size = base_size)+ plot_layout(2,3,axes="collect",byrow=TRUE,guides="collect")
dev.off()

pdf("simulations/spatialStudyplotDens_noLangevin.pdf",height=11,width=8)
spatialStudyplotDens_noLangevin$low$stepPlot +theme_gray(base_size = base_size)+ xlim(0,15)+ ggtitle("low") + spatialStudyplotDens_noLangevin$low$anglePlot +theme_gray(base_size = base_size)+
  spatialStudyplotDens_noLangevin$moderate$stepPlot +theme_gray(base_size = base_size)+ xlim(0,15)+ ggtitle("moderate")+ spatialStudyplotDens_noLangevin$moderate$anglePlot +theme_gray(base_size = base_size)+
  spatialStudyplotDens_noLangevin$high$stepPlot +theme_gray(base_size = base_size)+ xlim(0,15)+ ggtitle("high")+ spatialStudyplotDens_noLangevin$high$anglePlot +theme_gray(base_size = base_size)+ plot_layout(2,3,axes="collect",byrow=TRUE,guides="collect")
dev.off()
