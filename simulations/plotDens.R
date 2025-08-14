plotDens <- function(plotDat,scenario){
  
  library(ggplot2)

  stepdens1 <- density(plotDat$step[which(plotDat$states==1)],na.rm=TRUE,adjust=2,from=0)
  stepdens2 <- density(plotDat$step[which(plotDat$states==2)],na.rm=TRUE,adjust=2,from=0)
  stepdens1$y <- stepdens1$y * mean(plotDat$states==1)
  stepdens2$y <- stepdens2$y * mean(plotDat$states==2)
  stepdf <- rbind(data.frame(x=stepdens1$x,y=stepdens1$y,state="state 1",col="#E69F00"),data.frame(x=stepdens2$x,y=stepdens2$y,state="state 2",col="#56B4E9"))
  stepdf$state <- factor(stepdf$state)
  COLORS <- c("state 1" = "#E69F00", "state 2" = "#56B4E9")
  
  
  stepPlot <- ggplot(plotDat) + geom_histogram(aes(x=step,y=after_stat(density)),binwidth=0.25,fill = "lightgrey", color = "black",na.rm=TRUE) +
    geom_line(data = stepdf[which(stepdf$state=="state 1"),], aes(x = x, y = y, color="state 1"),linewidth=1.5)+
    geom_line(data = stepdf[which(stepdf$state=="state 2"),], aes(x = x, y = y, color="state 2"),linewidth=1.5)+
    scale_color_manual(values = COLORS, name = "")+xlab("step length")
  
  
  angledens1 <- density(plotDat$angle[which(plotDat$states==1)],na.rm=TRUE,from=-pi,to=pi,adjust=2)
  angledens2 <- density(plotDat$angle[which(plotDat$states==2)],na.rm=TRUE,from=-pi,to=pi,adjust=2)
  angledens1$y <- angledens1$y * mean(plotDat$states==1)
  angledens2$y <- angledens2$y * mean(plotDat$states==2)
  
  angledf <- rbind(data.frame(x=angledens1$x,y=angledens1$y,state="state 1",col="#E69F00"),data.frame(x=angledens2$x,y=angledens2$y,state="state 2",col="#56B4E9"))
  angledf$state <- factor(angledf$state)
  
  anglePlot <- ggplot(plotDat) + geom_histogram(aes(x=angle,y=after_stat(density)),breaks=seq(-pi,pi,length=30),fill = "lightgrey", color = "black",na.rm=TRUE) +
    geom_line(data = angledf[which(angledf$state=="state 1"),], aes(x = x, y = y, color="state 1"),linewidth=1.5)+
    geom_line(data = angledf[which(angledf$state=="state 2"),], aes(x = x, y = y, color="state 2"),linewidth=1.5)+
    scale_color_manual(values = COLORS, name = "")+scale_x_continuous(breaks=c(-pi, -pi/2, 0, pi/2, pi),labels=expression(-pi, -pi/2, 0, pi/2, pi))+xlab("turn angle (radians)")
  
  return(list(stepPlot=stepPlot,anglePlot=anglePlot))
}

plotDens3 <- function(plotDat,scenario){
  
  library(ggplot2)
  
  nbStates <- length(unique(plotDat$states))
  
  stepdens1 <- density(plotDat$step[which(plotDat$states==1)],na.rm=TRUE,adjust=2,from=0,to=max(plotDat$step,na.rm=TRUE))
  stepdens2 <- density(plotDat$step[which(plotDat$states==2)],na.rm=TRUE,adjust=2,from=0,to=max(plotDat$step,na.rm=TRUE))
  stepdens1$y <- stepdens1$y * mean(plotDat$states==1)
  stepdens2$y <- stepdens2$y * mean(plotDat$states==2)
  stepdf <- rbind(data.frame(x=stepdens1$x,y=stepdens1$y,state="state 1",col="#E69F00"),data.frame(x=stepdens2$x,y=stepdens2$y,state="state 2",col="#56B4E9"))
  COLORS <- c("state 1" = "#E69F00", "state 2" = "#56B4E9")
  tot <- stepdens1$y + stepdens2$y
  linetype <- c("solid","solid","dashed")
  
  if(nbStates==3){
    stepdens3 <- density(plotDat$step[which(plotDat$states==3)],na.rm=TRUE,adjust=2,from=0,to=max(plotDat$step,na.rm=TRUE))
    stepdens3$y <- stepdens3$y * mean(plotDat$states==3)
    stepdf <- rbind(stepdf,data.frame(x=stepdens3$x,y=stepdens3$y,state="state 3",col="#009E73"))
    COLORS <- c(COLORS,"state 3" = "#009E73")
    tot <- stepdens1$y + stepdens2$y + stepdens3$y
    linetype <- c("solid","solid","solid","dashed")
  }
  
  stepdf <- rbind(stepdf,data.frame(x=stepdens1$x,y=tot,state="total",col="black"))
  COLORS <- c(COLORS,"total" = "black")
  
  stepdf$state <- factor(stepdf$state)
  
  stepPlot <- ggplot(plotDat) + geom_histogram(aes(x=step,y=after_stat(density)),binwidth=25,fill = "lightgrey", color = "black",na.rm=TRUE,boundary=0) +
    geom_line(data = stepdf[which(stepdf$state=="state 1"),], aes(x = x, y = y, color="state 1",linetype="state 1"),linewidth=1)+
    geom_line(data = stepdf[which(stepdf$state=="state 2"),], aes(x = x, y = y, color="state 2",linetype="state 2"),linewidth=1)+
    {if(nbStates==3) geom_line(data = stepdf[which(stepdf$state=="state 3"),], aes(x = x, y = y, color="state 3",linetype="state 3"),linewidth=1)}+
    geom_line(data = stepdf[which(stepdf$state=="total"),], aes(x = x, y = y, color="total",linetype="total"),linewidth=1)+
    scale_linetype_manual(values=linetype) + guides(color=guide_legend("linetype"), fill = "none")+
    theme_gray(base_size=14) +scale_color_manual(values = COLORS, name = "")+xlab("step length")+
    theme(legend.title=element_blank())

  
  
  angledens1 <- density(plotDat$angle[which(plotDat$states==1)],na.rm=TRUE,adjust=2,from=-pi,to=pi)
  angledens2 <- density(plotDat$angle[which(plotDat$states==2)],na.rm=TRUE,adjust=2,from=-pi,to=pi)
  angledens1$y <- angledens1$y * mean(plotDat$states==1)
  angledens2$y <- angledens2$y * mean(plotDat$states==2)
  angledf <- rbind(data.frame(x=angledens1$x,y=angledens1$y,state="state 1",col="#E69F00"),data.frame(x=angledens2$x,y=angledens2$y,state="state 2",col="#56B4E9"))
  tot <- angledens1$y + angledens2$y
  
  if(nbStates==3){
    angledens3 <- density(plotDat$angle[which(plotDat$states==3)],na.rm=TRUE,adjust=2,from=-pi,to=pi)
    angledens3$y <- angledens3$y * mean(plotDat$states==3)
    angledf <- rbind(angledf,data.frame(x=angledens3$x,y=angledens3$y,state="state 3",col="#009E73"))
    tot <- angledens1$y + angledens2$y + angledens3$y
  }
  
  angledf <- rbind(angledf,data.frame(x=angledens1$x,y=tot,state="total",col="black"))
  
  angledf$state <- factor(angledf$state)
  
  anglePlot <- ggplot(plotDat) + geom_histogram(aes(x=angle,y=after_stat(density)),binwidth=.25,fill = "lightgrey", color = "black",na.rm=TRUE,boundary=0) +
    geom_line(data = angledf[which(angledf$state=="state 1"),], aes(x = x, y = y, color="state 1",linetype="state 1"),linewidth=1)+
    geom_line(data = angledf[which(angledf$state=="state 2"),], aes(x = x, y = y, color="state 2",linetype="state 2"),linewidth=1)+
    {if(nbStates==3) geom_line(data = angledf[which(angledf$state=="state 3"),], aes(x = x, y = y, color="state 3",linetype="state 3"),linewidth=1)}+
    geom_line(data = angledf[which(stepdf$state=="total"),], aes(x = x, y = y, color="total",linetype="total"),linewidth=1)+
    scale_linetype_manual(values=linetype) + guides(color=guide_legend("linetype"), fill = "none")+
    theme_gray(base_size=14) +scale_color_manual(values = COLORS, name = "")+scale_x_continuous(limits=c(-pi,pi),breaks=c(-pi, -pi/2, 0, pi/2, pi),labels=expression(-pi, -pi/2, 0, pi/2, pi))+xlab("turn angle (radians)")+
    theme(legend.title=element_blank())
  
  return(list(stepPlot=stepPlot,anglePlot=anglePlot))
}


