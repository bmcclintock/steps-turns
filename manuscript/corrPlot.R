beta <- seq(-0.99,0.99,length=21)

x<-matrix(0,nrow=100001,ncol=2)
for(i in 1:100000){
  x[i+1,] <- mvtnorm::rmvnorm(1,x[i,],sigma=matrix(c(1,0.5,0.5,1),2,2))                     
}
dat <- prepData(data.frame(ID=1,x=x[,1],y=x[,2]))
cor(dat$step,abs(dat$angle),use="complete")

corr2 <- data.frame(beta=beta,cor=NA)
for(k in 1:length(beta)){
  x<-matrix(0,nrow=100001,ncol=2)
  x[2,] <- mvtnorm::rmvnorm(1,x[1,],sigma=matrix(c(1,0.5,0.5,1),2,2))     
  for(i in 2:100000){
    x[i+1,] <- mvtnorm::rmvnorm(1,x[i,]+(x[i,]-x[i-1,])*beta[k],sigma=matrix(c(1,0.5,0.5,1),2,2))                     
  }
  dat <- momentuHMM::prepData(data.frame(ID=1,x=x[,1],y=x[,2]))
  corr2[k,"cor"] <- cor(dat$step,abs(dat$angle),use="complete")
}

library(ggplot2)
ggplot(corr2) + geom_line(aes(x=beta,y=cor))+xlab(expression(beta))+ylab("correlation between step length and turn angle")+theme_gray(base_size=14)

corr3 <- data.frame(beta=beta,cor=NA)
for(k in 1:length(beta)){
  x<-matrix(0,nrow=100001,ncol=2)
  x[2,] <- mvtnorm::rmvnorm(1,x[1,],sigma=matrix(c(1,0,0,1),2,2))     
  for(i in 2:100000){
    x[i+1,] <- mvtnorm::rmvnorm(1,x[i,]+(x[i,]-x[i-1,])*beta[k],sigma=matrix(c(1,0,0,1),2,2))                     
  }
  dat <- momentuHMM::prepData(data.frame(ID=1,x=x[,1],y=x[,2]))
  corr3[k,"cor"] <- cor(dat$step,abs(dat$angle),use="complete")
}

ggplot(corr3) + geom_line(aes(x=beta,y=cor))+xlab(expression(beta))+ylab("correlation between step length and turn angle")+theme_gray(base_size=14)

corr4 <- data.frame(beta=beta,cor=NA)
for(k in 1:length(beta)){
  x<-matrix(0,nrow=100001,ncol=2)
  x[2,] <- mvtnorm::rmvnorm(1,x[1,],sigma=matrix(c(1,-0.5,-0.5,1),2,2))     
  for(i in 2:100000){
    x[i+1,] <- mvtnorm::rmvnorm(1,x[i,]+(x[i,]-x[i-1,])*beta[k],sigma=matrix(c(1,-0.5,-0.5,1),2,2))                     
  }
  dat <- momentuHMM::prepData(data.frame(ID=1,x=x[,1],y=x[,2]))
  corr4[k,"cor"] <- cor(dat$step,abs(dat$angle),use="complete")
}

ggplot(corr4) + geom_line(aes(x=beta,y=cor))+xlab(expression(beta))+ylab("correlation between step length and turn angle")+theme_gray(base_size=14)

combcorr <- rbind(data.frame(corr4,rho="-0.5"),data.frame(corr3,rho="0"),data.frame(corr2,rho="0.5"))

pdf("corrPlot.pdf",width=8,height=8)
ggplot(combcorr) +geom_vline(xintercept = 0,linetype=2)+geom_hline(yintercept = 0,linetype=2)+ geom_line(aes(x=beta,y=cor,col=rho),linewidth=1.5,alpha=.75)+
  xlab(expression(beta))+ylab("correlation between step length and turn angle")+theme_gray(base_size=14)+guides(color = guide_legend(title = expression(rho)))
dev.off()
