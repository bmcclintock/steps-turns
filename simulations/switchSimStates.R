switchStates <- function(pars,idx){
  nam <- colnames(pars$est)
  newpars <- lapply(pars,function(x) x[,idx,drop=FALSE])
  lapply(newpars,function(x) {colnames(x) <- nam;return(x)})
}

for(sc in c("low","moderate","high")){
  load(paste0("spatialSimStudy_",sc,".RData"))
  for(isim in 1:nsims){
    for(mod in 2:4){
      if(results[[isim]][[mod]]$CIreal$step$est[1,1]>results[[isim]][[mod]]$CIreal$step$est[1,2]){
        if(results[[isim]][[mod]]$CIreal$angle$est[1+ifelse(mod==4,1,0),1]>results[[isim]][[mod]]$CIreal$angle$est[1+ifelse(mod==4,1,0),2]){
          if(mod==2){
            results[[isim]][[mod]]$CIbeta$step <- switchStates(results[[isim]][[mod]]$CIbeta$step,c(2,1,4,3))
            results[[isim]][[mod]]$CIbeta$angle <- switchStates(results[[isim]][[mod]]$CIbeta$angle,c(2,1))
            results[[isim]][[mod]]$CIbeta$beta <- lapply(results[[isim]][[mod]]$CIbeta$beta,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$CIreal$step <- lapply(results[[isim]][[mod]]$CIreal$step,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$CIreal$angle <- lapply(results[[isim]][[mod]]$CIreal$angle,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$viterbi <- (!results[[isim]][[mod]]$viterbi==2)+1
            results[[isim]][[mod]]$stateProbs[,c(1,2)] <- results[[isim]][[mod]]$stateProbs[,c(2,1)]
          } else if(mod==3){
            results[[isim]][[mod]]$CIbeta$step <- switchStates(results[[isim]][[mod]]$CIbeta$step,c(6:10,1:5,12,11))
            results[[isim]][[mod]]$CIbeta$angle <- switchStates(results[[isim]][[mod]]$CIbeta$angle,c(6:10,1:5))
            results[[isim]][[mod]]$CIbeta$beta <- lapply(results[[isim]][[mod]]$CIbeta$beta,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$CIreal$step <- lapply(results[[isim]][[mod]]$CIreal$step,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$CIreal$angle <- lapply(results[[isim]][[mod]]$CIreal$angle,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$viterbi <- (!results[[isim]][[mod]]$viterbi==2)+1
            results[[isim]][[mod]]$stateProbs[,c(1,2)] <- results[[isim]][[mod]]$stateProbs[,c(2,1)]
          } else {
            results[[isim]][[mod]]$CIbeta$step <- switchStates(results[[isim]][[mod]]$CIbeta$step,c(6:10,1:5,12,11))
            results[[isim]][[mod]]$CIbeta$angle <- switchStates(results[[isim]][[mod]]$CIbeta$angle,c(6:10,1:5))
            results[[isim]][[mod]]$CIbeta$beta <- lapply(results[[isim]][[mod]]$CIbeta$beta,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$CIreal$step <- lapply(results[[isim]][[mod]]$CIreal$step,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$CIreal$angle <- lapply(results[[isim]][[mod]]$CIreal$angle,function(x) x[,c(2,1)])
            results[[isim]][[mod]]$viterbi <- (!results[[isim]][[mod]]$viterbi==2)+1
            results[[isim]][[mod]]$stateProbs[,c(1,2)] <- results[[isim]][[mod]]$stateProbs[,c(2,1)]
          }
        }
      }
    }
  }
  save(beta0,beta1,beta2,crw1,crw2,results,sd_1,sd_2,file=paste0("spatialSimStudy_",sc,".RData"))
}