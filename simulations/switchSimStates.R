switchStates <- function(pars,idx){
  nam <- colnames(pars$est)
  newpars <- lapply(pars,function(x) x[,idx,drop=FALSE])
  lapply(newpars,function(x) {colnames(x) <- nam;return(x)})
}

for(sc in c("low","moderate","high")){
  load(paste0("spatialSimStudy_",sc,".RData"))
  for(isim in 1:nsims){
    for(mod in 2:4){
      if(results_list[[isim]]$fit_results[[mod]]$CIreal$step$est[1,1]>results_list[[isim]]$fit_results[[mod]]$CIreal$step$est[1,2]){
        if(results_list[[isim]]$fit_results[[mod]]$CIreal$angle$est[1+ifelse(mod==4,1,0),1]>results_list[[isim]]$fit_results[[mod]]$CIreal$angle$est[1+ifelse(mod==4,1,0),2]){
          if(mod==2){
            results_list[[isim]]$fit_results[[mod]]$CIbeta$step <- switchStates(results_list[[isim]]$fit_results[[mod]]$CIbeta$step,c(2,1,4,3))
            results_list[[isim]]$fit_results[[mod]]$CIbeta$angle <- switchStates(results_list[[isim]]$fit_results[[mod]]$CIbeta$angle,c(2,1))
            results_list[[isim]]$fit_results[[mod]]$CIbeta$beta <- lapply(results_list[[isim]]$fit_results[[mod]]$CIbeta$beta,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$CIreal$step <- lapply(results_list[[isim]]$fit_results[[mod]]$CIreal$step,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$CIreal$angle <- lapply(results_list[[isim]]$fit_results[[mod]]$CIreal$angle,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$viterbi <- (!results_list[[isim]]$fit_results[[mod]]$viterbi==2)+1
            results_list[[isim]]$fit_results[[mod]]$stateProbs[,c(1,2)] <- results_list[[isim]]$fit_results[[mod]]$stateProbs[,c(2,1)]
          } else if(mod==3){
            results_list[[isim]]$fit_results[[mod]]$CIbeta$step <- switchStates(results_list[[isim]]$fit_results[[mod]]$CIbeta$step,c(6:10,1:5,12,11))
            results_list[[isim]]$fit_results[[mod]]$CIbeta$angle <- switchStates(results_list[[isim]]$fit_results[[mod]]$CIbeta$angle,c(6:10,1:5))
            results_list[[isim]]$fit_results[[mod]]$CIbeta$beta <- lapply(results_list[[isim]]$fit_results[[mod]]$CIbeta$beta,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$CIreal$step <- lapply(results_list[[isim]]$fit_results[[mod]]$CIreal$step,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$CIreal$angle <- lapply(results_list[[isim]]$fit_results[[mod]]$CIreal$angle,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$viterbi <- (!results_list[[isim]]$fit_results[[mod]]$viterbi==2)+1
            results_list[[isim]]$fit_results[[mod]]$stateProbs[,c(1,2)] <- results_list[[isim]]$fit_results[[mod]]$stateProbs[,c(2,1)]
          } else {
            results_list[[isim]]$fit_results[[mod]]$CIbeta$step <- switchStates(results_list[[isim]]$fit_results[[mod]]$CIbeta$step,c(6:10,1:5,12,11))
            results_list[[isim]]$fit_results[[mod]]$CIbeta$angle <- switchStates(results_list[[isim]]$fit_results[[mod]]$CIbeta$angle,c(6:10,1:5))
            results_list[[isim]]$fit_results[[mod]]$CIbeta$beta <- lapply(results_list[[isim]]$fit_results[[mod]]$CIbeta$beta,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$CIreal$step <- lapply(results_list[[isim]]$fit_results[[mod]]$CIreal$step,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$CIreal$angle <- lapply(results_list[[isim]]$fit_results[[mod]]$CIreal$angle,function(x) x[,c(2,1)])
            results_list[[isim]]$fit_results[[mod]]$viterbi <- (!results_list[[isim]]$fit_results[[mod]]$viterbi==2)+1
            results_list[[isim]]$fit_results[[mod]]$stateProbs[,c(1,2)] <- results_list[[isim]]$fit_results[[mod]]$stateProbs[,c(2,1)]
          }
        }
      }
    }
  }
  save(beta0,delta1,delta2,beta1,beta2,results_list,sd_1,sd_2,file=paste0("simulations/spatialSimStudy_",sc,".RData"))
}