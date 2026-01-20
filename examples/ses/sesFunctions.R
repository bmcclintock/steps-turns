simMod <- function(model,obsPerAnimal=1000,nbAnimals=6,covs=data.frame(time=(1:1000)/75),initialPosition,centers,spCovs=NULL,gradient=FALSE,states=TRUE,retrySims=10,matchModelObs=FALSE,plot=TRUE,mask=NULL){

  if(is.null(spCovs)) attr(model$data,"gradient") <- FALSE
  
  sim <- momentuHMM::simData(model=model,obsPerAnimal=obsPerAnimal,nbAnimals=nbAnimals,covs=covs,initialPosition=initialPosition,centers=centers,spatialCovs=spCovs,gradient=gradient,states=states,retrySims=retrySims,matchModelObs = matchModelObs,mask=mask)
  
  if(!is.null(model$conditions$mvnCoords)){
    sim <- prepData(sim,coordNames = c("mu.x","mu.y")) # get steps and turns
  }

  remInd <- NULL
  for(i in 1:length(unique(sim$ID))){
    colInd <- which((sim$ID==unique(sim$ID)[i] & sim$states==4) & ((sqrt((sim$x-center[1])^2+(sim$y-center[2])^2)<0.1)))
    if(length(colInd) && (colInd[1]+1)<(obsPerAnimal*i)) remInd <- c(remInd,(colInd[1]+1):(obsPerAnimal*i)) # truncate tracks once back to colony
  }
  if(length(remInd)) sim <- sim[-remInd,]
  names(spatialCovs$land) <- "land"
  if(plot) momentuHMM::plotSpatialCov(sim,spatialCovs$land,states=sim$states,stateNames=model$stateNames,color=viridis::viridis(10))
  return(sim)
}

plotTracks <- function(track){
  latlongtrack <- data.frame(track)
  latlongtrack$x <- 1000*latlongtrack$x + meanx
  latlongtrack$y <- 1000*latlongtrack$y + meany
  sp::coordinates(latlongtrack) <- c("x","y")
  proj4string(latlongtrack) <- myCRS
  latlongtrack <- data.frame(spTransform(latlongtrack,CRS=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")))
  
  par(mar = c(5, 5, 5, 2))
  
  plot(tracks@coords[,1], tracks@coords[,2], axes=FALSE, xlab="Longitude", ylab="Latitude", cex.lab = 1.8, col="white",xlim=c(min(tracks@coords[,1],latlongtrack$coords.x1),max(tracks@coords[,1],latlongtrack$coords.x1)),ylim=c(min(tracks@coords[,2],latlongtrack$coords.x2),max(tracks@coords[,2],latlongtrack$coords.x2)))
  degAxis(1, cex.axis = 1.5) # x-axis numbers
  degAxis(2, cex.axis = 1.5) # y-axis numbers
  raster::plot(bathdat, add=TRUE,
               image=TRUE,
               step = 10000,
               deepest=-5000,
               shallowest=2000,
               bpal=colorRampPalette(c("dodgerblue3","lightsteelblue2"))(100),
               col=grey(0,0.1),
               lwd=0.1)
  maps::map('worldHires', add=TRUE, fill=TRUE, col='white')
  lines(tracks@coords[,1], tracks@coords[,2],lwd=0.5,col=alpha("black",0.35))
  
  pal <- c("#78c679","#F0E442","#E31A1C","#88419d")
  
  for(id in unique(track$ID)) {
    ind <- which(track$ID==id)
    segments(x0 = latlongtrack$coords.x1[ind[-length(ind)]], y0 = latlongtrack$coords.x2[ind[-length(ind)]], x1 = latlongtrack$coords.x1[ind[-1]], y1 = latlongtrack$coords.x2[ind[-1]], 
             col = pal[latlongtrack$states[ind[-length(ind)]]], lwd=2)
  }
  
  legend("topleft", legend = stateNames,
         col = pal, lwd=2, bg="white", cex = 1.2)
}

plotCor <- function(data,model,sim,patchwork=TRUE){
  st <- viterbi(model)
  p <- list()
  p[["origData"]] <- cylcop::plot_joint_box(x=data$step*1000,theta=data$angle)
  p[["simData"]] <- cylcop::plot_joint_box(x=sim$step*1000,theta=sim$angle)
  p[["origData_state1"]] <- cylcop::plot_joint_box(x=data$step[which(st==1)]*1000,theta=data$angle[which(st==1)])
  p[["simData_state1"]] <- cylcop::plot_joint_box(x=sim$step[which(sim$states==1)]*1000,theta=sim$angle[which(sim$states==1)])
  p[["origData_state2"]] <- cylcop::plot_joint_box(x=data$step[which(st==2)]*1000,theta=data$angle[which(st==2)])
  p[["simData_state2"]] <- cylcop::plot_joint_box(x=sim$step[which(sim$states==2)]*1000,theta=sim$angle[which(sim$states==2)])
  p[["origData_state3"]] <- cylcop::plot_joint_box(x=data$step[which(st==3)]*1000,theta=data$angle[which(st==3)])
  p[["simData_state3"]] <- cylcop::plot_joint_box(x=sim$step[which(sim$states==3)]*1000,theta=sim$angle[which(sim$states==3)])
  p[["origData_state4"]] <- cylcop::plot_joint_box(x=data$step[which(st==4)]*1000,theta=data$angle[which(st==4)])
  p[["simData_state4"]] <- cylcop::plot_joint_box(x=sim$step[which(sim$states==4)]*1000,theta=sim$angle[which(sim$states==4)])
  if(patchwork) p$origData + p$simData + p$origData_state1 + p$simData_state1 + p$origData_state2 + p$simData_state2 + p$origData_state3 + p$simData_state3 + p$origData_state4 + p$simData_state4 + patchwork::plot_layout(nrow=5,byrow=TRUE)+plot_annotation(tag_levels="a")
  else return(p)
}

plotDens <- function(model,binwidth=5,stateNames=NULL,
                      title_sz = 18,  # Axis titles (e.g., "step length")
                      text_sz  = 14,  # Axis tick numbers
                      leg_sz   = 14  # Legend state names
                     ){
  
  plotDat <- model$data
  plotDat$step <- plotDat$step * 1000
  plotDat$states <- viterbi(model)
  if(is.null(stateNames)) stateNames <- model$stateNames
  stateNames <- c(stateNames,"total")
  
  stepdens <- list()
  
  if((model$conditions$dist[[1]] %in% momentuHMM:::rwdists) || isFALSE(model$conditions$DMind[[1]])){
    
    for(i in 1:(length(stateNames)-1)){
      stepdens[[i]] <- density(plotDat$step[which(plotDat$states==i)],na.rm=TRUE,adjust=2,from=0) 
      stepdens[[i]]$y <- stepdens[[i]]$y * mean(plotDat$states==i)
    }
    steptotal <- density(plotDat$step,na.rm=TRUE,adjust=2,from=0, bw = "SJ")
  
  } else {
    grid <- seq(min(model$data$step,na.rm=TRUE),max(model$data$step,na.rm=TRUE),length=1000) * 1000
    steptotal <- data.frame(x=grid,y=0)
    if(model$conditions$dist$step=="gamma"){
      for(i in 1:(length(stateNames)-1)){
        stepdens[[i]] <- data.frame(x=grid,y=NA)
        shape <- model$CIreal$step$est[1,i]^2/model$CIreal$step$est[2,i]^2
        scale <- 1000 * model$CIreal$step$est[2,i]^2/model$CIreal$step$est[1,i]
        stepdens[[i]]$y <- dgamma(grid,shape=shape,scale=scale) * mean(plotDat$states==i)
        steptotal$y <- steptotal$y + stepdens[[i]]$y
      }
    } else {
      genArgs <- list()
      genArgs[[1]] <- grid
      for(i in 1:(length(stateNames)-1)){
        stepdens[[i]] <- data.frame(x=grid,y=NA)
        genFun <- momentuHMM:::intdcrwrice
        genArgs[[2]] <- model$CIreal$step$est[1,i]
        genArgs[[3]] <- model$CIreal$step$est[2,i] * 1000
        genArgs[[4]] <- min(plotDat$step[which(plotDat$states==i)],na.rm=TRUE)
        genArgs[[5]] <- max(plotDat$step[which(plotDat$states==i)],na.rm=TRUE)
        d <- density(plotDat$step[which(plotDat$states==i)],na.rm=TRUE)
        genArgs[[6]] <- stats::approxfun(d$x,d$y)
        stepdens[[i]]$y <- do.call(genFun,genArgs) * mean(plotDat$states==i)
        steptotal$y <- steptotal$y + stepdens[[i]]$y
      }
    }
  }
  
  pal <- c("#78c679","#F0E442","#E31A1C","#88419d","black")
  
  stepdf <- rbind(data.frame(x=stepdens[[1]]$x,y=stepdens[[1]]$y,state=stateNames[1],col=pal[1]),
                  data.frame(x=stepdens[[2]]$x,y=stepdens[[2]]$y,state=stateNames[2],col=pal[2]),
                  data.frame(x=stepdens[[3]]$x,y=stepdens[[3]]$y,state=stateNames[3],col=pal[3]),
                  data.frame(x=stepdens[[4]]$x,y=stepdens[[4]]$y,state=stateNames[4],col=pal[4]),
                  data.frame(x=steptotal$x,y=steptotal$y,state=stateNames[5],col=pal[5]))
  stepdf$state <- factor(stepdf$state,levels=stateNames)
  COLORS <- pal
  names(COLORS) <- factor(stateNames,levels=stateNames)
  
  
  stepPlot <- ggplot(plotDat) + geom_histogram(aes(x=step,y=after_stat(density)),binwidth=binwidth,fill = "lightgrey", color = "black",na.rm=TRUE) +
    geom_line(data = stepdf[which(stepdf$state==stateNames[1]),], aes(x = x, y = y, color=stateNames[1]),linewidth=1.5)+
    geom_line(data = stepdf[which(stepdf$state==stateNames[2]),], aes(x = x, y = y, color=stateNames[2]),linewidth=1.5)+
    geom_line(data = stepdf[which(stepdf$state==stateNames[3]),], aes(x = x, y = y, color=stateNames[3]),linewidth=1.5)+
    geom_line(data = stepdf[which(stepdf$state==stateNames[4]),], aes(x = x, y = y, color=stateNames[4]),linewidth=1.5)+
    geom_line(data = stepdf[which(stepdf$state==stateNames[5]),], aes(x = x, y = y, color=stateNames[5]),linewidth=1.5,linetype=2)+
    scale_color_manual(values = COLORS, breaks=stateNames, name = "")+xlab("step length")+scale_x_continuous(breaks=seq(0,100,length=5),limits=c(0,100))
  
  angledens <- list()
  
  if((model$conditions$dist[[1]] %in% momentuHMM:::rwdists) || isFALSE(model$conditions$DMind[[1]])){
    
    for(i in 1:(length(stateNames)-1)){
      angledens[[i]] <- density(plotDat$angle[which(plotDat$states==i)],na.rm=TRUE,from=-pi,to=pi,adjust=2) 
      angledens[[i]]$y <- angledens[[i]]$y * mean(plotDat$states==i)
    }
    angletotal <- density(plotDat$angle,na.rm=TRUE,from=-pi,to=pi,adjust=2, bw = "SJ")
    
  } else {
    grid <- seq(-pi,pi,length=1000)
    angletotal <- data.frame(x=grid,y=0)
    if(model$conditions$dist$angle=="vm"){
      for(i in 1:(length(stateNames)-1)){
        angledens[[i]] <- data.frame(x=grid,y=NA)
        angledens[[i]]$y <- dvm(grid,0,model$CIreal$angle$est[1,i]) * mean(plotDat$states==i)
        angletotal$y <- angletotal$y + angledens[[i]]$y
      }
    } else {
      genArgs <- list()
      genArgs[[1]] <- grid
      for(i in 1:(length(stateNames)-1)){
        angledens[[i]] <- data.frame(x=grid,y=NA)
        genFun <- momentuHMM:::intdcrwvm
        genArgs[[2]] <- model$CIreal$angle$est[1,i]
        genArgs[[3]] <- model$CIreal$angle$est[2,i] * 1000
        genArgs[[4]] <- min(plotDat$step[which(plotDat$states==i)],na.rm=TRUE)
        genArgs[[5]] <- max(plotDat$step[which(plotDat$states==i)],na.rm=TRUE)
        d <- density(plotDat$step[which(plotDat$states==i)],na.rm=TRUE)
        genArgs[[6]] <- stats::approxfun(d$x,d$y)
        angledens[[i]]$y <- do.call(genFun,genArgs) * mean(plotDat$states==i)
        angletotal$y <- angletotal$y + angledens[[i]]$y
      }
    }
  }
  
  pal <- c("#78c679","#F0E442","#E31A1C","#88419d","black")
  
  angledf <- rbind(data.frame(x=angledens[[1]]$x,y=angledens[[1]]$y,state=stateNames[1],col=pal[1]),
                  data.frame(x=angledens[[2]]$x,y=angledens[[2]]$y,state=stateNames[2],col=pal[2]),
                  data.frame(x=angledens[[3]]$x,y=angledens[[3]]$y,state=stateNames[3],col=pal[3]),
                  data.frame(x=angledens[[4]]$x,y=angledens[[4]]$y,state=stateNames[4],col=pal[4]),
                  data.frame(x=angletotal$x,y=angletotal$y,state=stateNames[5],col=pal[5]))
  angledf$state <- factor(angledf$state,levels=stateNames)
  
  anglePlot <- ggplot(plotDat) + geom_histogram(aes(x=angle,y=after_stat(density)),breaks=seq(-pi,pi,length=30),fill = "lightgrey", color = "black",na.rm=TRUE) +
    geom_line(data = angledf[which(angledf$state==stateNames[1]),], aes(x = x, y = y, color=stateNames[1]),linewidth=1.5)+
    geom_line(data = angledf[which(angledf$state==stateNames[2]),], aes(x = x, y = y, color=stateNames[2]),linewidth=1.5)+
    geom_line(data = angledf[which(angledf$state==stateNames[3]),], aes(x = x, y = y, color=stateNames[3]),linewidth=1.5)+
    geom_line(data = angledf[which(angledf$state==stateNames[4]),], aes(x = x, y = y, color=stateNames[4]),linewidth=1.5)+
    geom_line(data = angledf[which(angledf$state==stateNames[5]),], aes(x = x, y = y, color=stateNames[5]),linewidth=1.5,linetype=2)+
    scale_color_manual(values = COLORS, breaks=stateNames, name = "")+scale_x_continuous(breaks=c(-pi, -pi/2, 0, pi/2, pi),labels=expression(-pi, -pi/2, 0, pi/2, pi))+xlab("turn angle (radians)")
  
  stepPlot <- stepPlot + theme(
    axis.title = element_text(size = title_sz),
    axis.text  = element_text(size = text_sz),
    legend.text = element_text(size = leg_sz)
  )
  
  anglePlot <- anglePlot + theme(
    axis.title = element_text(size = title_sz),
    axis.text  = element_text(size = text_sz),
    legend.text = element_text(size = leg_sz)
  )
  
  return(list(stepPlot=stepPlot,anglePlot=anglePlot))
}

prPlot <- function(m,what="step"){
  pr <- momentuHMM::pseudoRes(m)
  prdf <- data.frame(do.call(cbind,pr))
  colnames(prdf) <- what
  pl <- vector('list',length(what))
  names(pl) <- what
  for(i in what){
    
    pl[[i]] <- list()
    ylim <- range(pr[[paste0(i,"Res")]][which(is.finite(pr[[paste0(i,"Res")]]))],na.rm=TRUE)
    pl[[i]]$pr <- ggplot(data=prdf) + geom_line(aes(x=1:nrow(prdf),y=.data[[i]])) + ylab(paste0(i," pseudo-residuals")) + xlab("Observation index") + theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + scale_y_continuous(breaks=seq(-7,7,by=3.5),limits=c(-7,7))
    
    n <- length(prdf[[i]][which(!is.na(prdf[[i]]))])
    qq_df <- data.frame(sample_quantiles = sort(prdf[[i]]))
    qq_df$theoretical_quantiles <- qnorm(ppoints(n))
    q_sample <- quantile(qq_df$sample_quantiles, probs = c(0.25, 0.75), na.rm = TRUE)
    q_theoretical <- qnorm(c(0.25, 0.75))
    
    b <- (q_sample[2] - q_sample[1]) / (q_theoretical[2] - q_theoretical[1])
    a <- q_sample[1] - b * q_theoretical[1]
    
    qq_df$fit <- a + b * qq_df$theoretical_quantiles
    
    p <- ppoints(n)
    conf <- 0.95 # Standard 95% confidence level
    zz <- qnorm(1 - (1 - conf) / 2)
    SE <- (b / dnorm(qq_df$theoretical_quantiles)) * sqrt(p * (1 - p) / n)
    
    qq_df$upper_bound <- qq_df$fit + zz * SE
    qq_df$lower_bound <- qq_df$fit - zz * SE
    
    pl[[i]]$qq <- ggplot(qq_df, aes(x = theoretical_quantiles, y = sample_quantiles)) +
      
      geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "gray", alpha = 0.4) +
      
      geom_line(aes(y = fit), color = "black", linewidth = 1) +
      
      geom_point(color = "red") +
      
      labs(
        x = "Theoretical Quantiles",
        y = "Sample Quantiles",
        title = ""
      ) +
      
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + scale_y_continuous(breaks=seq(-7,7,by=3.5),limits=c(-7,7))
    
    acf_data <- stats::acf(prdf[[i]], na.action = na.pass, plot = FALSE)
    
    acf_df <- data.frame(
      lag = acf_data$lag,
      acf = acf_data$acf
    )
    
    conf_level <- 0.95
    ci <- qnorm((1 + conf_level) / 2) / sqrt(acf_data$n.used)
    
    pl[[i]]$acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
      
      geom_segment(data = acf_df[acf_df$lag > 0,], aes(xend = lag, yend = 0), linewidth = 1) +
      geom_hline(yintercept = c(ci, -ci), color = "steelblue", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "black") +
      
      labs(
        x = "Lag",
        y = "ACF"#,
        #title = paste("Series:", deparse(substitute(pr)))
      ) +
      
      theme_bw() +
      theme(
        # Remove the default grid lines for a cleaner look like the base plot.
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + coord_cartesian(ylim = c(-0.1, 1))
    
  }
  return(pl)
}

plotHist <- function(model,what="bath",state,binwidth=0.5,xl=c(-6,0),yl=c(0,550),xla="sea floor depth (km)",title=waiver()){
  st <- viterbi(model)
  df <- data.frame(cov=model$data[[what]][which(st==state)])
  ggplot(df) + geom_histogram(aes(x=cov),binwidth=binwidth,boundary=xl[1],fill = "lightgrey", color = "black",na.rm=TRUE) + xlim(xl) + ylim(yl) + xlab(xla) + ylab("Frequency") + ggtitle(title) + theme(text = element_text(size = 15))
}

create_latex_table <- function(dfmu, dfbeta) {
  
  format_param_name <- function(name) {
    parts <- str_split(name, ":", n = 2, simplify = TRUE)
    main_param <- str_trim(parts[1])
    covariate <- str_trim(parts[2])
    
    if (str_detect(main_param, "^mean") || str_detect(main_param, "^sd")) {
      
      if (str_detect(main_param, "^mean")) {
        subscripts <- str_match(main_param, "^mean\\.([xy])_(\\d+)$")
        p <- paste0("$\\mu_{", subscripts[2], ",", subscripts[3], "}$")
        greek_letter <- "delta"
      } else { # sd case
        subscripts <- str_match(main_param, "^sd\\.[xy]_(\\d+)$")
        p <- paste0("$\\sigma_{", subscripts[2], "}$")
        greek_letter <- "beta"
      }
      
      if (str_detect(covariate, "crw")) {
        c <- ": $\\beta_{crw}$"
      } else if (covariate == "(Intercept)") {
        c <- paste0(": $\\", greek_letter, "_0$")
      } else {
        covariate_clean <- ifelse(str_detect(covariate, "boundary"), "bound", str_remove(covariate, "\\.[xy]$"))
        c <- paste0(": $\\", greek_letter, "_{", covariate_clean, "}$")
      }
      return(paste0(p, c))
      
    } else {
      trans_parts <- str_extract_all(main_param, "\\d+")[[1]]
      p <- paste0("$\\gamma_{", trans_parts[1], trans_parts[2], "}$")
      
      if (covariate == "(Intercept)") {
        c <- ": $\\beta_0$"
      } else if (str_detect(covariate, "I\\(")) {
        inner_cov <- str_match(covariate, "I\\((.*)\\)")[, 2]
        c <- paste0(": $\\beta_{", inner_cov, "}$")
      } else {
        c <- paste0(": $\\beta_{", covariate, "}$")
      }
      return(paste0(p, c))
    }
  }
  
  # Define the caption and header string
  header_string <- "\\caption{Southern elephant seal parameter estimates from the 4-state position-based potential function model. Parameters are on the real scale for the state-specific mean paramters ($\\mu_{c,s}; c \\in \\{x,y\\}, s \\in \\{1,2,3,4\\}$), the log scale for the state-specific speed parameters ($\\sigma_s$), and the logit scale for the state transition probability parameters ($\\gamma_{ij}; i,j \\in \\{1,2,3,4\\}$). For $\\mu_{c,s}$, $\\beta_{crw}$ is the lag 1 autocorrelation parameter and $\\delta_k$ is the potential function coefficient for covariate $k$. For $\\sigma_s$ and $\\gamma_{ij}$, $\\beta_0$ is the intercept term and $\\beta_k$ is the coefficent for covariate $k$. Covariates included Euclidean distance to coast (``d2coast''), anisotropically-scaled Euclidean distance from colony (``dfrcol''), study area boundary (``bound''), Euclidean distance to colony (``d2col''), sea floor depth (``bath''), Euclidean distance inland (``land''), and time since departure from the colony (``time'').} \\label{tab:model_params} \\\\\n\n\\toprule\nParameter & Est & SE & LCI & UCI \\\\\n\\bottomrule\n"
  
  latex_string <- paste0("\\begin{longtable}{@{}lrrrr@{}}\n", header_string)
  
  # Process dfmu
  for (i in 1:nrow(dfmu)) {
    param_latex <- format_param_name(rownames(dfmu)[i])
    vals <- sprintf("%.2f", as.numeric(dfmu[i,]))
    row_string <- paste(param_latex, paste(vals, collapse = " & "), sep = " & ")
    latex_string <- paste0(latex_string, row_string, " \\\\\n")
    if (i < nrow(dfmu) && str_extract(rownames(dfmu)[i], "^[^:]+") != str_extract(rownames(dfmu)[i+1], "^[^:]+")) {
      latex_string <- paste0(latex_string, "\\hdashline\n")
    }
  }
  
  latex_string <- paste0(latex_string, "\\hdashline\n")
  
  # Process dfbeta
  for (i in 1:nrow(dfbeta)) {
    param_latex <- format_param_name(rownames(dfbeta)[i])
    vals <- sprintf("%.2f", as.numeric(dfbeta[i,]))
    row_string <- paste(param_latex, paste(vals, collapse = " & "), sep = " & ")
    latex_string <- paste0(latex_string, row_string, " \\\\\n")
    if (i < nrow(dfbeta) && str_extract(rownames(dfbeta)[i], "^[^:]+") != str_extract(rownames(dfbeta)[i+1], "^[^:]+")) {
      latex_string <- paste0(latex_string, "\\hdashline\n")
    }
  }
  
  latex_string <- paste0(latex_string, "\\bottomrule\n\\end{longtable}\n")
  
  cat(latex_string)
}

cm_to_latex <- function(cm,
                        caption = "My Caption",
                        label = "tab:my-label",
                        pred_name = "Prediction",
                        ref_name = "Reference",
                        row_names = NULL,
                        col_names = NULL) {
  
  # Extract the confusion matrix table
  mat <- as.matrix(cm$table)
  
  # Use provided names or default to original dimnames
  if (is.null(row_names)) {
    row_names <- rownames(mat)
  }
  if (is.null(col_names)) {
    col_names <- colnames(mat)
  }
  
  # Calculate row and column totals
  row_totals <- rowSums(mat)
  col_totals <- colSums(mat)
  
  # Combine matrix with totals
  mat_ext <- rbind(cbind(mat, Total = row_totals), Total = c(col_totals, sum(col_totals)))
  
  # Update row and column names for the extended matrix
  final_row_names <- c(row_names, "Total")
  final_col_names <- c(col_names, "Total")
  
  # --- Start building the LaTeX string ---
  
  # Header
  latex_string <- paste0(
    "\\begin{table}[h!]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{@{}l", paste(rep("c", ncol(mat_ext)), collapse = ""), "@{}}\n",
    "\\toprule\n"
  )
  
  # Column titles
  num_cols <- length(final_col_names)
  latex_string <- paste0(
    latex_string,
    "& \\multicolumn{", num_cols, "}{c}{\\textbf{", ref_name, "}} \\\\\n",
    "\\cmidrule(l){2-", num_cols + 1, "}\n"
  )
  
  # Sub-column titles
  latex_string <- paste0(
    latex_string,
    "\\textbf{", pred_name, "} & ",
    paste(final_col_names, collapse = " & "),
    " \\\\\n",
    "\\midrule\n"
  )
  
  # Table body
  for (i in 1:nrow(mat_ext)) {
    # The last row's final cell should be empty
    if (i == nrow(mat_ext)) {
      mat_ext[i, ncol(mat_ext)] <- ""
    }
    
    row_content <- paste(mat_ext[i, ], collapse = " & ")
    latex_string <- paste0(
      latex_string,
      "\\text{", final_row_names[i], "} & ",
      row_content,
      " \\\\\n"
    )
  }
  
  # Footer
  latex_string <- paste0(
    latex_string,
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  # Print the result to the console
  cat(latex_string)
  
  # Invisibly return the string
  invisible(latex_string)
}
