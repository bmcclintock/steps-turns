# Add this line to the top of your script
library(ggrastr)

prPlot <- function(pr, what = "step") {
  prdf <- data.frame(do.call(cbind, pr))
  colnames(prdf) <- what
  pl <- vector('list', length(what))
  names(pl) <- what
  for (i in what) {
    
    pl[[i]] <- list()
    ylim <- range(pr[[paste0(i, "Res")]][which(is.finite(pr[[paste0(i, "Res")]]))], na.rm = TRUE)
    # This plot is simple and likely doesn't need rasterization
    pl[[i]]$pr <- ggplot(data = prdf) + geom_line(aes(x = 1:nrow(prdf), y = .data[[i]])) + ylab(paste0(i, " pseudo-residuals")) + xlab("Observation index") + theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + scale_y_continuous(breaks = seq(-7, 7, by = 3.5), limits = c(-7, 7))
    
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
      
      # MODIFICATION: Rasterize the complex ribbon and point layers
      rasterise(geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "gray", alpha = 0.4), dpi = 300) +
      rasterise(geom_point(color = "red"), dpi = 300) +
      
      # Keep the fit line as a sharp vector
      geom_line(aes(y = fit), color = "black", linewidth = 1) +
      
      labs(
        x = "theoretical quantiles",
        y = paste(i, "sample quantiles"),
        title = ""
      ) +
      
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + scale_y_continuous(breaks = seq(-7, 7, by = 3.5), limits = c(-7, 7))
    
    acf_data <- stats::acf(prdf[[i]], na.action = na.pass, plot = FALSE)
    
    acf_df <- data.frame(
      lag = acf_data$lag,
      acf = acf_data$acf
    )
    
    conf_level <- 0.95
    ci <- qnorm((1 + conf_level) / 2) / sqrt(acf_data$n.used)
    
    pl[[i]]$acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
      
      # MODIFICATION: Rasterize the segments which can be numerous
      rasterise(geom_segment(data = acf_df[acf_df$lag > 0, ], aes(xend = lag, yend = 0), linewidth = 1), dpi = 300) +
      
      geom_hline(yintercept = c(ci, -ci), color = "steelblue", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "black") +
      
      labs(
        x = "Lag",
        y = "ACF"
      ) +
      
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + coord_cartesian(ylim = c(-0.3, 0.75))
    
  }
  return(pl)
}

compressPDF <- function(input_pdf, output_pdf, quality = "/ebook") {
  
  if (missing(output_pdf)) output_pdf <- input_pdf
  
  if (!file.exists(input_pdf)) {
    stop("Error: Input file not found at the specified path!")
  }
  
  # This command includes the -dPDFSETTINGS flag for compression
  gs_command <- sprintf(
    "gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=%s -dNOPAUSE -dQUIET -dBATCH -sOutputFile=%s %s",
    quality,
    shQuote(output_pdf),
    shQuote(input_pdf)
  )
  
  print(gs_command)
  system(gs_command)
}

prcrwplot <- prcrwSTplot <- prSTplot <- list()
for(sc in 1:length(scenarios)){
  prcrwplot[[scenarios[sc]]] <- list()
  prcrwplot[[scenarios[sc]]]$crwpos <- list()
  prcrwplot[[scenarios[sc]]]$crwST <- list()
  prcrwplot[[scenarios[sc]]]$ST <- list()
  
  prcrwSTplot[[scenarios[sc]]] <- list()
  prcrwSTplot[[scenarios[sc]]]$crwpos <- list()
  prcrwSTplot[[scenarios[sc]]]$crwST <- list()
  prcrwSTplot[[scenarios[sc]]]$ST <- list()
  
  prSTplot[[scenarios[sc]]] <- list()
  prSTplot[[scenarios[sc]]]$crwpos <- list()
  prSTplot[[scenarios[sc]]]$crwST <- list()
  prSTplot[[scenarios[sc]]]$ST <- list()
  
  load(paste0("simulations/crwSims_",scenarios[sc],".RData"))  
  for(i in 1:length(crwpos)){
    prcrwplot[[scenarios[sc]]]$crwpos[[i]] <- prPlot(crwpos[[i]]$pseudoRes,what="mu")
    prcrwplot[[scenarios[sc]]]$crwST[[i]] <- prPlot(crwST[[i]]$pseudoRes,what=c("step","angle"))
    prcrwplot[[scenarios[sc]]]$ST[[i]] <- prPlot(ST[[i]]$pseudoRes,what=c("step","angle"))
  }
  
  load(paste0("simulations/crwSTSims_",scenarios[sc],".RData"))  
  for(i in 1:length(crwpos)){
    prcrwSTplot[[scenarios[sc]]]$crwpos[[i]]  <- prPlot(crwpos[[i]]$pseudoRes,what="mu")
    prcrwSTplot[[scenarios[sc]]]$crwST[[i]]  <- prPlot(crwST[[i]]$pseudoRes,what=c("step","angle"))
    prcrwSTplot[[scenarios[sc]]]$ST[[i]]  <- prPlot(ST[[i]]$pseudoRes,what=c("step","angle"))
  }
  
  load(paste0("simulations/STSims_",scenarios[sc],".RData"))  
  for(i in 1:length(crwpos)){
    prSTplot[[scenarios[sc]]]$crwpos[[i]] <- prPlot(crwpos[[i]]$pseudoRes,what="mu")
    prSTplot[[scenarios[sc]]]$crwST[[i]] <- prPlot(crwST[[i]]$pseudoRes,what=c("step","angle"))
    prSTplot[[scenarios[sc]]]$ST[[i]] <- prPlot(ST[[i]]$pseudoRes,what=c("step","angle"))
  }
}
  
ind <- 70 #sample.int(100,1) # 20, 70, 89 #
for(pl in c("pr","qq","acf")){
  for(what in c("step","angle")){
    
    pdf(paste0("simulations/prcrw_",what,"_",pl,".pdf"),width=11,height=11)
    print(
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "crwpos", angle = 0, size = 5) + 
        theme_void() + prcrwplot$low$crwpos[[ind]]$mu[[pl]] + ggtitle("low") + prcrwplot$moderate$crwpos[[ind]]$mu[[pl]] + ggtitle("moderate") + prcrwplot$high$crwpos[[ind]]$mu[[pl]] + ggtitle("high") +
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "crwST", angle = 0, size = 5) + 
        theme_void() + prcrwplot$low$crwST[[ind]][[what]][[pl]] + prcrwplot$moderate$crwST[[ind]][[what]][[pl]] + prcrwplot$high$crwST[[ind]][[what]][[pl]] +
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "ST", angle = 0, size = 5) + 
        theme_void() + prcrwplot$low$ST[[ind]][[what]][[pl]] + prcrwplot$moderate$ST[[ind]][[what]][[pl]] + prcrwplot$high$ST[[ind]][[what]][[pl]] + plot_layout(nrow=3,ncol=4,byrow=TRUE,guides="collect",axes="collect", widths = c(1, 5, 5, 5))
    )
    dev.off()
    
    pdf(paste0("simulations/prcrwST_",what,"_",pl,".pdf"),width=11,height=11)
    print(
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "crwpos", angle = 0, size = 5) + 
        theme_void() + prcrwSTplot$low$crwpos[[ind]]$mu[[pl]] + ggtitle("low") + prcrwSTplot$moderate$crwpos[[ind]]$mu[[pl]] + ggtitle("moderate") + prcrwSTplot$high$crwpos[[ind]]$mu[[pl]] + ggtitle("high") +
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "crwST", angle = 0, size = 5) + 
        theme_void() + prcrwSTplot$low$crwST[[ind]][[what]][[pl]] + prcrwSTplot$moderate$crwST[[ind]][[what]][[pl]] + prcrwSTplot$high$crwST[[ind]][[what]][[pl]] +
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "ST", angle = 0, size = 5) + 
        theme_void() + prcrwSTplot$low$ST[[ind]][[what]][[pl]] + prcrwSTplot$moderate$ST[[ind]][[what]][[pl]] + prcrwSTplot$high$ST[[ind]][[what]][[pl]] + plot_layout(nrow=3,ncol=4,byrow=TRUE,guides="collect",axes="collect", widths = c(1, 5, 5, 5))
    )
    dev.off()
    
    pdf(paste0("simulations/prST_",what,"_",pl,".pdf"),width=11,height=11)
    print(
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "crwpos", angle = 0, size = 5) + 
        theme_void() + prSTplot$low$crwpos[[ind]]$mu[[pl]] + ggtitle("low") + prSTplot$moderate$crwpos[[ind]]$mu[[pl]] + ggtitle("moderate") + prSTplot$high$crwpos[[ind]]$mu[[pl]] + ggtitle("high") +
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "crwST", angle = 0, size = 5) + 
        theme_void() + prSTplot$low$crwST[[ind]][[what]][[pl]] + prSTplot$moderate$crwST[[ind]][[what]][[pl]] + prSTplot$high$crwST[[ind]][[what]][[pl]] +
      ggplot() + 
        annotate(geom = "text", x = 0, y = 0, label = "ST", angle = 0, size = 5) + 
        theme_void() + prSTplot$low$ST[[ind]][[what]][[pl]] + prSTplot$moderate$ST[[ind]][[what]][[pl]] + prSTplot$high$ST[[ind]][[what]][[pl]] + plot_layout(nrow=3,ncol=4,byrow=TRUE,guides="collect",axes="collect", widths = c(1, 5, 5, 5))
    )
    dev.off()
  }
}

for(m in c("crwpos","crwST","ST")){
  if(m=="crwpos") what <- "mu"
  else what <- c("step","angle")
  for(j in what){
    for(pl in c("pr","qq","acf")){
      for(sc in scenarios){
        

        pdf(paste0("simulations/prcrw_",sc,"_",m,"_",j,"_",pl,".pdf"),width=11,height=11)
        print(wrap_plots(lapply(prcrwplot[[sc]][[m]],function(x) x[[j]][[pl]]),10,10,guides="collect",axes="collect"))
        dev.off()
        
        pdf(paste0("simulations/prcrwST_",sc,"_",m,"_",j,"_",pl,".pdf"),width=11,height=11)
        print(wrap_plots(lapply(prcrwSTplot[[sc]][[m]],function(x) x[[j]][[pl]]),10,10,guides="collect",axes="collect"))
        dev.off()
        
        pdf(paste0("simulations/prST_",sc,"_",m,"_",j,"_",pl,".pdf"),width=11,height=11)
        print(wrap_plots(lapply(prSTplot[[sc]][[m]],function(x) x[[j]][[pl]]),10,10,guides="collect",axes="collect"))
        dev.off()
        
        
      }
    }
  }
}

## compile appendix from cammand line to generate .aux file for references
# $ pdflatex Appendix_S1.tex
## compress PDF
# tools::compactPDF("manuscript/Appendix_S1.pdf", gs_quality = "ebook")
