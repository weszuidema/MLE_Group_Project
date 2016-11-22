## Binned prediction plots and ROC plots for binary GLMs
## Christopher Adolph    faculty.washington.edu/cadolph
## 23 October 2016
##
## plot.binPredict is general but requires the tile package;
## the method binPredict could be extended to non-GLM binary classifiers
##
## TO DO: extend to count models, linear regression
##
## Typical usage:
##
## plot(binPredict(glm.m1, col="blue", label="M1"),
##      binPredict(glm.m2, col="red", label="M2"))
##
######################################################################

binPredict <- function(x,...)
  UseMethod("binPredict")

binPredict.glm <- function(x,  # glm() model object
                           newdata=NULL, # a data frame with alternative data for out of sample tests
                           newpred=NULL, # a vector with alternative predictions, e.g. from LOO-CV
                           col="black", # color of this model in plot
                           label=NULL, # a name for this model (to label plot)
                           bins=20, # a vector of breakpoints or number of quantiles
                           quantiles=FALSE, # force all bins to same number of obs per bin
                           sims=1000  # if 0, use point estimates for bin simulation (ROC always uses PE's)                    
) {
  require(verification)  
  
  # Preliminaries
  if (!any(class(x)=="glm")) stop("x must be a glm object.")
  #if ((overlap<0)||(overlap>=1)) stop("overlap must be >=0 and <1.")
  if (length(bins)==1) {
    bins <- floor(bins)
    if (bins<1)
      stop("bins must be a scalar >1 or a vector with values between 0 and 1.")
  } else {
    if (any((bins<0)|(bins>1)))
      stop("bins must be a scalar >1 or a vector with values between 0 and 1.")
  }    
  model <- x$formula
  if (is.null(newdata)) {
    data <- na.omit(x$data)
    outcome <- x$y
  }
  else {
    data <- na.omit(newdata)
    outcome <- data[[all.vars(model)[1]]]
  }
  n <- nrow(data)
  pe <- coef(x)
  fm <- x$family$family
  fn <- x$family$linkinv
  if (sims) {
    require(MASS)
    vc <- vcov(x)
    simpe <- mvrnorm(sims,pe,vc)
  } else {
    simpe <- t(as.matrix(pe))
    sims <- 1
  }
  
  # Simulate from model posterior using simcf-style functions and code
  xhyp <- model.matrix(model, data)
  if (is.null(newpred)) {
    ypred <- as.vector(sapply(split(xhyp, row(xhyp)),
                              function(x,simpe,fn)fn(simpe%*%x), simpe, fn))
  } else {
    ypred <- newpred
    sims <- 1
  }
  yact <- as.vector(sapply(outcome, function(x,sims)rep(x,sims), sims)) 
  
  # Simple predictions from point estimates
  if (is.null(newpred)) {
    ypred0 <- as.vector(sapply(split(xhyp, row(xhyp)),
                               function(x,simpe,fn)fn(simpe%*%x), t(as.matrix(pe)), fn))
  } else {
    ypred0 <- newpred
  }
  
  # Create bounds on each bin
  # to FINISH
  lbins <- length(bins)
  if ((lbins==1)&&(bins>=1)) {
    limits <- seq(0,1,length.out=bins+1)
    if (quantiles) limits <- quantile(ypred, probs=limits)
    bounds <- cbind(limits[1:bins], limits[2:(bins+1)])
    #print(bounds)
  } else {
    #print(bins)
    if (is.unsorted(bins,strictly=TRUE)) {
      stop("bins must be a positive integer greater than one or an increasing vector in [0,1].")
    } else {
      limits <- bins
      if (!identical(bins[lbins],1)) limits <- c(limits,1)
      if (!identical(bins[1],0)) limits <- c(0,limits)
      nbins <- length(limits)-1
      bounds <- cbind(limits[1:nbins], limits[2:(nbins+1)])
      #print(bounds)
    }
  }
  
  # Bin simulations according to bounds  
  res <- sapply(split(bounds,row(bounds)),
                function(x,ypred,yact,sims){cond <- (ypred<=x[2])&(ypred>=x[1]);
                c(mean(ypred[cond]),
                  mean(yact[cond]),
                  sum(cond)/sims)
                },
                ypred, yact, sims)
  res[is.nan(res)] <- NA
  
  # Compute prediction error:  actual/prediction
  if (any(!is.na(match(fm, c("gaussian","poisson","quasipoisson"))))) {
    yerr <- res[2,]/res[1,]
    yerr[yerr==Inf] <- NA
    yerr[yerr==-Inf] <- NA
  }
  if (any(!is.na(match(fm, c("binomial","quasibinomial"))))) {
    yerr <- res[2,]/res[1,]
    yerr[yerr==Inf] <- NA
    yerr[yerr==-Inf] <- NA
  }
  
  # receiver operator curve  (avoid simulation by default -- TO ADD)
  ypredS <- rev(ypred0[order(ypred0)])
  yactS <- rev(outcome[order(ypred0)])
  ctr <- xROC <- yROC  <- 0
  tROC <- 1
  for (i in 2:length(ypredS)) {
    ctr <- ctr+1
    if (!identical(yactS[i],yactS[i-1])) {
      if (yactS[i-1]) {
        xROC <- c(xROC, xROC[length(xROC)])
        yROC <- c(yROC, yROC[length(yROC)] + ctr)
      } else {
        xROC <- c(xROC, xROC[length(xROC)] + ctr)
        yROC <- c(yROC, yROC[length(yROC)])
      }
      tROC <- c(tROC, ypredS[i])
      ctr <- 0
    }
  }
  xROC <- c(xROC/sum(1-yactS), 1)
  yROC <- c(yROC/sum(yactS), 1)
  tROC <- c(tROC, 0)
  auc <- roc.area(yactS, ypredS)$A
  
  # Return as binPredict object
  res <- list(ypred=res[1,], yact=res[2,], yerr=yerr, roc=cbind(xROC, yROC, tROC), auc=auc, n=res[3,], bounds=bounds, col=col, label=label)
  class(res) <- "binPredict"
  res
}


giantListVector <- function(res,...) {
  x <- list(...)
  e <- 0
  for (i in 1:length(x)) {
    if (length(x[[i]])) {
      s <- e+1
      e <- (s-1)+length(x[[i]])
      for (j in 1:length(x[[i]])) 
        res[[(s:e)[j]]] <- x[[i]][[j]]
    }
  }  
  res
}

auc <- function(x) x$auc


plot.binPredict <- function(x, #binPredict object
                            ..., # additional binPredict objects (optional)
                            together=TRUE,  # logical, plot models on same plot?
                            display=c("avp", "evp", "roc"), # character or character vector,
                            # avp:  plot predicted actual vs predicted probs
                            # evr:  plot actual/predicted vs predicted probs
                            # roc:  plot receiver operator characteristic curves 
                            thresholds=NULL, # numeric vector, show these thresholds on ROC plot
                            hide=NULL, # logical, hide number of observations in each bubble?
                            ignore=5, # numeric, ignore bins with fewer obs than this
                            totalarea=0.1, # numeric, a model's total circle area relative to plot
                            labx=0.1, # numeric, absolute x location of model labels
                            cex=NULL, # numeric, override automatic character size
                            border=FALSE, # logical, draw outlines around each bubble
                            alpha=0.5, # numeric, transparency of bubbles (1 for opaque)
                            clip=FALSE, # logical, clip bubbles to plotting area?
                            showbins=FALSE, # logical, show boundaries of bins
                            bottom=FALSE, # logical, when showing #obs/bin, print this at the bottom
                            file=NULL  # character, save result to this pdf file instead of screen
) {
  
  # Preliminaries
  require(tile)
  #display <- display[1]
  traceData <- list(x,...)
  nTraces <- length(traceData)
  nPlots <- ((!together)*nTraces + together)*length(display)
  nicegray <- "gray45"
  
  # Size of characters
  if (is.null(cex)) {
    if (is.null(file)) {
      cex <- 0.9
    } else {
      cex <- 0.75
    }
  }
  
  # Allow clipping?
  if (clip) clip <- "on" else clip <- "off"
  
  # Create traces from each binPredict object
  bubbleTraces <- numberTraces <- labelTraces <-
    bubble2Traces <- number2Traces <- label2Traces <-
    bubble3Traces <- label3Traces <- auc3Traces <- border3Traces <- thresh3Trace <- threshL3Trace <- 
    vector("list", nTraces)
  
  if (together) {
    binTraces <- bin2Traces <- lineTraces <- line2Traces <- line3Traces <- auctitle3Traces <- vector("list", 1)
  } else {
    binTraces <- bin2Traces <- lineTraces <- line2Traces <- line3Traces <- auctitle3Traces <- vector("list", nTraces)
  }
  
  # Quick loop to get plotting area for error plots, number of obs, etc.
  rmax <- 0
  rmin <- 1
  nbub <- NULL
  for (i in 1:nTraces) {
    
    # Grab current model
    ct <- traceData[[i]]
    
    # Ignore bins with too few cases
    ct$yerr[ct$n<=ignore] <- ct$ypred[ct$n<=ignore] <- ct$yact[ct$n<=ignore] <- NA
    
    # Warn if ignoring small bins eliminates all plotting for this model
    if (all(ct$n<=ignore))
      warning(paste0("Model \"", ct$label, "\" has no bins larger than the \"ignore\" threshold.
                     Perhaps you should lower \"ignore\", widen the bins, or increase the number of bins."))
    
    # Get running max error
    rmax <- max(c(rmax, ct$yerr), na.rm=TRUE)
    rmin <- min(c(rmin, ct$yerr), na.rm=TRUE)
    
    # Collect n
    nbub <- c(nbub, ct$n)
  }
  if (is.null(hide)) {
    if (length(unique(floor(nbub)))==1) hide <- TRUE else hide <- FALSE
  }
  
  # In case all elements eliminated by ignore, set a default rmax
  if (length(rmax)==0||rmax==0) rmax <- 2
  if (length(rmin)==0||rmin==0) rmin <- 0.5
  
  if (rmin>0) {
    logAxis <- TRUE
    rmax <- max(abs(log(c(rmin,rmax))))
    rmin <- exp(-rmax)
    rmax <- exp(rmax)
  } else {
    logAxis <- FALSE
  }
  
  # Loop again to do plotting
  for (i in 1:nTraces) {
    
    # Grab current model
    ct <- traceData[[i]]
    
    # Ignore bins with too few cases
    ct$yerr[ct$n<=ignore] <- ct$ypred[ct$n<=ignore] <- ct$yact[ct$n<=ignore] <- NA
    
    # Adjust radii of circles
    radius <- ct$n*sqrt(totalarea/(pi*sum(ct$n^2)))
    
    # Option to show border of circles
    colOpt <- NA
    if (border) colOpt <- ct$col
    
    # Add actual vs predicted
    if (any(display=="avp")) {
      
      # Find current column and plot number  
      curcol <- which(display=="avp")
      curplot <- (!together)*(i-1)*length(display) + curcol
      
      # Bubble plot
      bubbleTraces[[i]] <- circleTile(x=ct$ypred,
                                      y=ct$yact,
                                      r=radius,
                                      col=colOpt,
                                      fill=ct$col,
                                      alpha=alpha,
                                      layer=10,
                                      clip=clip,
                                      plot=curplot
      )
      
      if (hide) nlab <- rep("",length(ct$ypred)) else nlab <- round(ct$n,0)
      if (bottom) ynum <- rep(c(0.03,0.07), length(ct$ypred))[1:length(ct$ypred)] else ynum <- ct$yact
      numberTraces[[i]] <- textTile(x=ct$ypred,
                                    y=ynum,
                                    labels=nlab,
                                    col=ct$col,
                                    cex=cex-0.3,
                                    layer=1,
                                    clip=clip,
                                    plot=curplot
      )
      
      # Ref Line
      if (together) {
        lineTraces[[1]] <- linesTile(x=c(0,1), y=c(0,1),
                                     layer=15,
                                     plot=curplot                                   
        )
      } else {
        lineTraces[[i]] <- linesTile(x=c(0,1), y=c(0,1),
                                     layer=15,
                                     plot=curplot                                   
        )
      }
      
      # Model Label
      if (together) {yloc <- 1-(i/12)} else {yloc <- 0.95}
      labelTraces[[i]] <- textTile(labels=ct$label,
                                   x=labx,
                                   y=yloc,
                                   col=ct$col,
                                   cex=cex+0.1,
                                   layer=0,
                                   clip=clip,
                                   plot=curplot
      )
      # Bin boundaries
      if (together) {
        if (i==1) {
          if (showbins) {
            cb <- unique(vec(ct$bounds))
            binTraces[[1]] <- polylinesTile(x=as.numeric(vec(t(cbind(cb,cb)))),
                                            y=rep(c(-0.05,1.05), length(cb)),
                                            id=as.numeric(vec(t(cbind(1:length(cb), 1:length(cb))))),
                                            lty="dotted",
                                            col=nicegray,
                                            layer=20,
                                            plot=curplot
            )
          }
        }
      } else {
        if (showbins) {
          cb <- unique(vec(ct$bounds))
          binTraces[[i]] <- polylinesTile(x=as.numeric(vec(t(cbind(cb,cb)))),
                                          y=rep(c(-0.05,1.05), length(cb)),
                                          id=as.numeric(vec(t(cbind(1:length(cb), 1:length(cb))))),
                                          lty="dotted",
                                          col=nicegray,
                                          layer=20,
                                          plot=curplot
          )
        }
      }
    }
    
    
    
    # Add error rate vs predicted
    if (any(display=="evp")) {
      
      # Find current column and plot number  
      curcol <- which(display=="evp")
      curplot <- (!together)*(i-1)*length(display) + curcol
      
      # Bubble plot 2
      bubble2Traces[[i]] <- circleTile(x=ct$ypred,
                                       y=ct$yerr,
                                       r=unit(radius,"npc"),
                                       col=colOpt,
                                       fill=ct$col,
                                       alpha=alpha,
                                       clip=clip,
                                       layer=10,
                                       plot=curplot
      )
      
      if (hide) nlab <- rep("",length(ct$ypred)) else nlab <- round(ct$n,0)
      number2Traces[[i]] <- textTile(x=ct$ypred,
                                     y=ct$yerr,
                                     labels=nlab,
                                     col=ct$col,
                                     cex=cex,
                                     layer=1,
                                     clip=clip,
                                     plot=curplot
      )
      
      
      # Ref Line 2
      if (together) {
        line2Traces[[1]] <- linesTile(x=c(0,1), y=c(1,1),
                                      layer=15,
                                      plot=curplot                                   
        )
      } else {
        line2Traces[[i]] <- linesTile(x=c(0,1), y=c(1,1),
                                      layer=15,
                                      plot=curplot                                   
        )
      }
      
      # Model Label 2
      if (together) {yloc <- 1-(i/12)} else {yloc <- 0.95}
      if (logAxis)
        ys <- exp(log(rmax) - (log(rmax)-log(rmin))*(1-yloc))
      else
        ys <- rmax - (rmax-rmin)*(1-yloc)
      label2Traces[[i]] <- textTile(labels=ct$label,
                                    x=labx,
                                    y=ys,
                                    col=ct$col,
                                    cex=cex+0.1,
                                    layer=0,
                                    clip=clip,
                                    plot=curplot
      )
      
      # Bin boundaries
      if (together) {
        if (i==1) {
          if (showbins) {
            cb <- unique(vec(ct$bounds))
            bin2Traces[[1]] <- polylinesTile(x=as.numeric(vec(t(cbind(cb,cb)))),
                                             y=rep(c(rmin,rmax), length(cb)),
                                             id=as.numeric(vec(t(cbind(1:length(cb), 1:length(cb))))),
                                             lty="dotted",
                                             col=nicegray,
                                             layer=20,
                                             plot=curplot
            )
          }
        }
      } else {
        if (showbins) {
          cb <- unique(vec(ct$bounds))
          bin2Traces[[i]] <- polylinesTile(x=as.numeric(vec(t(cbind(cb,cb)))),
                                           y=rep(c(rmin,rmax), length(cb)),
                                           id=as.numeric(vec(t(cbind(1:length(cb), 1:length(cb))))),
                                           lty="dotted",
                                           col=nicegray,
                                           layer=20,
                                           plot=curplot
          )
        }
      }
      
    }  
    
    
    # add ROC
    if (any(display=="roc")) {
      
      # Find current column and plot number  
      curcol <- which(display=="roc")
      curplot <- (!together)*(i-1)*length(display) + curcol
      
      # Add ROC curve -- as shaded region or line
      if (together&&(i>1)) {
        bubble3Traces[[i]] <- linesTile(x=ct$roc[,1],
                                        y=ct$roc[,2],
                                        col=ct$col,
                                        layer=10,
                                        plot=curplot
        )
      } else {
        bubble3Traces[[i]] <- polygonTile(x=c(ct$roc[,1],0),
                                          y=c(ct$roc[,2],0),
                                          col=NA,
                                          fill=ct$col,
                                          alpha=alpha,
                                          layer=50,
                                          plot=curplot
        )
        if (border)
          border3Traces[[i]] <- linesTile(x=ct$roc[,1],
                                          y=ct$roc[,2],
                                          col=ct$col,
                                          layer=10,
                                          plot=curplot
          )
      }
      
      # Optionally mark thresholds on ROC
      if (!is.null(thresholds)) {
        xThresh <- yThresh <- rep(NA, length(thresholds))
        for (j in 1:length(thresholds)) {
          threshL <- which(abs(ct$roc[,3]-thresholds[j])==min(abs(ct$roc[,3]-thresholds[j])))
          xThresh[j] <- ct$roc[threshL[1],1]
          yThresh[j] <- ct$roc[threshL[1],2]
        }
        if (together) {
          if (nTraces==2) {
            if (ct$auc==max(c(traceData[[1]]$auc,traceData[[2]]$auc))) {labdirection <- -1} else {labdirection <- 1}                      
          } else {
            labdirection <- (-1)^(i%%2)
          }
        } else {
          labdirection <- -1
        }
        thresh3Trace[[i]] <- pointsTile(x=xThresh, y=yThresh,
                                        col=ct$col, clip=clip, cex=cex-0.15, pch=16,
                                        plot=curplot)
        threshL3Trace[[i]] <- textTile(x=xThresh+0.05*labdirection, y=yThresh, labels=thresholds,
                                       col=ct$col, clip=clip, cex=cex-0.15,
                                       plot=curplot)
      }
      
      # Ref Line 3
      if (together) {
        line3Traces[[1]] <- linesTile(x=c(0,1), y=c(0,1),
                                      layer=15,
                                      col=nicegray,
                                      lty="dashed",
                                      plot=curplot                                   
        )
      } else {
        line3Traces[[i]] <- linesTile(x=c(0,1), y=c(0,1),
                                      layer=15,
                                      col=nicegray,
                                      lty="dashed",
                                      plot=curplot                                   
        )
      }
      
      
      # Model Label 3
      if (together) {yloc <- 1-(i/12)} else {yloc <- 0.95}
      label3Traces[[i]] <- textTile(labels=ct$label,
                                    x=labx,
                                    y=yloc,
                                    col=ct$col,
                                    cex=cex+0.1,
                                    layer=0,
                                    clip=clip,
                                    plot=curplot
      )
      
      # AUC title 3
      if (together) {
        auctitle3Traces[[1]] <- textTile(labels="AUC",
                                         x=0.8,
                                         y=0.5,
                                         col="black",
                                         cex=cex,
                                         layer=0,
                                         clip=clip,
                                         plot=curplot
        )
      } else {
        auctitle3Traces[[i]] <- textTile(labels="AUC",
                                         x=0.8,
                                         y=0.45+(1/12),
                                         col="black",
                                         cex=cex,
                                         layer=0,
                                         clip=clip,
                                         plot=curplot
        )
      }
      
      # AUC label 3 
      if (together) {yloc <- 0.5-(i/12)} else {yloc <- 0.45}
      auc3Traces[[i]] <- textTile(labels=sprintf("%.3f", ct$auc),
                                  x=0.8,
                                  y=yloc,
                                  col=ct$col,
                                  cex=cex,
                                  layer=0,
                                  clip=clip,
                                  plot=curplot
      )
      
    }
  }
  
  
  # Count required traces
  totTraces <- nTraces*6*length(display)
  totTraces <- totTraces + 4*(as.numeric(together) + (!together)*nTraces)
  
  # Collect traces
  collectedTraces <- giantListVector(vector("list", totTraces),
                                     bubbleTraces, bubble2Traces, bubble3Traces,
                                     numberTraces, number2Traces, 
                                     lineTraces, line2Traces, line3Traces,
                                     labelTraces, label2Traces, label3Traces,
                                     auc3Traces, auctitle3Traces, border3Traces, thresh3Trace, threshL3Trace,
                                     binTraces, bin2Traces
  )
  
  # Compute RxC
  RxC <- c((as.numeric(together) + (!together)*nTraces), length(display))
  if (!together&&(length(display)==1)) RxC <- rev(RxC)
  
  # Compute ylog
  ylog <- rep(FALSE, length(display))
  if (any(display=="evp")&&logAxis) ylog[which(display=="evp")] <- TRUE
  if (!together) ylog <- rep(ylog, nTraces)
  
  # Compute limits
  avpLim <- c(-0.01,1.01,-0.01,1.01)
  evpLim <- c(-0.01,1.01,rmin,rmax)
  limits <- NULL
  for (i in 1:length(display)) {
    if ((display[i]=="avp")||(display[i]=="roc"))
      limits <- c(limits, avpLim)
    if (display[i]=="evp")
      limits <- c(limits, evpLim)
  }
  if (!together) limits <- rep(limits, nTraces)
  limits <- matrix(limits, byrow=TRUE, nrow=nPlots, ncol=4)
  
  
  # Compute yaxistitle
  yt1 <- "Actual Frequency (by bin)"
  yt2 <- "Actual/Predicted Probability"
  yt3 <- "True positive rate (sensitivity)"
  yaxistitle <- NULL
  for (i in 1:length(display)) {
    if ((display[i]=="avp")) yaxistitle <- c(yaxistitle, yt1)
    if ((display[i]=="evp")) yaxistitle <- c(yaxistitle, yt2)
    if ((display[i]=="roc")) yaxistitle <- c(yaxistitle, yt3)
  }
  
  # Compute xaxistitle
  xt1 <- xt2 <- "Predicted Probability (binned)"
  xt3 <- "False positive rate (1 - specificity)"
  xaxistitle <- NULL
  for (i in 1:length(display)) {
    if ((display[i]=="avp")) xaxistitle <- c(xaxistitle, xt1)
    if ((display[i]=="evp")) xaxistitle <- c(xaxistitle, xt2)
    if ((display[i]=="roc")) xaxistitle <- c(xaxistitle, xt3)
  }
  
  # Calculate pdf width
  pdfWidth <- RxC[2]*3.25
  
  # Plot traces using tile
  if (is.null(file)) {
    tile(collectedTraces,
         RxC=RxC,
         limits=limits,
         xaxis=list(cex=cex, major=FALSE),
         yaxis=list(cex=cex, major=FALSE, log=ylog),
         xaxistitle=list(labels=xaxistitle),
         yaxistitle=list(labels=yaxistitle),
         height=list(topborder=0.5,bottomborder=0.25,spacer=2.5),
         width=list(spacer=1,null=2.5,yaxistitle=3,yaxis.labelspace=-0.5,
                    leftborder=0.25, rightborder=0.5)
    )
  } else {
    tile(collectedTraces,
         RxC=RxC,
         limits=limits,
         xaxis=list(cex=cex+0.05, major=FALSE),
         yaxis=list(cex=cex+0.05, major=FALSE, log=ylog),
         xaxistitle=list(labels=xaxistitle, cex=cex),
         yaxistitle=list(labels=yaxistitle, cex=cex),
         height=list(topborder=0.25,bottomborder=0.25, xaxistitle=0.5),
         width=list(spacer=1.5,null=5,yaxistitle=0.5,yaxis.labelspace=-0.5,
                    leftborder=0.25, rightborder=0.25),
         output=list(file=file, width=pdfWidth, type="pdf")
    )
  }
}

