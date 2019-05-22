x<-teste
rescea
ptype=4
bk=0.2
n.breaks=100
plotnull=TRUE
avnull=TRUE
nullcol="black"

# function to plot decision matx
ptcea2=function(rescea, ptype=4, bk=0.2, n.breaks=100, plotnull=TRUE, avnull=TRUE, nullcol="black"){
  ptype=4
  if(!is.numeric(bk))bk=0.2
  bk=min(1,max(0.1,bk))
  n.breaks=as.integer(n.breaks)
  if(!is.numeric(ptype) && !is.integer(ptype))ptype=1
  if(!ptype%in%c(1,2,3,4,5))ptype=1
  if(!is.numeric(n.breaks) && !is.integer(n.breaks))n.breaks=100
  if(n.breaks<2)n.breaks=2
  #
  if(ptype<=3){
    null.dist=rescea$null.dist
    decision.dist=as.numeric(rescea$decision.mt)
    decision.neg=density(decision.dist,to=-bk)	
    decision.pos=density(decision.dist,from=bk)	
    maxdeci=max(decision.neg$y, decision.pos$y)
    #
    corr.dist=as.numeric(rescea$corr.mt)
    corr.dist=density(corr.dist)	
    maxcorr=max(corr.dist$y)
    #maxdata=min(c(maxdeci,maxcorr))
    maxdata=maxcorr
  } else {
    null.dist=rescea$null.dist
    decision.dist=as.numeric(rescea$decision.mt)
    #
    decision.hist=hist(decision.dist, plot=FALSE, breaks=n.breaks)
    maxdeci=max(decision.hist$density)
    #
    decision.neg=density(decision.dist,to=-bk)	
    decision.pos=density(decision.dist,from=bk)	
    #
    corr.dist=as.numeric(rescea$corr.mt)
    corr.dist=hist(corr.dist, plot=FALSE, breaks=n.breaks)	
    maxcorr=max(corr.dist$density)
    maxdata=maxcorr		
  }
  #
  if(ptype==1){
    plot(x=c(-1,1),y=c(0,1), col="black",type="n", main="Gene co-expression analysis", 
         xlab="Correlation coefficient", ylab="Density")
    #
    y=apply(null.dist$y,1,mean)
    y=y/max(y)
    x=apply(null.dist$x,1,mean)	
    polygon(x=x,y=y,col="grey", lty="blank")
    #
    x=decision.neg$x		
    y=decision.neg$y/maxdata
    polygon(x=x,y=y, col="blue",lty="blank")
    #
    x=decision.pos$x		
    y=decision.pos$y/maxdata
    polygon(x=x,y=y, col="red",lty="blank")		
    #  
    x=corr.dist$x  	
    y=corr.dist$y/maxdata
    lines(x=x,y=y, col="black", lty="dashed", lwd=2.0)
    #   
    legend("topleft",c("all associations","sig. neg. associations",
                       "sig. pos. association","null distribution"), pch=c(NA_integer_ ,15, 15, 15), 
           pt.cex=1, lty=c(2, 0, 0, 0), merge=TRUE, col=c('black','blue','red','grey'), bty="n", cex=0.7)  
  } 
  if(ptype==2){
    plot(x=c(-1,1),y=c(0,1), col="black",type="n", main="Gene co-expression analysis", 
         xlab="Correlation coefficient", ylab="Density")
    #  
    x=corr.dist$x  	
    y=corr.dist$y/maxdata
    polygon(x=x,y=y, col="grey", lty="blank")
    #
    x=decision.neg$x		
    y=decision.neg$y/maxdata
    polygon(x=x,y=y, col="blue",lty="blank")
    #
    x=decision.pos$x		
    y=decision.pos$y/maxdata
    polygon(x=x,y=y, col="red",lty="blank")
    #	
    if(plotnull){	
      y=null.dist$y/max(null.dist$y)
      x=null.dist$x
      if(avnull){
        x=apply(x,1,mean)
        y=apply(y,1,mean)				
        lines(x,y, lwd=3.0, col=nullcol, lty="dashed")
      } else {
        for(i in 1:min(100,ncol(x)))lines(x=x[,i],y=y[,i],lwd=1.0, col=i, lty="dashed")					
      }		
      legend("topleft",c("all associations","sig. neg. associations",
                         "sig. pos. associations","null distributions"), pch=c(15, 15, 15, NA_integer_), 
             pt.cex=1, lty=c(0, 0, 0, 2), merge=TRUE, col=c('grey','blue','red','black'), bty="n", cex=0.7)
    } else {
      legend("topleft",c("all association","sig. neg. associations",
                         "sig. pos. associations"), pch=c(15, 15, 15), 
             pt.cex=1, lty=c(0, 0, 0), col=c('grey','blue','red'), bty="n", cex=0.7)			
    }
  }
  if(ptype==3){
    plot(x=c(-1,1),y=c(0,1), col="black",type="n", main="Gene co-expression analysis", 
         xlab="Correlation coefficient", ylab="Density")
    #  
    x=corr.dist$x  	
    y=corr.dist$y/maxdata
    polygon(x=x,y=y, col="grey", lty="blank")
    #
    x=decision.neg$x		
    y=decision.neg$y/maxdata
    lines(x=x,y=y, col="blue",lty="dashed", lwd=3.0)
    #
    x=decision.pos$x		
    y=decision.pos$y/maxdata
    lines(x=x,y=y, col="red", lty="dashed", lwd=3.0)
    #
    if(plotnull){	
      y=null.dist$y/max(null.dist$y)
      x=null.dist$x
      if(avnull){
        x=apply(x,1,mean)
        y=apply(y,1,mean)				
        lines(x,y, lwd=3.0, col=nullcol, lty="dashed")
      } else {
        for(i in 1:min(100,ncol(x)))lines(x=x[,i],y=y[,i],lwd=1.0, col=i, lty="dashed")					
      }			
      legend("topleft",c("all associations","sig. neg. associations",
                         "sig. pos. associations","null distributions"), pch=c(15, 15, 15, NA_integer_), 
             pt.cex=1, lty=c(0, 0, 0, 2), merge=TRUE, col=c('grey','blue','red','black'), bty="n", cex=0.7)
    } else {
      legend("topleft",c("all association","sig. neg. associations",
                         "sig. pos. associations"), pch=c(15, 15, 15), 
             pt.cex=1, lty=c(0, 0, 0), col=c('grey','blue','red'), bty="n", cex=0.7)			
    }
  }	
  if(ptype==4){
    plot(x=c(-1,1),y=c(0,1), col="black",type="n", main="Gene co-expression analysis", 
         xlab="Correlation coefficient", ylab="Density")
    #  
    x=corr.dist$breaks  	
    y=corr.dist$density/maxdata	
    rect(x[-length(x)], 0, x[-1], y,col='grey90',border='grey30')
    #	
    xd=decision.hist$breaks		
    yd=decision.hist$density/maxdata
    idx=xd>bk;
    x=c(0,xd[idx]);y=c(0,yd[idx])
    rect(x[-length(x)], 0, x[-1], y,col='red',border='darkred')
    idx=xd<(-bk);
    x=c(xd[idx],0);y=c(yd[idx],0)	
    rect(x[-length(x)], 0, x[-1], y,col='blue',border='darkblue')
    #
    if(plotnull){	
      y=null.dist$y/max(null.dist$y)
      x=null.dist$x
      if(avnull){
        x=apply(x,1,mean)
        y=apply(y,1,mean)
        lines(x,y, lwd=3.0, col=nullcol, lty="dashed")
      } else {
        for(i in 1:min(100,ncol(x)))lines(x=x[,i],y=y[,i],lwd=1.0, col=i, lty="dashed")	
      }			
      legend("topleft",c("all associations","sig. neg. associations",
                         "sig. pos. associations","null distribution"), pch=c(15, 15, 15, NA_integer_), 
             pt.cex=1, lty=c(0, 0, 0, 2), merge=TRUE, col=c('grey','blue','red','black'), bty="n", cex=0.7)
    } else {
      legend("topleft",c("all associations","sig. neg. associations",
                         "sig. pos. associations"), pch=c(15, 15, 15), 
             pt.cex=1, lty=c(0, 0, 0), col=c('grey','blue','red'), bty="n", cex=0.7)			
    }
  }
  if(ptype==5){
    plot(x=c(-1,1),y=c(0,1), col="black",type="n", main="Gene co-expression analysis", 
         xlab="Correlation coefficient", ylab="Density")
    #  
    x=corr.dist$breaks  	
    y=corr.dist$density/maxdata	
    rect(x[-length(x)], 0, x[-1], y,col='grey90',border='grey70')
    #
    x=c(decision.neg$x,0)		
    y=c(decision.neg$y/maxdata,0)
    lines(x=x,y=y, col="blue", lty="dashed", lwd=3.0)
    #
    x=c(0,decision.pos$x)		
    y=c(0,decision.pos$y/maxdata)
    lines(x=x,y=y, col="red", lty="dashed", lwd=3.0)
    #
    if(plotnull){	
      y=null.dist$y/max(null.dist$y)
      x=null.dist$x
      if(avnull){
        x=apply(x,1,mean)
        y=apply(y,1,mean)
        lines(x,y, lwd=3.0, col=nullcol, lty="dashed")
      } else {
        for(i in 1:min(100,ncol(x)))lines(x=x[,i],y=y[,i],lwd=1.0, col=i, lty="dashed")	
      }			
      legend("topleft",c("all associations","sig. neg. associations",
                         "sig. pos. associations","null distribution"), pch=c(15, NA_integer_, NA_integer_, NA_integer_), 
             pt.cex=1, lty=c(0, 2, 2, 2), merge=TRUE, col=c('grey','blue','red','black'), bty="n", cex=0.7)
    } else {
      legend("topleft",c("all associations","sig. neg. associations",
                         "sig. pos. associations"), pch=c(15, NA_integer_, NA_integer_), 
             pt.cex=1, lty=c(0, 2, 2), col=c('grey','blue','red'), bty="n", cex=0.7)			
    }
  }
}

##-----------------------------------------------------------------------------
x<-teste
ceaMCore <- function (x, sig = 0.01, p.adj.method = "fdr", cor.method = "spearman", 
          nper = 1000, plotcea = TRUE, ncore="all", ...) 
{
  library(ff, quietly = TRUE)
  require(doMC)
  if(ncore=="all"){
    ncore = parallel::detectCores()
    registerDoMC(cores = ncore)
  } else{
    registerDoMC(cores = ncore)
  }
  
  
  if (is.data.frame(x)) {
    x = as.matrix(x)
  }else {
    if (!is.matrix(x)) 
      stop("NOTE: not a matrix!")
  }
  cat("Step 1 ...computing correlation", fill = TRUE)
  x = t(x)
  corrMt = cor(x, method = cor.method)
  diag(corrMt) = 0
  uniqueVec = unique(sort(corrMt))
  cat("Step 2 ...computing null distribution", fill = TRUE)
  ctsum = numeric(length(uniqueVec))
  nulldist = list()
  pb = txtProgressBar(min = 0, max = nper, initial = 0, char = "=", 
                      style = 1)
  results <- foreach(i = 1:nper) %dopar%  {
    nulldist2 = list()
    permt = matrix(sample(x), nrow = nrow(x), ncol = ncol(x))
    permt = cor(permt, method = cor.method)
    diag(permt) = NA
    permt = sort(permt)
    ct = findInterval(uniqueVec, permt)
    #ctsum = ctsum + ct
    nl = density(permt)
    nulldist2$x <- nl$x
    nulldist2$y <- nl$y
    nulldist2$size<- length(permt)
    nulldist2$ct<-ct
    #setTxtProgressBar(pb, value = i)
    print(nulldist2)
  }
  #return(results)
  x1<-matrix(nrow = length(results[[1]]$x), ncol = nper)
  x1<-data.frame(x1)
  y1<-matrix(nrow = length(results[[1]]$y), ncol = nper)
  y1<-data.frame(y1)
  for (i in 1:nper) {
    x1[,i]<- as.vector(results[[i]][1])
    y1[,i]<- as.vector(results[[i]][2])
    ctsum<- ctsum + as.vector(results[[i]]$ct)
  }
  nulldist$x<-matrix(x1)
  nulldist$y<-matrix(y1)
  #return(nulldist)
  permt<-results[[length(results)]]
  rm(x1,y1,results)
  cat("\n")
  cat("Step 3 ...computing probs", fill = TRUE)
  probs = ctsum/(length(permt) * nper)
  probs[probs > 0.5] = 1 - probs[probs > 0.5]
  probs = probs[match(as.numeric(corrMt), uniqueVec)]
  cat("Step 4 ...adjusting pvals", fill = TRUE)
  pvalAdj = p.adjust(probs, method = p.adj.method)
  pvalAdj = matrix(pvalAdj, nrow = nrow(corrMt), ncol = nrow(corrMt))
  decision = pvalAdj > sig
  decisionMt = corrMt
  decisionMt[decision] = 0
  rescea = list(corr.mt = corrMt, decision.mt = decisionMt, 
                pvalue.adj = pvalAdj, null.dist = nulldist)
  ptcea = function(rescea, ptype = 4, bk = 0.2, n.breaks = 100, 
                   plotnull = TRUE, avnull = TRUE, nullcol = "black") {
    ptype = 4
    if (!is.numeric(bk)) 
      bk = 0.2
    bk = min(1, max(0.1, bk))
    n.breaks = as.integer(n.breaks)
    if (!is.numeric(ptype) && !is.integer(ptype)) 
      ptype = 1
    if (!ptype %in% c(1, 2, 3, 4, 5)) 
      ptype = 1
    if (!is.numeric(n.breaks) && !is.integer(n.breaks)) 
      n.breaks = 100
    if (n.breaks < 2) 
      n.breaks = 2
    if (ptype <= 3) {
      null.dist = rescea$null.dist
      decision.dist = as.numeric(rescea$decision.mt)
      decision.neg = density(decision.dist, to = -bk)
      decision.pos = density(decision.dist, from = bk)
      maxdeci = max(decision.neg$y, decision.pos$y)
      corr.dist = as.numeric(rescea$corr.mt)
      corr.dist = density(corr.dist)
      maxcorr = max(corr.dist$y)
      maxdata = maxcorr
    }
    else {
      null.dist = rescea$null.dist
      decision.dist = as.numeric(rescea$decision.mt)
      decision.hist = hist(decision.dist, plot = FALSE, 
                           breaks = n.breaks)
      maxdeci = max(decision.hist$density)
      decision.neg = density(decision.dist, to = -bk)
      decision.pos = density(decision.dist, from = bk)
      corr.dist = as.numeric(rescea$corr.mt)
      corr.dist = hist(corr.dist, plot = FALSE, breaks = n.breaks)
      maxcorr = max(corr.dist$density)
      maxdata = maxcorr
    }
    if (ptype == 1) {
      plot(x = c(-1, 1), y = c(0, 1), col = "black", type = "n", 
           main = "Gene co-expression analysis", xlab = "Correlation coefficient", 
           ylab = "Density")
      y = apply(null.dist$y, 1, mean)
      y = y/max(y)
      x = apply(null.dist$x, 1, mean)
      polygon(x = x, y = y, col = "grey", lty = "blank")
      x = decision.neg$x
      y = decision.neg$y/maxdata
      polygon(x = x, y = y, col = "blue", lty = "blank")
      x = decision.pos$x
      y = decision.pos$y/maxdata
      polygon(x = x, y = y, col = "red", lty = "blank")
      x = corr.dist$x
      y = corr.dist$y/maxdata
      lines(x = x, y = y, col = "black", lty = "dashed", 
            lwd = 2)
      legend("topleft", c("all associations", "sig. neg. associations", 
                          "sig. pos. association", "null distribution"), 
             pch = c(NA_integer_, 15, 15, 15), pt.cex = 1, 
             lty = c(2, 0, 0, 0), merge = TRUE, col = c("black", 
                                                        "blue", "red", "grey"), bty = "n", cex = 0.7)
    }
    if (ptype == 2) {
      plot(x = c(-1, 1), y = c(0, 1), col = "black", type = "n", 
           main = "Gene co-expression analysis", xlab = "Correlation coefficient", 
           ylab = "Density")
      x = corr.dist$x
      y = corr.dist$y/maxdata
      polygon(x = x, y = y, col = "grey", lty = "blank")
      x = decision.neg$x
      y = decision.neg$y/maxdata
      polygon(x = x, y = y, col = "blue", lty = "blank")
      x = decision.pos$x
      y = decision.pos$y/maxdata
      polygon(x = x, y = y, col = "red", lty = "blank")
      if (plotnull) {
        y = null.dist$y/max(null.dist$y)
        x = null.dist$x
        if (avnull) {
          x = apply(x, 1, mean)
          y = apply(y, 1, mean)
          lines(x, y, lwd = 3, col = nullcol, lty = "dashed")
        }
        else {
          for (i in 1:min(100, ncol(x))) lines(x = x[, 
                                                     i], y = y[, i], lwd = 1, col = i, lty = "dashed")
        }
        legend("topleft", c("all associations", "sig. neg. associations", 
                            "sig. pos. associations", "null distributions"), 
               pch = c(15, 15, 15, NA_integer_), pt.cex = 1, 
               lty = c(0, 0, 0, 2), merge = TRUE, col = c("grey", 
                                                          "blue", "red", "black"), bty = "n", cex = 0.7)
      }
      else {
        legend("topleft", c("all association", "sig. neg. associations", 
                            "sig. pos. associations"), pch = c(15, 15, 
                                                               15), pt.cex = 1, lty = c(0, 0, 0), col = c("grey", 
                                                                                                          "blue", "red"), bty = "n", cex = 0.7)
      }
    }
    if (ptype == 3) {
      plot(x = c(-1, 1), y = c(0, 1), col = "black", type = "n", 
           main = "Gene co-expression analysis", xlab = "Correlation coefficient", 
           ylab = "Density")
      x = corr.dist$x
      y = corr.dist$y/maxdata
      polygon(x = x, y = y, col = "grey", lty = "blank")
      x = decision.neg$x
      y = decision.neg$y/maxdata
      lines(x = x, y = y, col = "blue", lty = "dashed", 
            lwd = 3)
      x = decision.pos$x
      y = decision.pos$y/maxdata
      lines(x = x, y = y, col = "red", lty = "dashed", 
            lwd = 3)
      if (plotnull) {
        y = null.dist$y/max(null.dist$y)
        x = null.dist$x
        if (avnull) {
          x = apply(x, 1, mean)
          y = apply(y, 1, mean)
          lines(x, y, lwd = 3, col = nullcol, lty = "dashed")
        }
        else {
          for (i in 1:min(100, ncol(x))) lines(x = x[, 
                                                     i], y = y[, i], lwd = 1, col = i, lty = "dashed")
        }
        legend("topleft", c("all associations", "sig. neg. associations", 
                            "sig. pos. associations", "null distributions"), 
               pch = c(15, 15, 15, NA_integer_), pt.cex = 1, 
               lty = c(0, 0, 0, 2), merge = TRUE, col = c("grey", 
                                                          "blue", "red", "black"), bty = "n", cex = 0.7)
      }
      else {
        legend("topleft", c("all association", "sig. neg. associations", 
                            "sig. pos. associations"), pch = c(15, 15, 
                                                               15), pt.cex = 1, lty = c(0, 0, 0), col = c("grey", 
                                                                                                          "blue", "red"), bty = "n", cex = 0.7)
      }
    }
    if (ptype == 4) {
      plot(x = c(-1, 1), y = c(0, 1), col = "black", type = "n", 
           main = "Gene co-expression analysis", xlab = "Correlation coefficient", 
           ylab = "Density")
      x = corr.dist$breaks
      y = corr.dist$density/maxdata
      rect(x[-length(x)], 0, x[-1], y, col = "grey90", 
           border = "grey30")
      xd = decision.hist$breaks
      yd = decision.hist$density/maxdata
      idx = xd > bk
      x = c(0, xd[idx])
      y = c(0, yd[idx])
      rect(x[-length(x)], 0, x[-1], y, col = "red", border = "darkred")
      idx = xd < (-bk)
      x = c(xd[idx], 0)
      y = c(yd[idx], 0)
      rect(x[-length(x)], 0, x[-1], y, col = "blue", border = "darkblue")
      if (plotnull) {
        y = null.dist$y/max(null.dist$y)
        x = null.dist$x
        if (avnull) {
          x = apply(x, 1, mean)
          y = apply(y, 1, mean)
          lines(x, y, lwd = 3, col = nullcol, lty = "dashed")
        }
        else {
          for (i in 1:min(100, ncol(x))) lines(x = x[, 
                                                     i], y = y[, i], lwd = 1, col = i, lty = "dashed")
        }
        legend("topleft", c("all associations", "sig. neg. associations", 
                            "sig. pos. associations", "null distribution"), 
               pch = c(15, 15, 15, NA_integer_), pt.cex = 1, 
               lty = c(0, 0, 0, 2), merge = TRUE, col = c("grey", 
                                                          "blue", "red", "black"), bty = "n", cex = 0.7)
      }
      else {
        legend("topleft", c("all associations", "sig. neg. associations", 
                            "sig. pos. associations"), pch = c(15, 15, 
                                                               15), pt.cex = 1, lty = c(0, 0, 0), col = c("grey", 
                                                                                                          "blue", "red"), bty = "n", cex = 0.7)
      }
    }
    if (ptype == 5) {
      plot(x = c(-1, 1), y = c(0, 1), col = "black", type = "n", 
           main = "Gene co-expression analysis", xlab = "Correlation coefficient", 
           ylab = "Density")
      x = corr.dist$breaks
      y = corr.dist$density/maxdata
      rect(x[-length(x)], 0, x[-1], y, col = "grey90", 
           border = "grey70")
      x = c(decision.neg$x, 0)
      y = c(decision.neg$y/maxdata, 0)
      lines(x = x, y = y, col = "blue", lty = "dashed", 
            lwd = 3)
      x = c(0, decision.pos$x)
      y = c(0, decision.pos$y/maxdata)
      lines(x = x, y = y, col = "red", lty = "dashed", 
            lwd = 3)
      if (plotnull) {
        y = null.dist$y/max(null.dist$y)
        x = null.dist$x
        if (avnull) {
          x = apply(x, 1, mean)
          y = apply(y, 1, mean)
          lines(x, y, lwd = 3, col = nullcol, lty = "dashed")
        }
        else {
          for (i in 1:min(100, ncol(x))) lines(x = x[, 
                                                     i], y = y[, i], lwd = 1, col = i, lty = "dashed")
        }
        legend("topleft", c("all associations", "sig. neg. associations", 
                            "sig. pos. associations", "null distribution"), 
               pch = c(15, NA_integer_, NA_integer_, NA_integer_), 
               pt.cex = 1, lty = c(0, 2, 2, 2), merge = TRUE, 
               col = c("grey", "blue", "red", "black"), bty = "n", 
               cex = 0.7)
      }
      else {
        legend("topleft", c("all associations", "sig. neg. associations", 
                            "sig. pos. associations"), pch = c(15, NA_integer_, 
                                                               NA_integer_), pt.cex = 1, lty = c(0, 2, 2), 
               col = c("grey", "blue", "red"), bty = "n", 
               cex = 0.7)
      }
    }
  }
  if (plotcea) {
    ptcea(rescea)
  }
  return(rescea$decision.mt)
}
