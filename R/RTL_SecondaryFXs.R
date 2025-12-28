
# Gram-Schmidt orthonormalization, vectors in rows
# implemented from an online matlab version of the code

orthnormal <- function(X)
{
  #X <- cbind(as.numeric(alg2_res$v_0[-length(alg2_res$v_0)]), as.numeric(alg2_res$U_robust[-length(alg2_res$U_robust)]))
  X<-as.matrix(X)
  n<-nrow(X)
  p<-ncol(X)

  W<-NULL
  if(p > 1) {
    W<-cbind(W, X[,1])
    for(k in 2:p) {
      #k = 2
      gw<-rep(0, n)

      for(i in 1:(k-1)) {
        #i=1
        gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
        gw<-gw + gki * W[,i]
      }
      W<-cbind(W, X[,k] - gw)
    }
  } else {
    W<-cbind(W, X[,1])
  }

  W <- replace(W, is.na(W), 0)


  W <- - apply(W, 2, function(x) norm(x, type="2"))
  W
}


#' Approximate minima by scanning mirrored windows
#'
#' @description
#' Uses paired sliding windows over an `x/y` density grid to approximate local
#' minima, returning the most representative minimum across windows.
#'
#' @param xyDF Data frame with `x` and `y` columns representing a density curve.
#' @param NWindows Integer number of window pairs to evaluate (must be even).
#'
#' @return Numeric position of the selected minimum.
#' @export
ClosingWindowAproxFunOptimizedMinima <- function(xyDF, NWindows = 50){
  #xyDF <- SmoothedXY.t[,c("x","y")]

  if (is.odd(NWindows)){
    NWindows <- NWindows + 1
    print("closing window, set length.out as even to pair up")
  }

  Windows <- seq(min(xyDF$x), max(xyDF$x), length.out = NWindows)
  Windows <- cbind(Windows[seq(1,NWindows/2)], Windows[seq(NWindows,(NWindows/2+1))])


  Minimas <- lapply(1:nrow(Windows), function(Win){
    xyDF.op <- optimize(approxfun(xyDF$x,xyDF$y),interval=as.numeric(Windows[Win,]))
    xyDF.op$minimum
  })

  Minimas[sapply(Minimas, is.null)] <- NULL

  Minimas <- as.numeric(Minimas)
  Minimas <- Minimas[(nrow(Windows)-5):nrow(Windows)]


  as.numeric(Minimas[which(round(as.numeric(unlist(Minimas))) == Mode(round(as.numeric(unlist(Minimas)))))])
  mean(as.numeric(Minimas[which(round(as.numeric(unlist(Minimas))) == Mode(round(as.numeric(unlist(Minimas)))))]))
}

#' Gradient-descent helper for identifying minima in 2D curves
#'
#' @description
#' Fits a linear model via stochastic average gradient descent over a subset of
#' the supplied density grid and returns the point of minimal deviation between
#' the model and observed curve.
#'
#' @param xyDF Data frame with `x` and `y` columns.
#' @param quants Optional quantile bounds used to trim the search region.
#' @param MinMaxRange Optional numeric range restricting `x` values.
#' @param learningRate Learning rate passed to `gradDescent::SAGD`.
#' @param maximumIters4GD Maximum number of iterations for gradient descent.
#' @param print2screen Logical; print optimization progress.
#' @param plot2screen Logical; produce diagnostic plots.
#'
#' @return A list with the optimized minimum, fitted coefficients, and
#'   preprocessing parameters.
#' @export
SpecializedGD2D <- function(xyDF,
                            quants            = c(0.1, 0.9),
                            MinMaxRange       = F,
                            learningRate      = 0.1,
                            maximumIters4GD   = 100,
                            print2screen = F, plot2screen = F){
  #library(gradDescent)

  # print2screen = F, plot2screen = F
  # x1 <- unique(round(c(rnorm(1000, 1), rnorm(600, 6)),3))
  # x1Dens <- density(x1, n=length(x1))
  # xyDF <- round(as.data.frame(cbind(x=x1Dens$x, y=x1Dens$y)),3)
  # plot(xyDF, type="l", lwd=2, col="plum");
  # remove(x1, x1Dens)
  # quants = c(0.1, 0.9) #F
  # MinMaxRange= F
  # learningRate = 0.001
  # maximumIters4GD = 100000
  # xyDF <- SmoothedXY.t.biPeaks

  #xyDF <- SmoothedXY.t
  #quants = F; MinMaxRange = c(-5,6); learningRate = 0.01; maximumIters4GD = 100000; print2screen = T; plot2screen = T


  # the result of gd is the coeffs i.e., slope m and intercept b
  f.hat.Lin = function(x, BM){
    x*BM[2] + BM[1]
  }


  xyDF.or <- xyDF

  if(class(quants) == "numeric"){
    if(length(quants) == 2) xyDF <- subset(subset(xyDF, x > quantile(xyDF$x, min(quants))), x < quantile(xyDF$x, max(quants)))
  } else {
    if(!quants){
      #quants F so no subset
      #if(print2screen) print("no quants in GD")
    } else {
      #err
      print("error @ quants")
    }
  }
  #  plot(xyDF, type="l", lwd=2, col="plum");


  if(class(MinMaxRange) == "numeric"){
    if(length(MinMaxRange) == 2) xyDF <- subset(subset(xyDF, x > min(MinMaxRange)), x< max(MinMaxRange))
  } else {
    if(!MinMaxRange){
      #MinMaxRange F so no subset
      #print("no minmax range")
    } else {
      #err
      print("error @ quants")
    }
  }
  #  plot(xyDF, type="l", lwd=2, col="plum");

  ## build GD model. Choosing the Stochastic average GD here heuristically performs better on tested data
  GDmodel <- gradDescent::SAGD(xyDF, alpha = learningRate, maxIter = maximumIters4GD);
  if(print2screen) print(GDmodel)


  f1.hat.Lin <- f.hat.Lin(xyDF$x, GDmodel)
  f1.MinLin <- as.data.frame(cbind(x=xyDF$x, y=f1.hat.Lin))

  if(plot2screen){
    plot(xyDF.or, type='l', col="plum", lwd=2); abline(v=3.6, lwd=2, lty=2, col="grey")
    abline(GDmodel, lty=2, col="powderblue", lwd=2)
    lines(f1.MinLin, col="red", lwd=2, lty=3)
  }



  xyDF$DistL2Y <- abs(round(xyDF$y - f1.MinLin$y))
  xyDF <- xyDF[order(xyDF$DistL2Y),]

  if(print2screen) head(xyDF)

  PossibleMins <- xyDF[which( xyDF$DistL2Y == min(xyDF$DistL2Y)),]
  temp.i = 0
  while(nrow(PossibleMins) < 4){
    temp.i = temp.i + 1
    PossibleMins <- rbind(PossibleMins,
                          xyDF[which( xyDF$DistL2Y == min(xyDF$DistL2Y) + temp.i ),])
    if(temp.i > 10) break
  }

  if(min(PossibleMins$x) == max(PossibleMins$x)) {
    PossibleMinsXmin <- optimize(approxfun(x=PossibleMins$x, y=PossibleMins$y),
                                 interval = range(min(PossibleMins$x) - 0.01,
                                                  max(PossibleMins$x) + 0.01))$minimum

  } else {
    PossibleMinsXmin <- optimize(approxfun(x=PossibleMins$x, y=PossibleMins$y),
                                 interval = range(min(PossibleMins$x),
                                                  max(PossibleMins$x)))$minimum
  }


  if(plot2screen) abline(v=PossibleMinsXmin)

  return(list(GDoptimMin = PossibleMinsXmin,
              allPossibleMins = PossibleMins,
              GDcoefs = GDmodel,
              call = list(quants            = quants,
                          MinMaxRange       = MinMaxRange,
                          learningRate      = learningRate,
                          maximumIters4GD   = maximumIters4GD,
                          xyDF = xyDF)))


}

#' Locate local maxima in a numeric series
#'
#' @param x Numeric vector to search.
#' @param m Neighborhood size for evaluating peak prominence.
#'
#' @return Integer indices of detected peaks.
#' @export
find_peaks <- function (x, m = 4){
  #https://github.com/stas-g/findPeaks
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#' Identify minima for bias correction curves
#'
#' @description
#' Combines gradient-descent searches, windowed minima, and peak/dip detection
#' to select a representative minimum from smoothed density estimates. The
#' method supports several selection modes to tailor behaviour to data type.
#'
#' @param SmoothedXY Data frame of smoothed `x/y` coordinates.
#' @param print2screen Logical; print progress updates.
#' @param datatyp Character flag indicating data type (e.g., `"FC"`).
#' @param learnRate Learning rate passed to gradient-descent helper routines.
#' @param mode Character specifying the minima selection strategy (e.g., `"NZ"`,
#'   `"MedianAll"`).
#' @param cleanQuantsX Logical; trim extreme `x` quantiles before searching.
#' @param cleanQuantsY Logical; trim low `y` density values before searching.
#'
#' @return A list containing selected peaks and minima positions.
#' @export
SmartMinimaAk <- function(SmoothedXY, print2screen=F, datatyp="FC", learnRate=0.1, mode="NZ", cleanQuantsX=T, cleanQuantsY=F){

  # cleanQuantsX = T; cleanQuantsY = F; print2screen = T; learnRate = 0.0000001; mode="MedianAll"
  # SmoothedXY <- Gaus_Ker_Smooth_sj_cj.DF

  # SmoothedXY <- Gaus_Ker_Smooth_ak_ck.DF
  #SmoothedXY <- xyDF
  #SmoothedXY <- as.data.frame(Gaus_Ker_Smooth_diff2)
  #SmoothedXY <- Gaus_Ker_Smooth_ak_ck.DF[,c(1,4)]*c(rep(1,nrow(Gaus_Ker_Smooth_ak_ck.DF)), rep(-1,nrow(Gaus_Ker_Smooth_ak_ck.DF)))

  #SmoothedXY <- s_j.varBWkernSmooth.ls[[smoothDFN]][,1:2]


  SmoothedXY.t <- SmoothedXY[complete.cases(SmoothedXY), ]
  #plot(SmoothedXY.t)


  if(cleanQuantsX) SmoothedXY.t <- subset(subset(SmoothedXY.t, x > quantile(SmoothedXY.t$x, .1)), x < quantile(SmoothedXY.t$x, .9))

  if(print2screen) plot(SmoothedXY.t, type ="l", col="plum", lwd=2, lty=3)

  ThreshY <- quantile(SmoothedXY.t$y, .1)*.1
  # plot(SmoothedXY.t, type="l", lwd=2, col="plum"); abline(h=ThreshY, lty=3)
  if(cleanQuantsY) SmoothedXY.t <- SmoothedXY.t[which(SmoothedXY.t$y > ThreshY),]




  # three different ways to get critical optima
  allPeaks <- try(SmoothedXY.t[find_peaks(SmoothedXY.t$y, m=round(length(SmoothedXY.t$y)/200)),])
  #abline(v=allPeaks$x)


  if(class(allPeaks) == "try-error" || nrow(allPeaks)==0 )  {
    allPeaks <- try(SmoothedXY.t[find_peaks(SmoothedXY.t$y, m=1),])
  }

  if(class(allPeaks) != "try-error")  {
    allPeaks <-rbind(allPeaks, SmoothedXY.t[unlist(CriticalPointsX(SmoothedXY.t$y, n=3, plot2screen=F)$top),])
  } else{
    allPeaks <- SmoothedXY.t[unlist(CriticalPointsX(SmoothedXY.t$y, n=2, plot2screen=F)$top),]
  }

  allPeaksV3 <- argmax(SmoothedXY.t$x, SmoothedXY.t$y, 1)$i
  if(class(allPeaksV3) != "try-error")  {
    allPeaks <- rbind(allPeaks, SmoothedXY.t[allPeaksV3,])
  }; remove(allPeaksV3)


  #
  #
  #   if(nrow(allPeaks) <= 5)  {
  #
  #
  #     allPeaks <- rbind(allPeaks,  SmoothedXY.t[as.numeric(lapply(quantile(SmoothedXY.t$x, c(.1, .4, .6, .9)), function(quantV){
  #       which.min(abs(SmoothedXY.t$x - quantV))
  #     })),])
  #
  #
  #   #allPeaks <- rbind(allPeaks, SmoothedXY.t[round(seq(2, nrow(SmoothedXY.t), length.out = (5 - nrow(allPeaks)))),])
  #
  #   }




  SmoothedXinvY <- as.data.frame(cbind(x=SmoothedXY.t$x, y=-1 * SmoothedXY.t$y))
  allDips <- try(SmoothedXinvY[find_peaks(SmoothedXinvY$y, m=2),])


  # plot(SmoothedXY.t)
  # abline(v=allPeaks$x)


  allPeaks <- allPeaks[order(allPeaks$x),]
  allPeaks <- allPeaks[rownames(unique(allPeaks)),]


  if(nrow(allPeaks) >= 2) {
    allPeaks.or <- allPeaks
  } else {
    allPeaks <- as.data.frame(rbind(allPeaks, SmoothedXY.t[which.max(SmoothedXY.t$y)+1, ]))
    allPeaks.or <- allPeaks
  }

  PeakFound <- allPeaks.or[which.max(allPeaks.or$y),]
  allPeaks <- allPeaks.or[-which.max(allPeaks.or$y),]
  while(nrow(allPeaks)>4){
    allPeaks <- allPeaks[-which.min(allPeaks$y),]
  }

  # plot(SmoothedXY.t)
  # abline(v=allPeaks$x)



  #add the min and max for windowed approach from outside to inside
  allPeaks <- rbind(SmoothedXY.t[which.min(SmoothedXY.t$x),], allPeaks, SmoothedXY.t[which.max(SmoothedXY.t$x),])
  #if(is.odd(nrow(allPeaks)))    allPeaks <- rbind(SmoothedXY.t[round(length(SmoothedXY.t$x)/2),], allPeaks)
  allPeaks <- allPeaks[order(allPeaks$x),]

  if(nrow(allPeaks)<4){
    #add a few points in the center
    allPeaks <-  as.data.frame(rbind(allPeaks, SmoothedXY.t[round(nrow(SmoothedXY.t)/2)-2,]))
    allPeaks <- as.data.frame(rbind(allPeaks, SmoothedXY.t[round(nrow(SmoothedXY.t)/2)+2,]))
  }



  #middle point to base the windows
  #windowBase <- allPeaks[round(length(allPeaks$x)/2),]

  #windowBase <- allPeaks[which.max(allPeaks$y),]

  #allPeaks <- allPeaks[-which(allPeaks$x == windowBase$x),]

  #print(paste("Doing GD xTimes windows:", nrow(allPeaks)/2))
  print(paste("Doing GD xTimes windows:", nrow(allPeaks)))
  #if(nrow(nrow(allPeaks)))
  #1:(nrow(allPeaks))/2
  GDls <- lapply(1:(nrow(allPeaks)), function(rowN){
    #rowN = 1
    #learnRate = 0.0000001
    #print(rowN)


    #MMR = c(allPeaks[rowN, ]$x, allPeaks[nrow(allPeaks) - rowN+1, ]$x)
    MMR = c(allPeaks[rowN, ]$x, PeakFound$x)
    #plot(SmoothedXY.t)
    #abline(v=MMR)

    if(nrow(subset(subset(SmoothedXY.t, x>min(MMR)), x<max(MMR)))> 2){

      sp2DGD.res.ls <- try(SpecializedGD2D(xyDF = SmoothedXY.t,
                                                quants = F,
                                                MinMaxRange= c(min(MMR),max(MMR)),
                                                learningRate = learnRate,
                                                maximumIters4GD = 100000,
                                                print2screen = T, plot2screen = T), silent = T); #sp2DGD.res.ls$GDoptimMin
      #abline(v=windowBase$x, col="gold", lwd=2)

      #if(class(sp2DGD.res.ls) == "try-error") sp2DGD.res.ls <- NA

      CWAFOM <-try(ClosingWindowAproxFunOptimizedMinima(subset(subset(SmoothedXY.t, x>min(MMR)), x<max(MMR))), silent = T); CWAFOM

      # SmoothedXinvY <- as.data.frame(cbind(x=subset(subset(SmoothedXY.t, x>min(MMR)), x<max(MMR))$x, y=-1 * subset(subset(SmoothedXY.t, x>min(MMR)), x<max(MMR))$y))

      # invFindPeak <- try(SmoothedXinvY[find_peaks(SmoothedXinvY$y, m=2),])

      #print(CWAFOM)

      if(class(CWAFOM)!="try-error") {
        if(class(sp2DGD.res.ls) == "try-error") {
          sp2DGD.res.ls <- list(GDoptimMin = CWAFOM)
          # if(class(invFindPeak)!="try-error" && nrow(invFindPeak)!=0 ){
          #   sp2DGD.res.ls$GDoptimMin <- (invFindPeak$x + sp2DGD.res.ls$GDoptimMin)/2
          # }

        } else{
          sp2DGD.res.ls$GDoptimMin <- mean(c(CWAFOM, sp2DGD.res.ls$GDoptimMin))
          # if(class(invFindPeak)!="try-error" && nrow(invFindPeak)!=0){
          #   sp2DGD.res.ls$GDoptimMin <- (invFindPeak$x + sp2DGD.res.ls$GDoptimMin)/2
          # }
        }
      }


    } else {
      sp2DGD.res.ls <- list(GDoptimMin = NULL)
    }



    sp2DGD.res.ls$GDoptimMin

  })
  GDls <- lapply(GDls, function(xN){
    if(!is.null(xN)) xN[[1]] else NA
  })

  if(length(GDls[sapply(GDls, is.na)])>0) GDls[sapply(GDls, is.na)] <- NULL







  MinimasFound <- SmoothedXY.t[unlist(lapply(1:length(GDls), function(GDlsN){
    #GDlsN=1
    if(class(GDls[[GDlsN]])=="list") GDls[[GDlsN]] = unlist(GDls[[GDlsN]])
    which.min(abs(SmoothedXY.t$x - GDls[[GDlsN]]))
  })),]

  # if(class(allDips)!="try-error"){
  #   allDips$y <- -1*allDips$y
  #   colnames(allDips) <- colnames(MinimasFound)
  #   MinimasFound <- rbind(MinimasFound, allDips)
  #   MinimasFound <- MinimasFound[order(MinimasFound$x),]
  # }

  #plot(SmoothedXY.t, typ="l", col="dodgerblue")
  #abline(v=MinimasFound$x, col="red", lwd=2, lty=2)

  MinimasFound <- as.data.frame(MinimasFound[which(MinimasFound$y <= quantile(MinimasFound$y, .8)),])
  MinimasFound <- MinimasFound[rownames(unique(MinimasFound)),]

  #PeakA <- allPeaks[which.max(allPeaks$y),]
  PeakA <- PeakFound
  allPeaks.or$dist2PeakA <- abs(allPeaks.or$x - PeakA$x)

  if(nrow(subset(allPeaks.or, dist2PeakA>1))>2) {
    tempPeakSet <- subset(allPeaks.or, dist2PeakA>1)
    PeakB <- tempPeakSet[-which.max(tempPeakSet$y),][which.max(tempPeakSet[-which.max(tempPeakSet$y),]$x),]
    remove(tempPeakSet)
  } else{
    if(nrow(subset(allPeaks.or, dist2PeakA>0.1))>2) {
      tempPeakSet <- subset(allPeaks.or, dist2PeakA>0.1)
      PeakB <- tempPeakSet[-which.max(tempPeakSet$y),][which.max(tempPeakSet[-which.max(tempPeakSet$y),]$x),]
      remove(tempPeakSet)
    } else{
      if(nrow(subset(allPeaks.or, dist2PeakA>0.05))>2) {
        tempPeakSet <- subset(allPeaks.or, dist2PeakA>0.05)
        PeakB <- tempPeakSet[-which.max(tempPeakSet$y),][which.max(tempPeakSet[-which.max(tempPeakSet$y),]$x),]
        remove(tempPeakSet)
      } else PeakB <- allPeaks.or[-which.max(allPeaks.or$y),][which.max(allPeaks.or[-which.max(allPeaks.or$y),]$x),]
    }
  }



  MinimasFound.or <- MinimasFound
  # if(PeakA$x > PeakB$x) {
  #   #OK if the Pos Pop freq < 0.5
  #   if(length(which(MinimasFound$x > PeakB$x))>1) MinimasFound <- MinimasFound[which(MinimasFound$x > PeakB$x),]
  # } else {
  #   #OK if the Pos Pop freq < 0.5
  #   if(length(which(MinimasFound$x > PeakA$x))>1) MinimasFound <- MinimasFound[which(MinimasFound$x > PeakA$x),]
  # }



  if(mode=="LeftOfPeak") {
    #Refers to peak B when bimodal and A if unimodal
    #Minimum to the left of peak
    if(nrow(subset(MinimasFound.or, x<PeakA$x))>0) {
      chosenMinima <- rbind(subset(MinimasFound.or, x<PeakA$x), subset(MinimasFound.or, x<PeakB$x))
      chosenMinima <- chosenMinima[which.min(chosenMinima$y),]

    } else {
      chosenMinima <- subset(MinimasFound.or, x<PeakB$x)
      chosenMinima <- chosenMinima[which.min(chosenMinima$y),]
    }
  }

  if(mode=="LeftOfPeakAonly") {
    #Refers to peak B when bimodal and A if unimodal
    #Minimum to the left of peak
    if(nrow(subset(MinimasFound.or, x<PeakA$x))>0) {
      chosenMinima <- rbind(subset(MinimasFound.or, x<PeakA$x))
      chosenMinima <- chosenMinima[which.min(chosenMinima$y),]

    } else {
      chosenMinima <- SmoothedXY.t[which.min(SmoothedXY.t$x),]
      chosenMinima$x <- mean(c(PeakA$x, chosenMinima$x))


    }
  }
  #closest to 0
  if(mode=="NZ") chosenMinima <- MinimasFound.or[which.min(abs(MinimasFound.or$x - 0)) ,]
  #smallest y
  if(mode=="Min") chosenMinima <- MinimasFound.or[which.min(MinimasFound.or$y) ,]
  #closest to peak
  if(mode=="NP") chosenMinima <- MinimasFound[which.min(abs(MinimasFound$x - PeakA$x)),]

  if(mode=="NPZMinAvg"){
    chosenMinima <- (MinimasFound.or[which.min(abs(MinimasFound.or$x - 0)) ,] +
                       MinimasFound[which.min(abs(MinimasFound$x - PeakFound$x)),] +
                       MinimasFound[which.min(MinimasFound$y) ,])/3
    # chosenMinima <- MinimasFound[which.min(order( (abs(MinimasFound$x - 0) +
    #                    order(abs(MinimasFound$x - PeakFound$x)) +
    #                    order(MinimasFound$y))/3 )),]
  }
  if(mode=="NPZAvg"){
    chosenMinima <- (MinimasFound.or[which.min(abs(MinimasFound.or$x - 0)) ,] +
                       MinimasFound[which.min(abs(MinimasFound$x - PeakFound$x)),])/2
  }
  if(mode=="MinNPZAvg"){
    chosenMinima <- (MinimasFound.or[which.min(abs(MinimasFound.or$x - 0)) ,] +
                       MinimasFound[which.min(abs(MinimasFound$x - PeakFound$x)),] +
                       MinimasFound[which.min(MinimasFound$y) ,])/3
  }

  if(mode=="NPMinAvg"){
    chosenMinima <- (MinimasFound.or[which.min(MinimasFound.or$y) ,] +
                       MinimasFound[which.min(abs(MinimasFound$x - PeakFound$x)),])/2
  }

  if(mode=="MedianAll"){
    chosenMinima <- as.data.frame(t(colMeans(MinimasFound.or)))
    chosenMinima$x <- median(MinimasFound.or$x)
  }
  if(mode=="MeanAll"){
    chosenMinima <- as.data.frame(t(colMeans(MinimasFound.or)))
  }

  if(mode=="MedianAllMinNP"){
    chosenMinima <- as.data.frame(t(colMeans(MinimasFound.or)))
    chosenMinima$x <- median(MinimasFound.or$x)

    chosenMinima <- (chosenMinima + MinimasFound.or[which.min(MinimasFound.or$y) ,] +
                       MinimasFound[which.min(abs(MinimasFound$x - PeakFound$x)),])/3
  }
  if(mode=="MedianAllMin"){
    chosenMinima <- as.data.frame(t(colMeans(MinimasFound.or)))
    chosenMinima$x <- median(MinimasFound.or$x)
    chosenMinima <- (chosenMinima + MinimasFound.or[which.min(MinimasFound.or$y) ,])/2
  }
  if(mode=="MedianAllNZ"){
    chosenMinima <- as.data.frame(t(colMeans(MinimasFound.or)))
    chosenMinima$x <- median(MinimasFound.or$x)
    chosenMinima <- (chosenMinima + MinimasFound.or[which.min(abs(MinimasFound.or$x - 0)) ,])/2
  }

  if(mode=="LeftOfPeakNPavg") {
    #Refers to peak B when bimodal and A if unimodal
    #Minimum to the left of peak
    if(nrow(subset(MinimasFound.or, x<PeakA$x))>0) {
      chosenMinima <- rbind(subset(MinimasFound.or, x<PeakA$x), subset(MinimasFound.or, x<PeakB$x))
      chosenMinima <- chosenMinima[which.min(chosenMinima$y),]

    } else {
      chosenMinima <- subset(MinimasFound.or, x<PeakB$x)
      chosenMinima <- chosenMinima[which.min(chosenMinima$y),]
    }

    chosenMinima <- (chosenMinima + MinimasFound[which.min(abs(MinimasFound$x - PeakA$x)),])/2

  }


  # if(nrow(MinimasFound)>1){
  #   MinimasFound <- MinimasFound[which(MinimasFound$x == chosenMinima$x),]
  #   chosenMinima <- rbind(chosenMinima, MinimasFound[which.min(abs(MinimasFound$x - PeakFound$x)),])
  #   chosenMinima <- chosenMinima[which.min(chosenMinima$y),]
  # }





  #if(chosenMinima$x < PeakA$x) PeakB <- allPeaks.or[which(allPeaks.or$x < chosenMinima$x),][which.min(allPeaks.or[which(allPeaks.or$x < chosenMinima$x),]$x),]
  #if(chosenMinima$x > PeakA$x) PeakB <- allPeaks.or[which(allPeaks.or$x > chosenMinima$x),][which.min(allPeaks.or[which(allPeaks.or$x > chosenMinima$x),]$x),]

  # plot(SmoothedXY.t, typ="l", col="dodgerblue")
  # abline(v=c(PeakA$x, PeakB$x, chosenMinima$x), col=c("plum", "plum", "red"))


  #return(list(PeakA=PeakA$x, PeakB=PeakB$x, Minima=sp2DGD.res.ls$GDoptimMin, Minimawin=ModeWindowsMinima))
  return(list(PeakA=PeakA$x, PeakB=PeakB$x, Minima=chosenMinima))

}




#' Find inflection points in a numeric sequence
#'
#' @param InputNumVec Numeric vector to evaluate.
#' @param n Integer specifying the sensitivity threshold.
#' @param plot2screen Logical; plot detected minima/maxima.
#'
#' @return A list with `bottom` and `top` indices corresponding to minima and
#'   maxima.
#' @export
CriticalPointsX<- function(InputNumVec, n = 3, plot2screen=T){



  #InputNumVec <- 100 + cumsum(rnorm(50, 0.2, 1)) # climbs upwards most of the time
  bottoms <- lapply(1:n, function(x) inflect(InputNumVec, threshold = x)$minima)

  tops <- lapply(1:n, function(x) inflect(InputNumVec, threshold = x)$maxima)
  # Color functions
  cf.1 <- grDevices::colorRampPalette(c("pink","red"))
  cf.2 <- grDevices::colorRampPalette(c("cyan","blue"))
  if(plot2screen){
    plot(InputNumVec, type = 'l', main = "Minima & Maxima\nVariable Thresholds")
    for(i in 1:n){
      points(bottoms[[i]], InputNumVec[bottoms[[i]]], pch = 16, col = cf.1(n)[i], cex = 1)
    }
    for(i in 1:n){
      points(tops[[i]], InputNumVec[tops[[i]]], pch = 16, col = cf.2(n)[i], cex = 1)
    }
    legend("topleft", legend = c("Minima",1:n,"Maxima",1:n),
           pch = rep(c(NA, rep(16,n)), 2), col = c(1, cf.1(n),1, cf.2(n)),
           pt.cex =  c(rep(c(1, c(1:n) / 1.5), 2)), cex = .75, ncol = 2)

  }


  return(list(bottom = unlist(bottoms[[n]]), top = unlist(tops[[n]])))

}

#' Compute minima and maxima indices using inflection thresholds
#'
#' @param x Numeric vector to analyze.
#' @param threshold Integer window size for detecting changes in slope.
#'
#' @return A list with `minima` and `maxima` index vectors.
#' @export
inflect <- function(x, threshold = 1){
  up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
  down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
  a    <- cbind(x,up,down)
  list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
}


#' Locate local maxima after LOESS smoothing
#'
#' @param x Numeric vector of x-values.
#' @param y Numeric vector of y-values.
#' @param w Window size supplied to `zoo::rollapply`.
#' @param ... Additional arguments passed to `stats::loess`.
#'
#' @return A list with the x-positions of maxima, their indices, and smoothed
#'   fitted values.
#' @export
argmax <- function(x, y, w=1, ...) {
  #xSmoothedXY.t$y #SmoothedXY.t$y
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}
