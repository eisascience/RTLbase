

#Mahyari, Eisa: Implementation of Lee et al.'s TL-FC into R
#March, 2017



##### Bias Compensation (Algorithm 4) ##########
###
### Input: hyperplane (w,b), target task data T
##  1: Compute: z_i = <w,x_t_i> +b, for all i
##  2: Build a Grid: s_i = sort (z_i)
##  3: for j = 1 to N_t
##  4: cj = sum(II{abs(z_i - s_j)/||w|| < 1}
##  5: end for
##  6: h = kernelbandwidth ({(s_j,c_j)}_j)
##  7: Smooth: p(z) = sum(c_j k_h (z, s_j))_j
##  8: z∗ = gradientdescent(p(z), 0)
##  9: b_new = b - z∗
##  Output: b_new or f_b(x)= <w,x_i> + b_new
###
############################################################
### Lee, Gyemin, L Stoolman, and C Scott.
### "Transfer Learning for Auto-Gating of Flow Cytometry Data.”
### JMLR(workshop), 2012, 155–66.
############################################################
###
############################################################
###
############################################################
#########
######
###

#' Update classification bias on target tasks
#'
#' @description
#' Implements Algorithm 4 to refine decision boundaries on target tasks by
#' estimating density-informed bias adjustments. Supports optional data
#' manipulations and RCS (rare-cell subset) focusing.
#'
#' @param task_list List of target task feature matrices.
#' @param alg1_result Output list from [alg1_baselineClass()].
#' @param alg2_result Output list from [alg2_rob_meanNCov()].
#' @param alg3_result Vector of bias updates from [alg3_shiftComp()].
#' @param goodColumns Optional column indices to retain.
#' @param save2file Logical; persist diagnostic plots to disk.
#' @param Marg Margin parameter controlling neighborhood width.
#' @param alg4MinFx Mode selection for minima detection (`"gd"`, `"mean"`, `"win"`).
#' @param ADM Logical; apply additional data manipulations.
#' @param useMedian Logical; smooth minima search using medians.
#' @param ZnormMappingBL Logical; use z-normalized mapping for bias updates.
#' @param datatyp Character flag indicating data type (e.g., `"FC"`).
#' @param RCSmodeBL Logical; enable rare-cell subset adjustments.
#' @param RCSfreqSet Numeric vector specifying frequency thresholds for RCS mode.
#' @param CoreClassifier Core classifier identifier (e.g., `"LinSVM"`).
#'
#' @return A list with updated predictions, density summaries, and adjusted
#'   intercepts for each target task.
#' @export
alg4_BiasUpdate <- function(task_list, alg1_result, alg2_result,
                            alg3_result, goodColumns,
                            save2file, Marg, alg4MinFx, ADM=F,
                            useMedian = T, ZnormMappingBL=F, datatyp="FC",
                            RCSmodeBL = F, RCSfreqSet = c(0,0),
                            CoreClassifier = "LinSVM"){



  # #Generic inputs for Dev needs; see kfold code as well
  # task_list = lapply(TestXY.ls, function(x){x$X.test});
  # alg1_result = alg1_res;
  # alg2_result = alg2_res;
  # alg3_result = alg3_res;
  # goodColumns = "";
  # Marg = 1; #[.1; 1]
  # alg4MinFx = "gd"; #mean; gd; win
  # save2file = F; RCSmodeBL = F;CoreClassifier = "LinSVM"; ADM=F; ZnormMappingBL = F; datatyp="FC"





  n_testSets = length(task_list)


  print("starting Alg 4, Bias Update Base Version. No. of test sets:")
  print(n_testSets)


  task_hat_norm_upd <- vector(mode="list")
  task_hat_Z <- vector(mode="list")

  for (m in 1:n_testSets) {
    #if(!exists("m"))   m = 1

    TASK <- as.data.frame(task_list[[m]])
    if(anyNA(TASK)) stop("Task data contains missing values; please clean inputs before bias update.")

    #boxplot(TASK, las=2)
    #violinMyDF(TASK)

    if(ADM) TASK <- as.data.frame(AllDataManipulations(TASK,
                                                       Scale=ScalePerCh,
                                                       Center=MaxPerCh,
                                                       X_cols2Keep = goodColumns,
                                                       globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                       globalRange1010Scale = F,
                                                       globalScaleByVal = F))[,goodColumns]

    if(ncol(TASK)==1){
      TASK <- as.data.frame(cbind(TASK, rep(-1, length(TASK))))
      colnames(TASK) <- c("X", "interc_multp")
    } else {
      TASK$interc_multp <- rep(-1, nrow(TASK))
    }




    if(CoreClassifier == "LinSVM"){
      task_hat_norm_upd[[m]] <- as.matrix(TASK) %*% c(-alg2_result$U_robust[-length(alg2_result$U_robust)],
                                                      alg3_result[m])
    }


    #plot(density(task_hat_norm_upd[[m]]))
    #boxplot(task_hat_norm_upd[[m]])

    tempYhat.orgi <- task_hat_norm_upd[[m]]

    if(RCSmodeBL) {
      xyDF <- as.data.frame(cbind(density(task_hat_norm_upd[[m]], n = nrow(task_hat_norm_upd[[m]]))$x, density(task_hat_norm_upd[[m]], n = nrow(task_hat_norm_upd[[m]]))$y))
      colnames(xyDF) <- c("x", "y")
      #xyDF$x <- order(task_hat_norm_upd[[m]])
      rownames(xyDF) <- rownames(task_hat_norm_upd[[m]])[order(task_hat_norm_upd[[m]])]
      #plot(xyDF)
      DescTools::AUC(x=xyDF$x, y=xyDF$y)

      aprioriMedianFreq <- max(RCSfreqSet)# + mean(rnorm(10000,RCSfreqSet, 3))*2

      #RCS correction, not to over chop to still have landmarks
      aprioriMedianFreq <- ifelse(aprioriMedianFreq < 0.05, 0.05, aprioriMedianFreq)


      thresholds <- seq(min(xyDF$x), max(xyDF$x), length.out = 100)

      AUCthr <- as.numeric(lapply(thresholds, function(xcut){
        tempCutSubDF <- subset(xyDF, x > xcut)
        subAUC <- DescTools::AUC(x=tempCutSubDF$x, y=tempCutSubDF$y)
        #if(subAUC >= aprioriMedianFreq) plot(tempCutSubDF, type="l", lwd=2, col="gold", main=paste("cut @", xcut, "\n Sub AUC = ", subAUC, sep=""))
        subAUC
      }))

      # plot(thresholds, AUCthr, pch=20, col="skyblue")

      #thresholds[max(which(AUCthr >= aprioriMedianFreq))]

      #xyDF <- subset(xyDF, x>=thresholds[max(which(AUCthr >= aprioriMedianFreq))])
      #plot(xyDF, typ="l", main="estimated range by Gauss est. N = 1000")

      #if(class(task_hat_norm_upd[[m]])!="data.frame") task_hat_norm_upd[[m]] <- as.data.frame(task_hat_norm_upd[[m]])
      #colnames(task_hat_norm_upd[[m]]) <- c("zhat")
      #xyDF <- subset(as.data.frame(task_hat_norm_upd[[m]]), zhat >= thresholds[max(which(AUCthr >= aprioriMedianFreq))])
      xyDF.sub <- subset(xyDF, x >= thresholds[max(which(AUCthr >= aprioriMedianFreq))] )
      xyDF.sub$x.true <- task_hat_norm_upd[[m]][rownames(xyDF.sub),]

      if(length(xyDF.sub$y) > 3){
        if (save2file){
          fileID = paste(BaseFigDIR, "/alg4_RCS_DensFocPT_",m, ".png",sep="" )
          png(fileID, width = 1024*2, height = 768*2, units = "px")
          par(mfrow=c(2,2))
          plot(xyDF, xlab="yhat", ylab="density", main=paste("PT thresholded set in red\nMarginal Freq from training: %", round(aprioriMedianFreq*100,3), sep=""), pch=19, col="grey")
          lines(xyDF, col="navy")
          points(xyDF.sub, pch=20, col="red")
          abline(v=thresholds[max(which(AUCthr >= aprioriMedianFreq))], lw=2, lty=3)
          plot(thresholds, AUCthr, pch=20, col="skyblue", main="RCS thresholding", ylab="Marginal Freq.", ylim=range(0,1), xlab="yhat")
          abline(v=thresholds[max(which(AUCthr >= aprioriMedianFreq))], lw=2, lty=3)
          abline(h=aprioriMedianFreq, lw=2, lty=3)
          plot(xyDF.sub[,1:2], main="focused region", col="red", lwd=2, pch=20)
          plot(density(xyDF.sub$x.true), main="density-curve of focused region", col="red", lwd=2)


          dev.off()
        }
        #boxplot(xyDF.sub$x.true)
        task_hat_norm_upd[[m]] <- xyDF.sub$x.true
        remove(tempYhat.orgi)

      } else{
        task_hat_norm_upd[[m]] <- tempYhat.orgi
        remove(tempYhat.orgi)}



      #task_hat_norm_upd[[m]] <- subset(xyDF, x>thresholds[max(which(AUCthr >= aprioriMedianFreq))])$y

    } #######END of RCSmode AUC-based Density FOcusing







    #plot(density(task_hat_norm_upd[[m]] ))
    #z norm for optional testing/future expansion
    #mean(task_hat_norm_upd[[m]])
    #sd(task_hat_norm_upd[[m]])
    z.norm <- (task_hat_norm_upd[[m]]-mean(task_hat_norm_upd[[m]]))/sd(task_hat_norm_upd[[m]]) ## standardized data
    task_hat_Z[[m]] <- z.norm
    #plot(density(z.norm))
    #wilcox.test(rnorm(10000), z.norm)

  }##### END of mapping all test sets

  remove (TASK, m)

  h_bandwidth <- vector()
  c_j <- vector()
  b_updated_alg4 <- vector()

  SVM_b_all <- vector(mode="list")

  for (tn in 1:n_testSets) {
    # if(!exists("tn"))   tn = 1
    print(paste("Starting Bias Update on test set #", tn, sep=""))
    #with alg3 b

    if(!ZnormMappingBL) z_i <- task_hat_norm_upd[[tn]]
    if(ZnormMappingBL) z_i <- task_hat_Z[[tn]]



    #testing approx for classifying RCS especially if(doapprox)
    #worked for scRNASeq but problematic with FC, over smooths,
    #in scRNASeq the low examples do need a boost
    #if(datatyp == "scRNASeqLogNormCount") z_i  <- approx(z_i, n=5000)$y


    #on the range of y_hat
    s_j <- sort(z_i) #s_i ...grid of biases
    names(s_j) <- 1:length(s_j)

    if(save2file){
      # plot(x=s_j, y=(z_i), pch=20, col=factor(sign(z_i))); abline(h=0); abline(v=0)
    }




    w_t_euc_mag <- alg2_result$w_euc_mag #euclidean norm of base

    #this needs to be redone for speedup

    FigSet <- unique(c(round(1:9/10*length(s_j)),
                       sample(round(6:9/10.1*length(s_j)), replace=F, 4),
                       sample(round(6:9/10.2*length(s_j)), replace=F, 4)))
    if(length(FigSet)>13) FigSet <- sort(c(1, sample(FigSet,13,replace=F), length(s_j)))




    GridX = 1:length(s_j)

    #needs speeding up
    for (j in GridX){
      if(!exists("j"))   j=5900

      #abs(bla)/||w_t||

      densNorm <- abs(z_i - s_j[j])/w_t_euc_mag


      c_j[j] <- sum(as.numeric(lapply((densNorm < Marg), IndicatorFX)))


    }
    # if(save2file) dev.off()

    names(c_j) <- 1:length(c_j)
    print("counting within margin complete")

    #testing phase figs
    if(save2file){



    }

    # ShortLoop = 0
    # while(ShortLoop<3){
    #   ShortLoop = ShortLoop + 1



    h_bandwidth[tn] <- round(KBand_fx(s_k=s_j, c_k=c_j), 9)
    #h_bandwidth[tn] <- KBand_fx(s_k=s_j, c_k=c_j)

    OneD_OptiBW <- density(s_j, n=length(s_j))$bw
    MeanBW = mean(c( h_bandwidth[tn], OneD_OptiBW))


    #############Smooth the counts, with a defined minimum length

    #list of x and y
    Gaus_Ker_Smooth_sj_cj = ksmooth(s_j,c_j,kernel="normal", h_bandwidth[tn],
                                    n.points = max(c(1000L, length(s_j))))



    ###alternate bw
    Gaus_Ker_Smooth_sj_cj.optbw = ksmooth(s_j,c_j,kernel="normal", OneD_OptiBW,
                                          n.points = max(c(1000L, length(s_j))))
    Gaus_Ker_Smooth_sj_cj.meanbw = ksmooth(s_j,c_j,kernel="normal", MeanBW,
                                           n.points = max(c(1000L, length(s_j))))





    Gaus_Ker_Smooth_sj_cj.DF         <- as.data.frame(Gaus_Ker_Smooth_sj_cj)
    Gaus_Ker_Smooth_sj_cj.optbw.DF   <- as.data.frame(Gaus_Ker_Smooth_sj_cj.optbw)
    Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- as.data.frame(Gaus_Ker_Smooth_sj_cj.meanbw)


    Gaus_Ker_Smooth_sj_cj.DF        <- Gaus_Ker_Smooth_sj_cj.DF[complete.cases(Gaus_Ker_Smooth_sj_cj.DF), ]
    Gaus_Ker_Smooth_sj_cj.optbw.DF  <- Gaus_Ker_Smooth_sj_cj.optbw.DF[complete.cases(Gaus_Ker_Smooth_sj_cj.optbw.DF), ]
    Gaus_Ker_Smooth_sj_cj.meanbw.DF <- Gaus_Ker_Smooth_sj_cj.meanbw.DF[complete.cases(Gaus_Ker_Smooth_sj_cj.meanbw.DF), ]


    # Remove 0s
    Gaus_Ker_Smooth_sj_cj.DF         <- Gaus_Ker_Smooth_sj_cj.DF[ min( which ( Gaus_Ker_Smooth_sj_cj.DF$y != 0 )):max( which( Gaus_Ker_Smooth_sj_cj.DF$y != 0 )), ]
    Gaus_Ker_Smooth_sj_cj.optbw.DF   <- Gaus_Ker_Smooth_sj_cj.optbw.DF[ min( which ( Gaus_Ker_Smooth_sj_cj.optbw.DF$y != 0 )):max( which( Gaus_Ker_Smooth_sj_cj.optbw.DF$y != 0 )), ]
    Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- Gaus_Ker_Smooth_sj_cj.meanbw.DF[ min( which ( Gaus_Ker_Smooth_sj_cj.meanbw.DF$y != 0 )):max( which( Gaus_Ker_Smooth_sj_cj.meanbw.DF$y != 0 )), ]

    if(datatyp == "FC"){
      #1% threshold
      #MaxMinDT = 10
      MinDensThreshold = round(max(Gaus_Ker_Smooth_sj_cj.DF$y)*ifelse(max(Gaus_Ker_Smooth_sj_cj.DF$y)>100, 0.1, 0.01))


      Gaus_Ker_Smooth_sj_cj.DF         <- subset(Gaus_Ker_Smooth_sj_cj.DF, y > MinDensThreshold)
      Gaus_Ker_Smooth_sj_cj.optbw.DF   <- subset(Gaus_Ker_Smooth_sj_cj.optbw.DF, y > MinDensThreshold)
      Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- subset(Gaus_Ker_Smooth_sj_cj.meanbw.DF, y > MinDensThreshold)


    }
    if(datatyp == "scRNASeqLogNormCount"){

      #1% threshold
      if(RCSmodeBL) {
        MinDensThreshold = 1 #round( max(Gaus_Ker_Smooth_sj_cj.DF$y)*ifelse(max(Gaus_Ker_Smooth_sj_cj.DF$y)>1000, 0.1, 0.01))

      } else {
        MinDensThreshold = round( max(Gaus_Ker_Smooth_sj_cj.DF$y)*ifelse(max(Gaus_Ker_Smooth_sj_cj.DF$y)>1000, 0.1, 0.05))

      }
      #MaxMinDT = 1

      #if the data represents the line y=1, the entire dataset gets removed and causes NAs
      if(nrow(subset(Gaus_Ker_Smooth_sj_cj.DF, y > MinDensThreshold))>2) Gaus_Ker_Smooth_sj_cj.DF         <- subset(Gaus_Ker_Smooth_sj_cj.DF, y > MinDensThreshold)
      if(nrow(subset(Gaus_Ker_Smooth_sj_cj.optbw.DF, y > MinDensThreshold))>2) Gaus_Ker_Smooth_sj_cj.optbw.DF   <- subset(Gaus_Ker_Smooth_sj_cj.optbw.DF, y > MinDensThreshold)
      if(nrow(subset(Gaus_Ker_Smooth_sj_cj.meanbw.DF, y > MinDensThreshold))>2) Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- subset(Gaus_Ker_Smooth_sj_cj.meanbw.DF, y > MinDensThreshold)


    }




    Gaus_Ker_Smooth_sj_cj.DF        <- try(as.data.frame(approx(Gaus_Ker_Smooth_sj_cj.DF, n=nrow(Gaus_Ker_Smooth_sj_cj.DF))))
    Gaus_Ker_Smooth_sj_cj.optbw.DF  <- try(as.data.frame(approx(Gaus_Ker_Smooth_sj_cj.optbw.DF, n=nrow(Gaus_Ker_Smooth_sj_cj.optbw.DF))))
    Gaus_Ker_Smooth_sj_cj.meanbw.DF <- try(as.data.frame(approx(Gaus_Ker_Smooth_sj_cj.meanbw.DF, n=nrow(Gaus_Ker_Smooth_sj_cj.meanbw.DF))))




    if(class(Gaus_Ker_Smooth_sj_cj.DF)=="try-error"){
      if(class(Gaus_Ker_Smooth_sj_cj.optbw.DF)=="try-error"){
        if(class(Gaus_Ker_Smooth_sj_cj.meanbw.DF)=="try-error"){
          #all are error so get original DF
          Gaus_Ker_Smooth_sj_cj.DF         <- as.data.frame(Gaus_Ker_Smooth_sj_cj)
        } else{
          Gaus_Ker_Smooth_sj_cj.DF         <- Gaus_Ker_Smooth_sj_cj.meanbw.DF
        }
      } else {

        Gaus_Ker_Smooth_sj_cj.DF         <- Gaus_Ker_Smooth_sj_cj.optbw.DF

      }

    }



    # print(paste("Freq of rows lost due to smoothing: %", round((length(Gaus_Ker_Smooth_sj_cj$x) - nrow(Gaus_Ker_Smooth_sj_cj.DF))/length(Gaus_Ker_Smooth_sj_cj$x)*100,4), sep=""))

    SMmin.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.DF, print2screen = T, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.00001), mode=ifelse(RCSmodeBL, "MedianAllMin", ifelse(datatyp=="FC","MedianAll", "MedianAll")) ), silent = T)


    #try 3 times to get it
    if(class(SMmin.try)=="try-error"){
      SMmin.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.DF, print2screen = F, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.0001), mode=ifelse(RCSmodeBL, "LeftOfPeakNPavg", ifelse(datatyp=="FC","MedianAll", "MedianAll")) ), silent = T)
      if(class(SMmin.try)=="try-error"){
        SMmin.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.DF, print2screen = F, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.0001), mode=ifelse(RCSmodeBL, "NP", ifelse(datatyp=="FC","NP", "NP")) ), silent = T)
      }
    }

    #if still error, try else
    if(class(SMmin.try)=="try-error"){

      SMmin.meanbw.try  <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.meanbw.DF, print2screen = T, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.00001), mode=ifelse(RCSmodeBL, "LeftOfPeakNPavg", ifelse(datatyp=="FC","NPZAvg", "NP")) ), silent = T)

      if(class(SMmin.meanbw.try)=="try-error"){

        SMmin.optbw.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.optbw.DF, print2screen = F, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.00001), mode=ifelse(RCSmodeBL, "LeftOfPeak", ifelse(datatyp=="FC","NPZAvg", "NPZAvg")) ), silent = T)


        if(class(SMmin.optbw.try)=="try-error"){
          ###No GD Minima Found regardless of bw adjustments
          #minimaErr = T
          #heuristic attempt
          SMmin.heur <- try(Gaus_Ker_Smooth_sj_cj.DF[find_peaks(-1*Gaus_Ker_Smooth_sj_cj.DF$y),])

          if(class(SMmin.heur)=="try-error"){
            #find a place in teh center at least :(
            SMmin.heur <- Gaus_Ker_Smooth_sj_cj.DF[which.min(abs(Gaus_Ker_Smooth_sj_cj.DF$x - median(Gaus_Ker_Smooth_sj_cj.DF$x))),]
          } else {
            if(nrow(SMmin.heur)==0) {
              #find a place in teh center at least :(
              SMmin.heur <- Gaus_Ker_Smooth_sj_cj.DF[which.min(abs(Gaus_Ker_Smooth_sj_cj.DF$x - median(Gaus_Ker_Smooth_sj_cj.DF$x))),]
            }
          }

          SMmin.ls <- list(PeakA = Gaus_Ker_Smooth_sj_cj.DF[which.max(Gaus_Ker_Smooth_sj_cj.DF$y),],
                           PeakB = Gaus_Ker_Smooth_sj_cj.DF[which.max(Gaus_Ker_Smooth_sj_cj.DF$y),],
                           Minima = SMmin.heur)

        } else{
          #SMmin.optbw.try
          SMmin.ls <- SMmin.optbw.try
          #ShortLoop = 10 #end the loop
        }
      } else {
        #SMmin.meanbw.try
        SMmin.ls <- SMmin.meanbw.try
        #ShortLoop = 10 #end the loop
      }
    } else {
      #SMmin.try
      SMmin.ls <- SMmin.try
      #ShortLoop = 10 #end the loop
    }


    #}#end of while(ShortLoop)



    #SMmin.ls

    if(save2file){
      fileID = paste(BaseFigDIR, "/alg4_KeyOptima_", names(task_list)[tn], "_", round(h_bandwidth[tn],4), ".png",sep="" )
      png(fileID, width = 1024, height = 768, units = "px")

      plot(as.data.frame(cbind(s_j, c_j)), typ="l", lwd=2, col="grey")
      points((Gaus_Ker_Smooth_sj_cj.DF), col="red", pch=19)
      #points((Gaus_Ker_Smooth_sj_cj.optbw.DF), col="dodgerblue", pch=19)
      #points((Gaus_Ker_Smooth_sj_cj.meanbw.DF), col="plum", pch=19)
      abline(h=MinDensThreshold, lty=3, col="navy")

      lines(Gaus_Ker_Smooth_sj_cj.DF, type="l", col="dodgerblue", lwd=2,
            main=paste("Key Local Optima on test sample: ",names(task_list)[tn] , sep=""),
            xlim=range(min(Gaus_Ker_Smooth_sj_cj$x),max(Gaus_Ker_Smooth_sj_cj$x)));
      #abline(v=alg1_result$baselineSVM[tn,2], col="red", lwd=2);
      abline(v=-alg3_result[tn], col="dodgerblue", lwd=2)
      abline(v=-alg2_result$U_robust["b.int"], col="forestgreen", lty=2, lwd=2);
      abline(v=SMmin.ls$PeakA, col="orange", lwd=2);
      abline(v=SMmin.ls$PeakB, col="orange", lwd=2);
      #abline(v=MinimaS_j$x, col="cyan", lwd=2, lty=3)
      abline(v=SMmin.ls$Minima$x, col="purple", lwd=2, lty=3)


      legend("topleft",
             c("alg2_robustMean", "alg3_shifted", "Peak Positions", "Found Local minima", "Counts above 1% of total"),
             lwd=c(2,2,2,2,2),
             col=c("forestgreen", "dodgerblue","orange","purple", "red"),
             lty=c(2,2,2,3,2))

      dev.off()
    }


    s_j_KeyOptima <- c(LocMin = SMmin.ls$Minima, peakYB = SMmin.ls$PeakB, peakYA = SMmin.ls$PeakA )
    names(s_j_KeyOptima) <- c("LocMin","peakYB","peakYA")

    b_updated_alg4[tn] <- s_j_KeyOptima[c("LocMin")]


    par(mfrow=c(1,1))

    #SVM_b_all[[tn]] <- cbind(alg1_result$baselineSVM[tn,ncol(alg1_result$baselineSVM)], alg3_b_upd=alg3_result[tn], alg4_b_upd = (b_updated_alg4[tn]));
    SVM_b_all[[tn]] <- cbind(alg3_b_upd=alg3_result[tn], alg4_b_upd = (b_updated_alg4[tn]));

  }

  SVM_b_all.df <- rbindlist(lapply(SVM_b_all, as.data.frame))

  #correction based on position
  SVM_b_all.df$alg4_b_upd <- as.numeric(SVM_b_all.df$alg3_b_upd) + as.numeric(SVM_b_all.df$alg4_b_upd)

  SVM_b_all.df <- as.data.frame( cbind(rep(alg2_result$U_robust[length(alg2_result$U_robust)],
                                           nrow(SVM_b_all.df)), SVM_b_all.df))
  colnames(SVM_b_all.df) <- c("b_alg2_norm", "b_alg3_norm", "b_alg_norm")
  rownames(SVM_b_all.df) <- paste("Task", 1:nrow(SVM_b_all.df), sep="")

  return(SVM_b_all.df)


}



