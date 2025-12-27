





#Mahyari, Eisa: Implementation of Lee et al.'s TL-FC into R
#March, 2017



##### Normal Vector Update (Algorithm 6) ##########
###
###
############################################################
### Lee, Gyemin, L Stoolman, and C Scott.
### "Transfer Learning for Auto-Gating of Flow Cytometry Data.”
### JMLR(workshop), 2012, 155–66.
############################################################
###
### Input: hyperplane (w,b), target task data T, and
###          direction of change vt
###
### for a_k = -0.5 to 0.5 step 0.01 do
###        w_k <- w + ak vt
###        ck = sum(II{abs(<w,x_t> + b/||w||) < 1}
### end for
### h <- Kernel bandwidth({(ak,cK)}k)
### Smooth: g(a) <- sum(ck kh(a, ak))
###  at <- gradient decent (g(a), 0)
###  w_new <- w + at vt
###
###
### Output: w_new
############################################################
###
############################################################
#########
######
###



#' Update normal vectors using Algorithm 6
#'
#' @description
#' Rotates the decision boundary normal vectors toward directions that maximize
#' density support on target tasks. The routine evaluates candidate steps along
#' principal directions and returns updated hyperplanes.
#'
#' @param task_list List of target task feature matrices.
#' @param alg1_result Output list from [alg1_baselineClass()].
#' @param alg2_result Output list from [alg2_rob_meanNCov()].
#' @param alg3_result Vector of bias updates from [alg3_shiftComp()].
#' @param alg4_result Data frame of bias-adjusted intercepts from
#'   [alg4_BiasUpdate()].
#' @param X_feat_cols Optional feature indices to retain.
#' @param Marg Margin parameter controlling density evaluation width.
#' @param save2file Logical; save diagnostic figures.
#' @param ADM Logical; apply additional data manipulations.
#' @param datatyp Character flag indicating data type (e.g., `"FC"`).
#' @param RCSmodeBL Logical; enable rare-cell subset adjustments.
#' @param RCSfreqSet Numeric vector specifying frequency thresholds for RCS mode.
#' @param CoreClassifier Core classifier identifier (e.g., `"LinSVM"`).
#'
#' @return A list containing updated normal vectors, smoothing metadata, and key
#'   optima explored during the search.
#' @export
alg6_NormalVectorUpdate <- function(task_list, alg1_result,
                                    alg2_result, alg3_result,
                                    alg4_result, X_feat_cols,
                                    Marg, save2file, ADM, datatyp="FC",
                                    RCSmodeBL = F, RCSfreqSet = c(0,0),
                                    CoreClassifier = "LinSVM"){

  # task_list = lapply(TrainTest.ls$TestSetXYls, function(x){x$X.test});
  # if(!exists("alg1_result"))   alg1_result = alg1_res;
  # if(!exists("alg2_result"))   alg2_result = alg2_res;
  # if(!exists("(alg3_result"))  alg3_result = alg3_res;
  # if(!exists("alg4_result"))   alg4_result = alg4_res
  # X_feat_cols = "";
  # save2file = F;
  # Marg = .2; ADM=F; datatyp="FC";
  # RCSmodeBL = F


  # if(!exists("alg1_result"))   alg1_result = alg1_res;
  # if(!exists("alg2_result"))   alg2_result = alg2_res;
  # if(!exists("(alg3_result"))  alg3_result = alg3_res;
  # if(!exists("alg4_result"))   alg4_result = alg4_res
  # if(!exists("X_feat_cols"))   X_feat_cols = ImpFeats
  # if(!exists("task_list"))     task_list  = TestXls.t; save2file = F



  n_testSets = length(task_list)

  print("starting Alg 6, Norm. Vec. Update")



  #inputs
  #w_t_org <- alg2_result$U_robust_norm[-length(alg2_result$U_robust_norm)]

  w_t_org <- alg2_result$U_robust[-length(alg2_result$U_robust)]

  #1st PC for direction of change
  v_t_org <- alg2_result$v_t[1]

  #initializations
  alg_6_h_bandwidth <- vector()
  alg6_slope <- vector()
  task_hat_norm_upd <- vector(mode="list")
  alg6_w_new <- matrix(0, ncol=length(w_t_org), nrow=n_testSets)
  w_k_euc_mag_vec <- vector()
  KeyOptima <- vector(mode="list")

  for (m in 1:n_testSets) {
    # m=1

    b_t_alg4 <- alg4_result[m,"b_alg_norm"]#/alg2_result$w_euc_mag

    print(paste(" robust mean hyperplane with alg4's updated bias for task #",m,":",sep=""));
    #print(b_t_alg4)
    if(class(b_t_alg4)=="list") b_t_alg4 <- unlist(b_t_alg4)



    print(paste("...NormVecUpdate", m, sep=" "))

    TASK <- as.data.frame(task_list[[m]])


    #violinMyDF(TASK)
    if(ADM) TASK <- as.data.frame(AllDataManipulations(TASK,
                                                       Scale=ScalePerCh,
                                                       Center=MaxPerCh,
                                                       X_cols2Keep = X_cols2Keep,
                                                       globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                       globalRange1010Scale = F,
                                                       globalScaleByVal = F))[,X_feat_cols]

    CountInMarg = T
    iterCount = 0
    while(CountInMarg){
      iterCount = iterCount + 1

      #c(min(TASK),max(TASK))
      if(iterCount == 1){
        # a_k <- c(rnorm(2000, 0, 3), rnorm(1000, 5, 3), rnorm(1000, -5, 3))
        # a_k <- unique(sort(round(unique(c(a_k, 0,
        #                                   seq(from = -0.5, to = 0.5, length.out = 1000),
        #                                   seq(from = min(TASK)/4, to = max(TASK)/4, length.out = 1000),
        #                                   seq(from = min(TASK)/2, to = max(TASK)/2, length.out = 1000),
        #                                   seq(from = min(TASK), to = max(TASK), length.out = 1000))),3)))
        a_k <- sort(unique(round(c(seq(-0.55,0.55, length.out = 400), seq(-0.09, 0.09, length.out = 500), seq(-0.002, 0, length.out = 400), seq(0, 0.002, length.out = 400), seq(-0.009, 0.009, length.out = 400), seq(-0.000009, 0.000009, length.out = 400)),9) ))


        #FC tends to have larger samples sizes so reducing the number of point to quantify to speed up the process
        if(datatyp=="FC"){
          #a_k <- seq(min(a_k), max(a_k), length.out = round((abs(max(a_k)) + abs(min(a_k))+1)*1.5))
          #a_k <- sort(unique(c(sample(a_k, 100, replace=F), seq(-0.5,0.5, length.out = 50))))
          a_k <- sort(unique(round(c(seq(-0.55,0.55, length.out = 400), seq(-0.09, 0.09, length.out = 500), seq(-0.002, 0, length.out = 400), seq(0, 0.002, length.out = 400), seq(-0.009, 0.009, length.out = 400), seq(-0.000009, 0.000009, length.out = 400)),9) ))

        }# in round 3,4 try something different



      }



      if(iterCount %in% c(2, 3)){
        a_k <- unique(sort(c(rnorm(2000, 0, 2), rnorm(1000, 5, 2), rnorm(1000, -5, 2))))
        if(datatyp=="FC") a_k <- sort(unique(round(c(seq(-0.55,0.55, length.out = 500), seq(-0.009, 0.009, length.out = 400), seq(-0.00009, 0.00009, length.out = 400), seq(-0.00000009, 0.00000009, length.out = 400)),9) ))




      }
      if(iterCount %in% c(4, 5,6)){
        a_k <- unique(sort(c(rnorm(2000, 0, 1.5), rnorm(2000, 3, 1.5), rnorm(2000, -3, 1.5))))
        a_k <- unique(sort(c(a_k,  asinh(a_k))))
      }


      if(!(0 %in% a_k)) a_k <- sort(c(0, a_k))
      #a_k[round(length(a_k)/2)+1] = 0


      #plot(density(a_k))

      if(ncol(TASK)==1){
        TASK <- as.data.frame(cbind(TASK, rep(ifelse(CoreClassifier=="LinSVM", -1, 1), length(TASK))))
        colnames(TASK) <- c("X", "interc_multp")
      } else {
        TASK$interc_multp <- rep(ifelse(CoreClassifier=="LinSVM", -1, 1), nrow(TASK))
      }




      c_k <- vector()

      printSet <- sort(c(1, sample(unique(c(round(1:9/10*length(a_k)),
                                            sample(round(1:10/10*length(a_k)/4), replace=F, 7)))
                                   ,13,replace=F), length(a_k)))




      #length(printSet)

      if(save2file){
        fileID = paste(BaseFigDIR, "/alg6_SampShifts_",names(task_list)[m], ".png",sep="" )
        png(fileID, width = 1024*2, height = 768*2, units = "px")
        par(mfrow=c(5,6))
      }

      for (aks in 1:length(a_k)){
        # aks = which(a_k==0)

        if(CoreClassifier=="LinSVM") w_k <- -(w_t_org + a_k[aks]*v_t_org)
        if(CoreClassifier=="GLM")    w_k <- 1 * round(w_t_org + a_k[aks]*v_t_org/alg2_result$w_euc_mag,3)

        if(is.na(as.numeric(b_t_alg4))) b_t_alg4 <- unlist(b_t_alg4)

        if(CoreClassifier=="LinSVM") w_b_k <- as.numeric(c(w_k, as.numeric(b_t_alg4)))
        if(CoreClassifier=="GLM")    w_b_k <- as.numeric(c(w_k, -as.numeric(b_t_alg4)))
        #w_b_k[41] = -2
        z_hat_temp <- as.matrix(TASK) %*% (w_b_k);

        # plot(density(z_hat_temp))
        # temp.yhat <- as.data.frame(z_hat_temp)
        # temp.yhat$yorg <- Kfold.ValidationSets.ls$kfold1$y[[1]]
        # temp.yhat$dy <- density(temp.yhat$V1, n=nrow(temp.yhat))$y
        #
        # temp.yhat <- temp.yhat[order(temp.yhat$V1), ]
        # points(temp.yhat$V1, temp.yhat$dy, col=factor(temp.yhat$yorg), pch=19)

        w_k_euc_mag <- norm(w_k, type="2") #euclidean norm
        z_hat_temp_norm <- abs(z_hat_temp/w_k_euc_mag)



        # plot(as.data.frame(cbind(z_hat_temp, z_hat_temp_norm)))

        # plot(density(z_hat_temp_norm))

        #if(CoreClassifier=="GLM") z_hat_temp_norm <- scale(abs(z_hat_temp/w_k_euc_mag))

        if(save2file){
          if(aks %in% printSet){
            #fileID = paste(BaseFigDIR, "/alg6_SampShifts_",names(task_list)[m], ".png",sep="" )
            #png(fileID, width = 1024, height = 768, units = "px")

            # plot(density(sign(z_hat_temp)), main=paste("density z_hat\na_k = ",a_k[aks],sep=""))
            # abline(v=0, col="grey", lty=2, lwd=2)
            # plot(density(z_hat_temp_norm), main=paste("density z_hat_norm\na_k = ",a_k[aks],sep=""), xlim=range(-1,20))
            # abline(v=Marg, col="red", lty=2, lwd=2)

          }}

        c_k[aks] <- sum(as.numeric(lapply((z_hat_temp_norm < Marg), IndicatorFX))) +1

      }
      if(save2file){
        dev.off()
      }
      par(mfrow=c(1,1))

      # plot(cbind(x_a_k=a_k, y_c_k=c_k), type="l", xlim=range(-1,2)); abline(v=0, col="red")
      # abline(h=MinDensThreshold)

      if(datatyp == "FC"){
        MinDensThreshold = round( max(c_k)*0.1)
        if(nrow(subset(as.data.frame(cbind(x_a_k=a_k, y_c_k=c_k)), y_c_k >= MinDensThreshold))<10){
          MinDensThreshold = round( max(c_k)*0.1)
          if(nrow(subset(as.data.frame(cbind(x_a_k=a_k, y_c_k=c_k)), y_c_k >= MinDensThreshold))<10){
            MinDensThreshold = 2
            if(nrow(subset(as.data.frame(cbind(x_a_k=a_k, y_c_k=c_k)), y_c_k >= MinDensThreshold))<10){
              MinDensThreshold = 0
            }
          }
        }

        Gaus_Ker_Smooth_ak_ck.DF  <- subset(as.data.frame(cbind(x_a_k=a_k, y_c_k=c_k)), y_c_k >= MinDensThreshold)

      }

      if(datatyp == "scRNASeqLogNormCount"){
        MinDensThreshold = round( max(c_k)*ifelse(RCSmodeBL, 0.01, 0.1))

        Gaus_Ker_Smooth_ak_ck.DF         <- subset(as.data.frame(cbind(x_a_k=a_k, y_c_k=c_k)), y_c_k >= MinDensThreshold)
      }

      Gaus_Ker_Smooth_ak_ck.DF <- Gaus_Ker_Smooth_ak_ck.DF[complete.cases(Gaus_Ker_Smooth_ak_ck.DF), ]

      #plot(Gaus_Ker_Smooth_ak_ck.DF, typ="l", col="red", lwd=2)

      if(nrow(Gaus_Ker_Smooth_ak_ck.DF)>=3) CountInMarg = F
      if(iterCount == 6) {
        CountInMarg = F #Cant find an appropriate density post count Find 0 as solx; hack to prevent crash
        print("Could not find appropriate density for the counts; 0 solx given")
        Gaus_Ker_Smooth_ak_ck.DF <- as.data.frame(c(rnorm(1000,-2,1), rnorm(1000,2,1)))
        colnames(Gaus_Ker_Smooth_ak_ck.DF) <- c("x_a_k")
        Gaus_Ker_Smooth_ak_ck.DF$y_c_k <- density(Gaus_Ker_Smooth_ak_ck.DF$x_a_k, n=2000)$y
        Gaus_Ker_Smooth_ak_ck.DF$x_a_k <- density(Gaus_Ker_Smooth_ak_ck.DF$x_a_k, n=2000)$x
        Gaus_Ker_Smooth_ak_ck.DF <- subset(subset(Gaus_Ker_Smooth_ak_ck.DF, x_a_k >-4),  x_a_k <4)
        plot(Gaus_Ker_Smooth_ak_ck.DF)
      }
    }

    a_k = Gaus_Ker_Smooth_ak_ck.DF$x_a_k
    c_k = Gaus_Ker_Smooth_ak_ck.DF$y_c_k

    if(nrow(Gaus_Ker_Smooth_ak_ck.DF)>3){
      #JUMP! :)
      alg_6_h_bandwidth[m] <- KBand_fx(s_k=a_k, c_k=c_k)
      OneD_OptiBW <- density(a_k, n=length(a_k))$bw

      #plot(as.data.frame(cbind(a_k, c_k)))
      Gaus_Ker_Smooth_ak_ck = ksmooth(a_k,c_k, kernel="normal", bandwidth=alg_6_h_bandwidth[m], n.points = length(a_k));
      Gaus_Ker_Smooth_ak_ck.DF <- as.data.frame(Gaus_Ker_Smooth_ak_ck)
      Gaus_Ker_Smooth_ak_ck.DF <- Gaus_Ker_Smooth_ak_ck.DF[complete.cases(Gaus_Ker_Smooth_ak_ck.DF), ]
      #plot(Gaus_Ker_Smooth_ak_ck.DF, type="l")

      repOptiBW =T
      iterCount = 0
      while(repOptiBW){
        iterCount = iterCount +1
        Gaus_Ker_Smooth_ak_ck.OptiBW = ksmooth(a_k,c_k, kernel="normal", bandwidth=OneD_OptiBW, n.points = length(a_k));
        Gaus_Ker_Smooth_ak_ck.DF.OptiBW <- as.data.frame(Gaus_Ker_Smooth_ak_ck.OptiBW)
        Gaus_Ker_Smooth_ak_ck.DF.OptiBW <- Gaus_Ker_Smooth_ak_ck.DF[complete.cases(Gaus_Ker_Smooth_ak_ck.DF.OptiBW), ]
        #plot(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, type="l")
        if(nrow(Gaus_Ker_Smooth_ak_ck.DF.OptiBW)>2) repOptiBW = F else OneD_OptiBW <- OneD_OptiBW/5
        if(iterCount>10) {
          Gaus_Ker_Smooth_ak_ck.DF.OptiBW <- Gaus_Ker_Smooth_ak_ck.DF
          repOptiBW = F
        }
      }


      if(datatyp!="FC"){
        Gaus_Ker_Smooth_ak_ck.DF <- as.data.frame(approx(Gaus_Ker_Smooth_ak_ck.DF, n=max(c(750, nrow(Gaus_Ker_Smooth_ak_ck.DF)))))
        Gaus_Ker_Smooth_ak_ck.DF.OptiBW <- as.data.frame(approx(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, n=max(c(750, nrow(Gaus_Ker_Smooth_ak_ck.DF)))))

      } else{

        Gaus_Ker_Smooth_ak_ck.DF <- as.data.frame(approx(Gaus_Ker_Smooth_ak_ck.DF, n=max(c(750, nrow(Gaus_Ker_Smooth_ak_ck.DF)))))
        Gaus_Ker_Smooth_ak_ck.DF.OptiBW <- as.data.frame(approx(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, n=max(c(750, nrow(Gaus_Ker_Smooth_ak_ck.DF)))))

        CountThresh <- 3

        if(nrow(subset(Gaus_Ker_Smooth_ak_ck.DF, y>CountThresh))>10){
          Gaus_Ker_Smooth_ak_ck.DF <- subset(Gaus_Ker_Smooth_ak_ck.DF, y>CountThresh)
          Gaus_Ker_Smooth_ak_ck.DF.OptiBW <- subset(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, y>CountThresh)
          # plot(cbind(x_a_k=a_k, y_c_k=c_k), type="l"); abline(v=0, col="red")
          # points(Gaus_Ker_Smooth_ak_ck.DF, pch=20, col="red")
          # points(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, pch=20, col="dodgerblue")
          # abline(h=CountThresh, lwd=2)
        } else
        {
          CountThresh <- 1 #round(max(Gaus_Ker_Smooth_ak_ck.DF$y)*ifelse(max(Gaus_Ker_Smooth_ak_ck.DF$y)>1000, 0.1, 0.01))
          Gaus_Ker_Smooth_ak_ck.DF <- subset(Gaus_Ker_Smooth_ak_ck.DF, y>CountThresh)
          Gaus_Ker_Smooth_ak_ck.DF.OptiBW <- subset(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, y>CountThresh)
        }
      }



      # library(RTL)
      # tempDen <- density(c(rnorm(1000, -1), rnorm(500, 3), rnorm(600, 7), rnorm(1000, 12)))
      # Gaus_Ker_Smooth_ak_ck.DF <- as.data.frame(cbind(tempDen$x, tempDen$y))
      # colnames(Gaus_Ker_Smooth_ak_ck.DF) <- c("x", "y")
      # plot(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, type="l")
      #plot(Gaus_Ker_Smooth_ak_ck.DF, type="l")

      # print(paste("Freq of rows lost due to smoothing: %", round((length(Gaus_Ker_Smooth_sj_cj$x) - nrow(Gaus_Ker_Smooth_sj_cj.DF))/length(Gaus_Ker_Smooth_sj_cj$x)*100,4), sep=""))



      SMmin.try         <- try(SmartMinimaAk(Gaus_Ker_Smooth_ak_ck.DF, learnRate=ifelse(CoreClassifier=="GLM", 0.001, ifelse(datatyp=="FC", 0.0000001, 0.1)), print2screen = T, mode=ifelse(datatyp=="FC", "MedianAllNZ", "MedianAllNZ")), silent = T)
      MinimaErr = F

      if(class(SMmin.try)=="try-error"){
        SMmin.OptiBW.try  <- try(SmartMinimaAk(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, learnRate=ifelse(CoreClassifier=="GLM", 0.001, ifelse(datatyp=="FC", 0.000001, 0.1)), print2screen = T, mode=ifelse(datatyp=="FC", "MedianAllNZ", "MedianAllNZ")), silent = T)

        if(class(SMmin.OptiBW.try)=="try-error"){
          ###No GD Minima Found regardless of bw adjustments
          SMmin.heur <- try(Gaus_Ker_Smooth_ak_ck.DF[find_peaks(-1*Gaus_Ker_Smooth_ak_ck.DF$y),])

          if(class(SMmin.heur)=="try-error"){
            #Out of options :(
            SMmin.ls <- list(PeakA = Gaus_Ker_Smooth_ak_ck.DF[which.max(Gaus_Ker_Smooth_ak_ck.DF$x),]$x,
                             PeakB = Gaus_Ker_Smooth_ak_ck.DF[which.min(Gaus_Ker_Smooth_ak_ck.DF$x),]$x,
                             Minima = Gaus_Ker_Smooth_ak_ck.DF[which.min(Gaus_Ker_Smooth_ak_ck.DF$y),])
          } else {
            SMmin.ls <- list(PeakA = Gaus_Ker_Smooth_ak_ck.DF[which.max(Gaus_Ker_Smooth_ak_ck.DF$x),]$x,
                             PeakB = Gaus_Ker_Smooth_ak_ck.DF[which.min(Gaus_Ker_Smooth_ak_ck.DF$x),]$x,
                             Minima = SMmin.heur[which.min(SMmin.heur$y),])
          }



        } else{
          #SMmin.OptiBW.try
          SMmin.ls <- SMmin.OptiBW.try
        }
      } else {
        #SMmin.meanbw.try
        SMmin.ls <- SMmin.try
      }

      #end of if nrow()>3
    } else{
      MinimaErr=F; print("Smoothed count DF has low N")
      if(nrow(Gaus_Ker_Smooth_ak_ck.DF)==0) {
        #avoid crashing
        SMmin.ls <- list(PeakA = 0,
                         PeakB = 0,
                         Minima = c(0,0))

      } else{
        SMmin.ls <- list(PeakA = Gaus_Ker_Smooth_ak_ck.DF[which.max(Gaus_Ker_Smooth_ak_ck.DF$y),]$x,
                         PeakB = Gaus_Ker_Smooth_ak_ck.DF[which.min(Gaus_Ker_Smooth_ak_ck.DF$y),]$x,
                         Minima = mean(c(Gaus_Ker_Smooth_ak_ck.DF$x, 0)))
      }
      names(SMmin.ls$PeakA)  = c("x") #c("x", "y")
      names(SMmin.ls$PeakB)  = c("x") #c("x", "y")
      names(SMmin.ls$Minima)  = c("x", "y")
      SMmin.ls$Minima <- as.data.frame(t(SMmin.ls$Minima))

    }

    names(SMmin.ls$PeakA)  = c("x") #c("x", "y")
    names(SMmin.ls$PeakB)  = c("x") #c("x", "y")
    names(SMmin.ls$Minima)  = c("x", "y")



    ak_KeyOptima <- (SMmin.ls)
    #names(ak_KeyOptima) <- c("PeakA",  "PeakB",  "LocMin"/"Minima")

    if(!MinimaErr){

      KeyOptima[[m]] <- ak_KeyOptima



      if(save2file && (nrow(Gaus_Ker_Smooth_ak_ck.DF)>=3)  ){
        fileID = paste(BaseFigDIR, "/alg6_NormVecUpdate_",m, ".png",sep="" )
        png(file = fileID, bg = "transparent", width = 1500, height = 1500, units = "px", res=130)

        #plot(cbind(a_k,c_k),type="n",main=paste("minimizing a_k = ", round(ak_KeyOptima["LocMin"],3), "\n Smoothed Density with Gaus. kernel bw =", round(alg_6_h_bandwidth[m], 4)),xlab="x=a_k",ylab="y=count",
        #     ylim=range(-5:(round(max(c_k))+round(max(c_k)*.2))), xlim=range(-5:5))

        plot(cbind(x_a_k=a_k, y_c_k=c_k), type="l", main=paste("minimizing a_k = ", round(ak_KeyOptima$Minima$x,3), "\n Smoothed Density with Gaus. kernel bw =", round(alg_6_h_bandwidth[m], 4)),xlab="x=a_k",ylab="y=count");
        abline(v=0, lwd=1, lty=3)
        points(Gaus_Ker_Smooth_ak_ck.DF, pch=19, col="red")
        points(Gaus_Ker_Smooth_ak_ck.DF.OptiBW, pch=20, col="plum")
        #abline(h=CountThresh, lwd=1, lty=3)
        #plot(Gaus_Ker_Smooth_ak_ck.DF,main=paste("minimizing a_k = ", round(ak_KeyOptima$Minima$x,3), "\n Smoothed Density with Gaus. kernel bw =", round(alg_6_h_bandwidth[m], 4)),xlab="x=a_k",ylab="y=count")

        #points(cbind(s_j,c_j),col="skyblue",pch=20,cex=0.9)
        #lines(Gaus_Ker_Smooth_ak_ck.DF[,1:2],col="skyblue",lwd=2)
        abline(v=0, lwd=.7, lty=2, col="grey");  abline(h=0, lwd=.7, lty=2, col="grey");
        abline(v=ak_KeyOptima$Minima$x, lty=1, lwd=1, col=col_vector[5])
        abline(v=ak_KeyOptima["PeakA"], lty=1, lwd=1, col=col_vector[1])
        abline(v=ak_KeyOptima["PeakB"], lty=1, lwd=1, col=col_vector[2])

        legend("bottomleft", legend = c("Axis", "Counts" , "Alg6.NormUpdate", "PeakA", "PeakB"),
               col = c("grey", "skyblue", col_vector[5], col_vector[1], col_vector[2]),
               ncol = 1, cex = .7, lwd = c(.7,2, 2,2,2), lty = c(2,1, 1,1,1))

        dev.off()

      }

      w_k_New <- -as.numeric(w_t_org + ak_KeyOptima$Minima$x*v_t_org)
      #w_k_New <- as.numeric(alg2_result$U_robust[1:(length(alg2_result$U_robust)-1)] + ak_KeyOptima["LocMin"]*v_t_org)

      w_b_k_New <- as.numeric(c(w_k_New, b_t_alg4))

      w_k_euc_mag_vec[m] <- norm(w_k_New, type="2") #euclidean norm

      alg6_slope[m] <-  ak_KeyOptima$Minima$x

      alg6_w_new[m, ] <- w_k_New

    }

  }

  remove (TASK);
  colnames(alg6_w_new) <- names(w_t_org)


  return(list(alg6_slope = as.numeric(alg6_slope), alg6_w_new = alg6_w_new, w_k_euc_mag_vec = w_k_euc_mag_vec, KeyOptima=KeyOptima))


}

