
#Mahyari, Eisa: Implementation of Lee et al.'s TL-FC into R
#March, 2017

#computing the kernel density estimate (KDE) of the projections onto w_0
#align the target data to each source data using maximum cross-correlation
# in a leave-one-out setting

##### Shift Compensation (Algorithm 3) ##########
###
### Input: hyperplane (w, b),
###    source task data {D_m} for m=1 to M, and
###    target task data Tau
###
###    z_t_i <- <w,x_m_i> + b, for all i
###
###    for m = 1 to M do
###        z_m_i <- <w,x_m_i> + b, for all i
###        e_m <- argmax_z KDE(z,z_t_i) *max_cross-cor* KDE(z, z_m_i)
###    end for
###
###    b <- b - median(e_m)
###
### Output: b
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

ccfmaxV3 <- function(a, b, e=0, maxLag, useAbsCor = T)
{
  # a = c(rnorm(5000,1,1), rnorm(5000,6,1))
  # b = c(rnorm(5000,2,1), rnorm(5000,8,1))
  # plot(density(a))
  # lines(density(b))
  # maxLag = 10


  #if(useAbsCor) print("AbsCor")
  #d <- ccf(a, b, plot = FALSE, lag.max = length(a)/2)
  d <- ccf(a, b, plot = FALSE, lag.max = maxLag)
  cor = round(d$acf[,,1], 3)
  if(useAbsCor) cor = abs(cor)
  lag = d$lag[,,1]
  res = data.frame(cor, lag)
  res_max = res[which.max(res$cor),]



  return(res_max)
}

#' Apply shift compensation using maximum cross-correlation
#'
#' @description
#' Aligns target task projections to source datasets by estimating robust lag
#' shifts via maximum cross-correlation between kernel density estimates. The
#' updated intercepts are returned for downstream bias correction.
#'
#' @param source_list List of source feature matrices used to train the
#'   classifiers.
#' @param task_list List of target task feature matrices to be aligned.
#' @param alg1_result Output list from [alg1_baselineClass()].
#' @param alg2_result Output list from [alg2_rob_meanNCov()].
#' @param print2screen Logical; print progress updates.
#' @param save2file Logical; save diagnostic plots to disk.
#' @param maximumLag Maximum lag to consider in cross-correlation.
#' @param ImpFeats Optional vector of feature indices to retain.
#' @param ADM Logical; apply additional data manipulations.
#' @param datatyp Character flag indicating data type (e.g., `"FC"`).
#' @param useAbsCor Logical; use absolute correlations when selecting lags.
#' @param medianMediansBL Logical; use the median of medians strategy for shift
#'   estimation.
#' @param CoreClassifier Core classifier identifier (e.g., `"LinSVM"`).
#'
#' @return A numeric vector of updated intercepts per task.
#' @export
alg3_shiftComp <- function(source_list, task_list, alg1_result, alg2_result,
                           print2screen, save2file, maximumLag, ImpFeats,
                           ADM=F, datatyp="FC", useAbsCor = T, medianMediansBL = F,
                           CoreClassifier){

  # task_list      = TestXls.t;
  # source_list    =  TrainXls.t;
  # alg2_result    =  alg2_res;
  # alg1_result    = alg1_res
  # print2screen   =  T;
  # save2file      =  F;
  # maximumLag     =  0;
  # ImpFeats = ImpFeats

  n_testSets = length(task_list)
  M = length(source_list)




  if(print2screen) print("starting Alg 3, shift compensation")

  #alg1task_hat <- vector(mode="list")
  #alg1source_hat <- vector(mode="list")

  w_vec.alg2 <- alg2_result$U_robust[1:(length(alg2_result$U_robust)-1)]
  bint.alg2  <- alg2_result$U_robust[length(alg2_result$U_robust)]



  ###########################CLASSIFIER Selection
  if(CoreClassifier=="LinSVM") hyp.alg2  <- -c(w_vec.alg2, b.int=bint.alg2)
  ###########################CLASSIFIER Selection



  task_hat <- vector(mode="list")
  source_hat <- vector(mode="list")

  for (m in 1:M) {

    # if(!exists("m")) m = 1

    if(print2screen) print(paste("Inference of class by the datasets's baseline hyperplane #", m,sep=""))

    if(is.na(ImpFeats)[1]){
      SOURCE <- as.data.frame(source_list[[m]])[,]
    } else {
      SOURCE <- as.data.frame(source_list[[m]])
    }

    #violinMyDF(SOURCE, title=paste("SOURCE data #", m, sep=""), y.title="value", x.title="feats")




    #boxplot(SOURCE)
    #change ADM to TransX
    if(ADM) SOURCE <- as.data.frame(AllDataManipulations(SOURCE,
                                                         Scale=ScalePerCh,
                                                         Center=MaxPerCh,
                                                         globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                         globalRange1010Scale = F,
                                                         globalScaleByVal = F))[,ImpFeats]

    #helps if 1D data vs nD
    if(ncol(SOURCE)==1){
      ###########################CLASSIFIER Selection
      SOURCE <- as.data.frame(cbind(SOURCE, rep(ifelse(CoreClassifier=="LinSVM", -1, 1), length(SOURCE))))
      colnames(SOURCE) <- c("X", "interc_multp")
    } else {
      ###########################CLASSIFIER Selection
      SOURCE$interc_multp <- rep(ifelse(CoreClassifier=="LinSVM", -1, 1), nrow(SOURCE))
    }
    source_hat[[m]] <- as.data.frame(as.matrix(SOURCE) %*% hyp.alg2)




    #-1 beacause the $rho returned from the svm() fx is the neg int

    # #z_m_i
    # w_vec.alg2  <- alg2_result$U_robust[1:(length(alg2_result$U_robust)-1)]
    # bint.alg2   <- alg2_result$U_robust[length(alg2_result$U_robust)]
    # hyp.alg2 <- -c(w_vec.alg2, b.int=bint.alg2)

    #hyp.alg2 <- c(w_vec.alg2, b.int=bint.alg2)



    # table(pred=factor(sign(source_hat[[m]])[,1], levels=c(-1,1)),tru=factor(Kfold.TrainingSets.ls[[1]]$y[[m]], levels=c(-1,1)))

    #plot(density(source_hat[[m]][,1]))


  }
  names(source_hat) <- names(source_list)

  for (j in 1:n_testSets) {

    #if(!exists("j")) j = 1

    if(print2screen) print(paste("Mapping to Y_hat for Task: ", j, sep=""))

    if(is.na(ImpFeats)[1]){
      TASK <- as.data.frame(task_list[[j]])[,]
    } else {
      TASK <- as.data.frame(task_list[[j]])
    }

    if(ADM) TASK <- as.data.frame(AllDataManipulations(TASK,
                                                       Scale=ScalePerCh,
                                                       Center=MaxPerCh,
                                                       globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                       globalRange1010Scale = F,
                                                       globalScaleByVal = F))[,ImpFeats]
    #boxplot(TASK)

    # ggplot(melt(TASK), aes(x = variable, y = value)) + geom_violin() +
    #   theme_bw() +
    #   theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=25, hjust=0.5)) +
    #   theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=20)) +
    #   labs(title = paste("TASK data #", j, sep=""), y = "value", x = "channels")+
    #   theme(axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
    #   scale_colour_manual(values = col_vector)



    if(ncol(TASK)==1){
      ###########################CLASSIFIER Selection
      TASK <- as.data.frame(cbind(TASK, rep(ifelse(CoreClassifier=="LinSVM", -1, 1), length(TASK))))
      colnames(TASK) <- c("X", "interc_multp")
    } else {
      ###########################CLASSIFIER Selection
      TASK$interc_multp <- rep(ifelse(CoreClassifier=="LinSVM", -1, 1), nrow(TASK))
    }


    task_hat[[j]] <- as.data.frame(as.matrix(TASK) %*% hyp.alg2)
    #z_t_i


  }
  names(task_hat) <- names(task_list)

  if(print2screen) print("step3 reached, task and source are mapped, starting Max-Cross Corr")

  e_med_mat <- vector(mode="list")
  cor_e_med_mat <- vector(mode="list")
  MS_med_mat <- vector(mode="list")



  for (t in 1:n_testSets){
    #default behavior
    #t=1
    if(print2screen) print(paste("starting cross-cor for Task: ", t, sep=""))
    if (maximumLag == 0){
      if(length(source_hat[[1]][,1])<50){
        maximumLag   =  length(source_hat[[1]][,1])
      } else{
        maximumLag   =  length(source_hat[[1]][,1])/10
      }

    }

    e <- vector()
    cor_e <- vector()
    MS <- vector()

    #M source, T tasks
    for (i in 1:M) {
      #i=1; t=1

      #need to oversmooth, else cross correlations can be spourious
      DensA <- density(task_hat[[t]][,1], n=length(task_hat[[t]][,1]))
      bwsig = DensA$bw
      if(print2screen) print(paste("with source: ", i, sep=""))
      #plot(DensA)

      CorThresh = 0.9
      #set rep4Rpba = 1 to match Lee et al. version
      rep4Roba = 3 #50

      #plot(density(source_hat[[i]][,1]))

      RobustCrossCorLagALL <- lapply(1:rep4Roba, function(id){
        #id = 1
        set.seed(abs(round(rnorm(id)*30000))[id])

        DensB1 <- density(source_hat[[i]][,1][sample(1:length(source_hat[[i]][,1]), min(c(length(source_hat[[i]][,1]), length(task_hat[[t]][,1]))))], bw=bwsig, n=length(task_hat[[t]][,1]))

        e_resamp <- ccfmaxV3(DensA$y,
                             DensB1$y,
                             maxLag = maximumLag, useAbsCor = useAbsCor);
        if(e_resamp$cor >= CorThresh) return(e_resamp)

      })

      RobustCrossCorLag <- as.numeric(as.character(lapply(RobustCrossCorLagALL, function(x){
        ifelse(is.null(x$lag), NA, x$lag)
      })))

      RobustCrossCorLag <- round(mean(RobustCrossCorLag[!is.na(RobustCrossCorLag)]))

      while(is.na(RobustCrossCorLag)){
        if(print2screen) print(paste("robust cross-cor failed w. threshold: ", CorThresh, sep=""))
        CorThresh <- CorThresh - 0.01
        rep4Robb = 10 #50
        RobustCrossCorLagALL <- lapply(1:rep4Robb, function(id){
          set.seed(abs(round(rnorm(1)*30000)))
          DensB1 <- density(source_hat[[i]][,1][sample(1:length(source_hat[[i]][,1]), min(c(length(source_hat[[i]][,1]), length(task_hat[[t]][,1]))))], bw=bwsig, n=length(task_hat[[t]][,1]))
          #plot(DensA, col="dodgerblue")
          #lines(DensB1)
          #InterpolatedDens <- approx(DensA$y, DensB1$y)
          e_resamp <- ccfmaxV3(DensA$y,
                               DensB1$y,
                               maxLag = maximumLag, useAbsCor = useAbsCor);
          if(e_resamp$cor > CorThresh) return(e_resamp)

        });
        if(print2screen) print(CorThresh)

        RobustCrossCorLag.t <- as.numeric(as.character(lapply(RobustCrossCorLagALL, function(x){
          ifelse(is.null(x$lag), NA,x$lag)
        })))

        RobustCrossCorLag <- Mode(RobustCrossCorLag.t[!is.na(RobustCrossCorLag.t)])

      }
      if(print2screen) print(paste("robust cross-cor found @ threshold: ", CorThresh, sep=""))



      RobustCrossCor <- as.numeric(as.character(lapply(RobustCrossCorLagALL, function(x){
        x$cor
      })))
      RobustCrossCor <- Mode(RobustCrossCor[!is.na(RobustCrossCor)])

      DensB1 <- density(source_hat[[i]][,1], n=length(source_hat[[i]][,1]))

      # plot(DensA, col="dodgerblue")
      # lines(DensB1)


      if(RobustCrossCorLag <  0) deltaLag = (DensB1$x[abs(RobustCrossCorLag)] - DensA$x[1])
      if(RobustCrossCorLag >  0) deltaLag = -(DensB1$x[abs(RobustCrossCorLag)] - DensA$x[1])
      if(RobustCrossCorLag == 0) deltaLag = 0

      # plot(x=DensA$x, y=DensA$y)
      # lines(x=DensB1$x, y=DensB1$y)


      #deltaLag <- - e_resamp[,2]
      MedianShift <- median(task_hat[[t]][,1]) - median(source_hat[[i]][,1])


      if (save2file) {

        fileID = paste(BaseFigDIR, "/alg3Shifts_", paste("t_",t,"_s_",i,sep=""), ".png",sep="" )

        png(fileID, width = 1024, height = 768, units = "px")

        par(mfrow=c(2,3))
        boxplot(list(task=task_hat[[t]][,1], source=source_hat[[i]][,1]),
                main=paste("task: ", t, " | source: ", i, "\n", names(task_hat)[t], " | ",
                           names(source_hat)[i], sep=""),
                ylab = "Prediction (y_hat)")

        ccfig <- ccf(DensA$y,
                     DensB1$y,
                     plot = F, lag.max = maximumLag)
        corTemp = round(ccfig$acf[,,1], 3)
        if(useAbsCor) corTemp = abs(corTemp)
        lagTemp = ccfig$lag[,,1]

        plot(data.frame(lagTemp, corTemp), main=paste("Maximum cross-correlation\n@ robust-lag:", RobustCrossCorLag, " | t_",t,"_s_",i,sep=""), xlab="Lag", ylab="abs(cor)")
        abline(v=RobustCrossCorLag, lwd=2, lty=2, col="red")

        plot(DensA, type='l',lwd=2, col="dodgerblue", main="pre-shift")
        #lines(density(task_hat[[t]][,1], bw=.1), type='l', col="red")
        lines(DensB1, lwd=2, col='firebrick')


        D2<- as.data.frame(cbind(x=DensA$x - deltaLag,
                                 y=DensA$y))

        plot(D2, lwd=2, type='l', col="dodgerblue",
             main=paste("post-shift w/ delta : ",round(deltaLag,3),
                        "\nMax(abs(Cor)) >= thr",
                        " and lag : ", round(RobustCrossCorLag,3),sep=""))
        lines(DensB1, lwd=2, col="firebrick")
        #lines(DensA, type='l',lwd=1, lty=3, col="black")


        D2<- as.data.frame(cbind(x=DensA$x - MedianShift, y=DensA$y))

        plot(D2, type='l', lwd=2, col="dodgerblue",
             main=paste("post-shift w/ delta : ",round(MedianShift,3),
                        "\nMedian-based"))
        lines(DensB1, lwd=2, col="firebrick")


        dev.off()

      }

      e[i] <- deltaLag
      cor_e[i] <- RobustCrossCor
      MS[i] <- MedianShift
      #e[i] <- e_resamp[,2]
    }
    e_med_mat[[t]] <- e
    cor_e_med_mat[[t]] <- cor_e
    MS_med_mat[[t]] <- MS

    maximumLag = 0
  }

  #names(e_med_mat) <- names(task_hat)
  #names(cor_e_med_mat) <- names(task_hat)
  #names(MS_med_mat) <- names(task_hat)

  medianShifts <- as.numeric(lapply(e_med_mat, function(x){
    median(x, na.rm = T)}))
  medianMedianShifts <- as.numeric(lapply(MS_med_mat, median))
  e_med_matDF <- as.data.frame(e_med_mat)
  #rownames(e_med_matDF) <- names(source_hat)
  cor_e_med_matDF <- as.data.frame(cor_e_med_mat)
  #rownames(cor_e_med_matDF) <- names(source_hat)
  MS_med_matDF <- as.data.frame(MS_med_mat)
  #rownames(MS_med_matDF) <- names(source_hat)

  remove(t)
  if (save2file) {


    fileID = paste(BaseFigDIR, "/alg3ShiftsCombo_", ".png",sep="" )

    png(fileID, width = 1024, height = 768, units = "px")
    par(mfrow=c(1,3))
    i=1
    D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=bwsig)$x, y=density(source_hat[[i]][,1], bw=.1)$y))
    plot(D2, col=col_vector[i], type='l', main="Task_hats pre-shift", lwd=2)
    if(length(task_hat)>1){
      for (i in 2:length(medianShifts)){
        D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=bwsig)$x, y=density(task_hat[[i]][,1], bw=.1)$y))
        lines(D2, col=col_vector[i], type='l', main="Task_hats\nmedian shifts", lwd=2)
      }
    }


    i=1
    D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=bwsig)$x-medianShifts[i], y=density(source_hat[[i]][,1], bw=.1)$y))

    plot(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of lag-shifts", lwd=2)
    if(length(task_hat)>1){
      for (i in 2:length(medianShifts)){
        D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=.5)$x-medianShifts[i], y=density(task_hat[[i]][,1], bw=.1)$y))
        lines(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of lag-shifts", lwd=2)
      }
    }


    i=1
    D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=.5)$x-medianMedianShifts[i], y=density(source_hat[[i]][,1], bw=.1)$y))
    plot(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of Median-shifts", lwd=2, xlim=range(-10:10))
    if(length(task_hat)>1){
      for (i in 2:length(medianMedianShifts)){
        D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=.5)$x-medianMedianShifts[i], y=density(task_hat[[i]][,1], bw=.1)$y))
        lines(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of Median-shifts", lwd=2)
      }
    }

    dev.off()

    e_med_matDF$median <- apply(e_med_matDF, 1, median)
    print(xtable(t(cor_e_med_matDF), caption='Bias Shift Correlations in Test (rows) Vs. Training (columns) ',
                 label='tab:cor_shifts', digits=3),
          file=paste(PermLogSaveDirV,
                     "/cor_shifts.txt", sep=""))
    print(xtable(t(e_med_matDF), caption='Bias Shift Values in Test (rows) Vs. Training (columns) ',
                 label='tab:cor_shifts', digits=3),
          file=paste(PermLogSaveDirV,
                     "/val_shifts.txt", sep=""))

    print(xtable(t(MS_med_matDF), caption='Bias Shift Values in Test (rows) Vs. Training (columns) ',
                 label='tab:medMed_shifts', digits=3),
          file=paste(PermLogSaveDirV,
                     "/val_Medshifts.txt", sep=""))




  }

  if(CoreClassifier=="LinSVM") {
    #b_updated <- alg2_result$U_robust_norm[length(alg2_result$U_robust_norm)] + -1*abs(as.numeric(lapply(e_med_mat, median)))
    b_updated <- hyp.alg2[length(alg2_result$U_robust)] - (as.numeric(lapply(e_med_mat, median)))
    names(b_updated) <- paste("Task", 1:length(b_updated), "_MedianOfTasknSourceLags", sep="")
    if(medianMediansBL){
      b_updated <- alg2_result$U_robust[length(alg2_result$U_robust)] - medianMedianShifts
    }
    ##as.numeric(lapply(e_med_mat, median)))
  }


  if(print2screen) print("Shift Compesation complete for all tasks")
  if(print2screen) print("the normed intercepts (b_t) have been updated by median lag in max cross-correlation")

  if(print2screen) print("the shifts for each task:")
  if(print2screen) print(as.numeric(lapply(e_med_mat, median)))
  if(print2screen) print("Final intercepts for all tasks:")
  if(print2screen) print(b_updated)
  return(b_updated)

}





