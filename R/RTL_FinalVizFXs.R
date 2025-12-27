


#' Visualize RTL hyperplane performance across algorithms
#'
#' @description
#' Compares baseline, bias-updated, and normal-vector–updated hyperplanes on
#' held-out target tasks, generating summary statistics and ggplot visualizations
#' for downstream reporting.
#'
#' @param TrainTestSet.ls List of train/test splits with `X.test` and `Y.test`
#'   elements.
#' @param alg1_result Output list from [alg1_baselineClass()].
#' @param alg2_result Output list from [alg2_rob_meanNCov()].
#' @param alg3_result Vector of bias updates from [alg3_shiftComp()].
#' @param alg4_result Data frame of bias-adjusted intercepts from
#'   [alg4_BiasUpdate()].
#' @param alg6_result List from [alg6_NormalVectorUpdate()] containing updated
#'   weights.
#' @param datatyp Character flag indicating data type (e.g., `"FC"`).
#' @param ADM Logical; apply additional data manipulations.
#'
#' @return A list of ggplot objects, combined statistics, and per-task
#'   predictions for each algorithm stage.
#' @export
FinalViz <- function(TrainTestSet.ls,
                     alg1_result,
                     alg2_result,
                     alg3_result,
                     alg4_result,
                     alg6_result,
                     datatyp,
                     ADM){

  #TrainTestSet.ls = TestXY.ls

  ChangeAcross <- lapply(1:length(TrainTestSet.ls), function(x.ds){
    #x.ds = 1
    x.ds <- as.numeric(x.ds[1])
    print(x.ds)

    X_test <- as.data.frame(TrainTestSet.ls[[x.ds]]$X.test)
    #X_train <- as.data.frame(TrainTestSet.ls[[x.ds]]$X.train)


    Y_test <- TrainTestSet.ls[[x.ds]]$Y.test
    if(!is.factor(Y_test)) Y_test <- factor(Y_test)
    if("TRUE" %in% levels(Y_test)) Y_test <- factor(ifelse(Y_test==T,1,-1))

    if(is.factor(Y_test)) if("Pos" %in% levels(Y_test)) Y_test <- factor(ifelse(Y_test=="Pos", 1, -1))
    #Y_train <- TrainTestSet.ls[[x.ds]]$Y.train



    # violinMyDF(X_test[which(Y_test==-1),])
    # violinMyDF(X_test[which(Y_test==1),])
    if(ADM) X_test <- as.data.frame(AllDataManipulations(X_test,
                                                 Scale=ScalePerCh,
                                                 Center=MaxPerCh,
                                                 X_cols2Keep = "",
                                                 globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                 globalRange1010Scale = F,
                                                 globalScaleByVal = F))[,ImpFeats]


    if(ncol(X_test)==1){
      X_test <- as.data.frame(cbind(X_test, rep(-1, length(X_test))))
      colnames(X_test) <- c("X", "interc_multp")
    } else {
      X_test$interc_multp <- rep(-1, nrow(X_test))
    }



    #alg2_result$U_robust_norm
    #robust mean as baseline
    w_vec.alg2 <- alg2_result$U_robust[1:(length(alg2_result$U_robust)-1)]
    bint.alg2  <- alg2_result$U_robust[length(alg2_result$U_robust)]

    hyp.alg2 <- -c(w_vec.alg2, b.int=bint.alg2)

    #hyp.alg2 <- c(-w_vec.alg2, b.int=bint.alg2)


    y_hat_alg2 <- factor(sign(as.matrix(X_test) %*% hyp.alg2), levels=c(-1,1))


    if(class(alg4_result$b_alg_norm[x.ds])=="list") Alg4_B <- unlist(alg4_result$b_alg_norm[x.ds])
    if(class(alg4_result$b_alg_norm[x.ds])=="numeric") Alg4_B <- alg4_res$b_alg_norm[x.ds]


    #hyp.alg4 <- -c(w_vec.alg2, unlist(alg4_result$b_alg_norm[x.ds])/alg2_res$w_euc_mag)

    hyp.alg4 <- c(-w_vec.alg2, Alg4_B)

    # plot(density(as.matrix(X_test) %*% hyp.alg4))



    y_hat_alg4 <- factor(sign(as.matrix(X_test) %*% hyp.alg4), levels=c(-1,1))
    # table(pred=y_hat_alg4, truth=Y_test)

    #... all about scaling



    #hyp.alg6 <- as.numeric(c(alg6_result$alg6_w_new[x.ds,], alg4_result$b_alg_norm[x.ds]))
    hyp.alg6 <- c(alg6_result$alg6_w_new[x.ds,], Alg4_B)
    y_hat_alg6 <- factor(sign(as.matrix(X_test) %*% hyp.alg6), levels=c(-1,1))
    #table(pred=y_hat_alg6, truth=Y_test)


    # plot(as.data.frame(cbind((as.matrix(X_test) %*% hyp.alg2),
    #                          (as.matrix(X_test) %*% hyp.alg4),
    #                          as.matrix(X_test) %*% hyp.alg6)), pch=19, col=Y_test)


    #confusionMatrix(table(pred = y_hat_alg6, truth=Y_test), positive = "1")
    #heatmap(cbind(hyp.alg2, hyp.alg4, hyp.alg6))

    if(class(Y_test) == "numeric") Y_test <- factor(Y_test, levels = c(-1, 1))

    if("1" %in% levels(Y_test)) POS = "1"
    if("Pos" %in% levels(Y_test)) POS = "Pos"


    levels(y_hat_alg2) <- levels(Y_test)
    levels(y_hat_alg4) <- levels(Y_test)
    levels(y_hat_alg6) <- levels(Y_test)



    ConfMats4Comp <- list(
      confusionMatrix(table(data = y_hat_alg2, reference=Y_test), positive = POS),
      confusionMatrix(table(data = y_hat_alg4, reference=Y_test), positive = POS),
      confusionMatrix(table(data = y_hat_alg6, reference=Y_test), positive = POS))

    compStats <- vector(mode="list")
    compStats[["accuracy"]] <- rbindlist(lapply(ConfMats4Comp, function(x){
      as.data.frame(t(x$overall[c(1,3,4)]))
    }))
    compStats[["rest"]] <- rbindlist(lapply(ConfMats4Comp, function(x){
      as.data.frame(t(x$byClass))
    }))

    compStats <- melt(compStats)
    compStats[is.na(compStats)] <- 0

    compStats$L2 <- factor(rep(c("2.baseline", "1.alg4", "0.alg6"), nrow(compStats)/3), ordered = T)



    gg <- ggplot(compStats,aes(x=factor(variable), y=(value*100), colour=factor(L2))) +
      geom_point(aes(size=factor(L2))) +
      theme_bw() +
      theme(plot.title = element_text( color="#666666", face="bold", size=25, hjust=0.5)) +
      theme(axis.title = element_text( color="#666666", face="bold", size=20)) + ylim(-5,105) +
      labs(title = paste("Classification by Hyperplanes\n", sep=""), y = "value(%)", x = "statistic")+
      theme(axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_colour_manual(values = col_vector)

    gg
    #plot(gg)


    fileID = paste(BaseFigDIR, "/alg6_ClassificationCompB_", names(TrainTestSet.ls)[x.ds], ".png",sep="" )

    ggsave(fileID, width = 30, height = 30, units = "cm", dpi=100)


    #pp <- plotly()
    #pp$ggplotly(gg, session="knitr")



    return(list(compStats=compStats, ggfigA=gg,
                y_hat_alg6=y_hat_alg6, hyp.alg6=hyp.alg6,
                y_hat_alg4=y_hat_alg4, hyp.alg4=hyp.alg4,
                y_hat_alg2=y_hat_alg2, hyp.alg2=hyp.alg2, y_test = Y_test))


  })

  comStats <- rbindlist(lapply(1:length(ChangeAcross), function(x){

    tempDF <- as.data.frame(ChangeAcross[[x]]$compStats)
    tempDF$id <- rep(x, nrow(tempDF))
    tempDF

  }))

  comStats.sub <- comStats[which(comStats$variable %in% c("Accuracy", "Specificity", "Precision", "Recall", "F1", "Prevalence", "Detection Rate")),]

  comStats.sub$variable <- as.character(comStats.sub$variable)
  comStats.sub[which(comStats.sub$variable == "Recall"),"variable"] <- "Recall/Sensitivity"
  comStats.sub$variable <- factor(comStats.sub$variable)


  gg3 <- ggplot(comStats.sub, aes(x=factor(id), y=(value*100), colour=factor(L2))) +
    geom_point(aes(size=factor(L2))) +
    facet_wrap(~factor(variable)) +
    theme_bw() + theme(plot.title = element_text( color="#666666", face="bold", size=25, hjust=0.5)) +
    theme(axis.title = element_text( color="#666666", face="bold", size=20)) + ylim(-5,105) +
    labs(title = paste("Classification by Hyperplanes\n", sep=""), y = "value(%)", x = "statistic") +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_colour_manual(values = col_vector)

  plot(gg3)

  fileID = paste(BaseFigDIR, "/alg6_ClassificationCompALL", ".png",sep="" )

  ggsave(fileID, width = 60, height = 60, units = "cm", dpi=120)

  #pp3 <- plotly()
  #pp3$ggplotly(gg3, session="knitr")



  return(list(gg.figs=list(A=gg3),
              comStats = comStats,
              comStats.sub = comStats.sub,
              ChangeAcross = ChangeAcross))

}




