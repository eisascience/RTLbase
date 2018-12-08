




alg1_baselineClass <- function(TrainXls,
                               TrainYls,
                               TestXls,
                               TestYls,
                               K_forCrossV,
                               svmGamma,
                               svmCost,
                               prnt2scr,
                               X_cols2Keep,
                               transX=F, sampleRed=F, doParalellSVM=F, datatyp="FC"){


  #initializations
  errorFlag = F


  #Source task data Dm for m = 1, ...., M datasets

  M = length(TrainXls); if(prnt2scr) print(paste("# of datasets = ", M, sep=""))

  if (length(TrainYls) != M) {
    errorFlag = T
    print("training X and Y lengths dont match")
  }



  #Assuming that all df will have the same # of cols;
  if(class(TrainXls[[1]]) == "numeric"){
    nDims = 1
    if(prnt2scr) print("datasets have 1 dim/feat each")

  } else {
    nDims = ncol(as.data.frame((TrainXls[[1]]))[,X_cols2Keep])
    if(prnt2scr) paste("datasets have ", nDims, " dims/feats each", sep="")
  }


  #this will hold the hyperplane parameters for each SVM
  baselineSVM <- matrix(0, nrow=M, ncol=nDims+1)

  #empty vectors for storage of results
  kfoldSVM_ls    <- vector(mode="list")
  results.all    <- vector(mode="list")

  #-------------------end of inputs and initializations


  if(prnt2scr) print(paste("Starting training with a SVM classifier || cost: ",
                           svmCost, ", gamma: ", svmGamma,
                           ", Kernel: linear" , ", cross: ", K_forCrossV, sep=""))

  for (m in 1:M) {
    #m=1

    if(prnt2scr) print(m)

    if(is.na(X_cols2Keep[1])) Dm.train  <- as.data.frame(TrainXls[[m]])
    if(!is.na(X_cols2Keep[1])) Dm.train  <- as.data.frame(TrainXls[[m]])



    if(transX){
      Dm.train <- as.data.frame(RTL::AllDataManipulations(Dm.train,
                                                          Scale=MaxPerCh,
                                                          Center=ScalePerCh,
                                                          X_cols2Keep = X_cols2Keep,
                                                          globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                          globalRange1010Scale = F,
                                                          globalScaleByVal = F))[,X_cols2Keep]
    } else{
      if(!is.na(X_cols2Keep[1])) {
        Dm.train <- Dm.train[,X_cols2Keep]
      } else {
        Dm.train <- as.data.frame(Dm.train)
      }
    }


    if(!is.factor(TrainYls[[m]])) {
      TrainYls[[m]] <- factor(TrainYls[[m]])
    }

    if("TRUE" %in% levels(TrainYls[[m]])) TrainYls[[m]] <- factor(ifelse(TrainYls[[m]]==T,1,-1))

    if("Pos" %in% levels(TrainYls[[m]])) POS = c("Neg", "Pos")
    if("Neg" %in% levels(TrainYls[[m]])) POS = c("Neg", "Pos")

    if("-1" %in% levels(TrainYls[[m]]))  POS = c(-1, 1)
    if("1" %in% levels(TrainYls[[m]]))  POS = c(-1, 1)


    if(ncol(Dm.train)==1){
      Dm.train <- as.data.frame(cbind(Dm.train, factor(TrainYls[[m]], levels = POS)))
      colnames(Dm.train) <- c("X", "labelSVM")
    } else {
      Dm.train$labelSVM        <- factor(TrainYls[[m]], levels = POS)
    }




    if(nDims == 1) colnames(Dm.train) <- c(paste(rep("X", nDims), 1:nDims, sep=""), "labelSVM")

    if(!(sampleRed==F)){
      if(nrow(Dm.train)>sampleRed){
        Dm.train <- Dm.train[sample(1:nrow(Dm.train), sampleRed),]
      }
    }

    #print(nrow(Dm.train))


    if(doParalellSVM==T){
      Dm.train <<- Dm.train
      #print("intiating paralell SVM")
      #for now running in to memory issues. for future dev
      Dm.train_model <- parallelSVM(labelSVM ~. , data=Dm.train,
                                    cost = svmCost , gamma = svmGamma,
                                    type="C-classification",
                                    kernel = 'linear',
                                    cross = K_forCrossV,
                                    scale = F, numberCores = 4)
    } else {
      #print("intiating non-paralell SVM")

      Dm.train_model <- svm(labelSVM ~. , data=Dm.train,
                            cost=svmCost, gamma=svmGamma,
                            type="C-classification",
                            kernel = 'linear',
                            cross = K_forCrossV,
                            scale = F)
    }


    #print("predicting on train set")

    Dm.train.pred  <- predict(Dm.train_model, Dm.train)

    if("Pos" %in% levels(Dm.train.pred)) POS = "Pos"
    if("1" %in% levels(Dm.train.pred)) POS = "1"


    results.all[[m]] <- list(train=conf.mat.stats(conf.mat.pred=Dm.train.pred, conf.mat.truth=Dm.train$labelSVM, POS),
                             DmTrainPred = Dm.train.pred, train_model=Dm.train_model)



    if(doParalellSVM){
      #Wm is the hyperplane coeffs
      Wm <- drop(t(Dm.train_model[[1]]$coefs) %*% Dm.train_model[[1]]$SV)
      #rho is the negative intercept
      bm <- drop(Dm.train_model[[1]]$rho)
    } else {
      #Wm is the hyperplane coeffs
      Wm <- drop(t(Dm.train_model$coefs) %*% Dm.train_model$SV)
      #rho is the negative intercept
      bm <- drop(Dm.train_model$rho)

    }

    baselineSVM[m,] <- as.vector(cbind(t(Wm), bm))

    remove(Dm.train)

  }



  colnames(baselineSVM) <- c(paste(rep("Wm", nDims), 1:nDims, sep=""), "b.int")


  rownames(baselineSVM) <- names(TrainXls)
  names(results.all) <- names(TrainXls)

  return(list(baselineSVM = baselineSVM,
              results.all = results.all,
              M = M,
              nDims = nDims))

}










