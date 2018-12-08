
ColorTheme <- function(){
  #scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
  scaleyellowred <- colorRampPalette(c("dodgerblue", "lightyellow", "red"), space = "rgb")(30)

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- col_vector[-4]
  return(list(col_vector=col_vector, scaleyellowred=scaleyellowred))
}

conf.mat.stats <- function(conf.mat.pred, conf.mat.truth, POS="1"){

  #conf.mat.pred   =  factor(sign(y_hat[,x_m]), levels=c(-1,1))
  #conf.mat.truth  =  TASK$TrueClass

  temp.confMatdat <- confusionMatrix(data=conf.mat.pred, reference=conf.mat.truth, positive=POS, mode="everything")
  return(list(overall= temp.confMatdat$overall, table=temp.confMatdat$table, byClass=temp.confMatdat$byClass))

}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

IndicatorFX <- function(x){
  if(x==T) return(1)
  if(x==F) return(0)
  if(!(x %in% c(T, F))) {
    print("Not T/F")}
}


is.odd <- function(x){
  x%%2 == 1
}
is.even <- function(x){
  x%%2 == 0
}






DF2TrainTestls <- function(X, Y){
  #X = BCC$PBMC
  #Y = BCC$PBMC.Class
  DF2KfoldBYClass <- function(tSDF, tSYV, Ksplits=10){
    #use factor to split first to make sure both cases are in each set

    # tSDF     = BCC$PBMC
    # tSYV     = BCC$PBMC.Class
    # Ksplits  = 10

    if(!is.factor(tSYV)) tSYV <- factor(tSYV)

    ClassA.DF <- subset(cbind(tSDF, tSYV), tSYV == levels(tSYV)[1])
    ClassB.DF <- subset(cbind(tSDF, tSYV), tSYV == levels(tSYV)[2])
    err = F
    if(nrow(ClassA.DF) < Ksplits) {
      print("cant split, too few A cases")
      err = T
    }
    if(nrow(ClassB.DF) < Ksplits) {
      print("cant split, too few B cases")
      err = T
    }
    if(!err){

      tempSampAIDs <- 1:nrow(ClassA.DF)
      tempSampBIDs <- 1:nrow(ClassB.DF)
      options(warn=-1)
      tempSampAIDs.sp <- split(sample(tempSampAIDs), 1:Ksplits)
      tempSampBIDs.sp <- split(sample(tempSampBIDs), 1:Ksplits)
      options(warn=0)
      Y.split.ls <- lapply(1:length(tempSampBIDs.sp), function(n){
        #n = 1
        tempYV <- c(ClassA.DF[tempSampAIDs.sp[[n]], "tSYV"],
                    ClassB.DF[tempSampBIDs.sp[[n]],"tSYV"])
        factor(tempYV, levels = c(1,2), labels = levels(tSYV))
      })

      X.split.ls <- lapply(1:length(tempSampBIDs.sp), function(n){
        #n = 1

        rbind(ClassA.DF[tempSampAIDs.sp[[n]],
                        setdiff(colnames(ClassA.DF), "tSYV")],
              ClassB.DF[tempSampBIDs.sp[[n]],
                        setdiff(colnames(ClassA.DF), "tSYV")])[,]
      })

      return(lapply(1:length(X.split.ls), function(leng){
        #leng = 1
        list(x=X.split.ls[[leng]],y=Y.split.ls[[leng]], KselectRows=c(tempSampAIDs.sp[[leng]], tempSampBIDs.sp[[leng]]))

      }))
    }

  }

  TestXY.ls <- DF2KfoldBYClass(X, Y)
  names(TestXY.ls) <- paste("TestSubSec", 1:length(TestXY.ls), sep="")

  Nminus1Train <- function(testXYlist){
    #testXYlist = TestXY.ls

    lapply(names(testXYlist), function(Nm){
      #Nm=names(testXYlist)[1]
      list(x = as.data.frame(rbindlist(lapply(names(testXYlist)[which(!(names(testXYlist) %in% Nm))], function(xNm){
        #xNm=names(testXYlist)[2]
        as.data.frame(testXYlist[[xNm]]$x)
      }))),
      y = unlist(lapply(names(testXYlist)[which(!(names(testXYlist) %in% Nm))], function(xNm){
        #xNm=names(testXYlist)[2]
        (testXYlist[[xNm]]$y)
      })))

    })

  }

  TrainXY.ls <- Nminus1Train(TestXY.ls)
  names(TrainXY.ls) <- paste("TrainSubSec", 1:length(TrainXY.ls), sep="")


  return(list(TrainXY.ls = TrainXY.ls, TestXY.ls = TestXY.ls))
}
