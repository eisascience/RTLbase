


DemoTestTrainList2D <- function(N=0){
  if(N!=0){
    tempLS <- lapply(1:N, function(n){
      XY.DF <- as.data.frame(cbind(X1=c(rnorm(500, 0), rnorm(500, 10)),
                                   X2=c(rnorm(800, 0), rnorm(200, 5))))

      XY.DF$Y <- c(rep(-1, 800), rep(1, 200))

      #plot(XY.DF[,1:2], pch=20, col=factor(XY.DF$Y))

      TrainID <- sample(1:1000, 800, replace = F)
      TestID <- setdiff(1:1000, TrainID)


      Train.XY.DF <- XY.DF[TrainID,]
      Test.XY.DF  <- XY.DF[TestID,]
      list(TrainX = Train.XY.DF[,1:2],
           TrainY = factor(Train.XY.DF[,3]),
           TestX  = Test.XY.DF[,1:2],
           TestY  = factor(Test.XY.DF[,3]))
    })
  } else {
    print("N > 0 please")
  }
  names(tempLS) <- paste("set", 1:N, sep="")
  return(tempLS)
}





