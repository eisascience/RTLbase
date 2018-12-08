
library(RTLbase)

col_vector <- RTLbase::ColorTheme()$col_vector

#a place to store figs and tables
BaseFigDIR = paste(getwd(), "results", "BCellClassification003", sep="/")
if(dir.exists(BaseFigDIR)==F) dir.create(BaseFigDIR, recursive = T)


#Create some training and testing data
#bimodal, 2 dims

# NoOfTestTrainSamps = 5
# source("./R/simData.R")


#A list of training and testing data
#each set in list has 4 segments
#TestX, TestY, TrainX, TrainY
#Hide for inputing Any other data in such a list
#TestTrain.ls <- DemoTestTrainList2D(NoOfTestTrainSamps)

BCC <- readRDS("./data/BcellClassification.rds")


# TestXY.ls <- DF2KfoldBYClass(BCC$PBMC, BCC$PBMC.Class, 10)
# names(TestXY.ls) <- paste("TestSubSec", 1:length(TestXY.ls), sep="")

TrainTestXY.ls <- DF2TrainTestls(BCC$PBMC, BCC$PBMC.Class)
names(TrainTestXY.ls)


##########
TrainXY.ls <- lapply(TrainTestXY.ls$TrainXY.ls, function(TTls){
  list(X.train=TTls$x, Y.train=TTls$y)
})
TestXY.ls <- lapply(TrainTestXY.ls$TestXY.ls, function(TTls){
  list(X.test=TTls$x, Y.test=TTls$y)
})
remove(TestTrain.ls)

#########end of data simulation

head(TrainXY.ls$TrainSubSec1$X.train)

par(mfrow=c(1,1))
plot(as.data.frame(TrainXY.ls$TrainSubSec1$X.train[,1:5]),
     pch=20, col=TrainXY.ls$TrainSubSec1$Y.train)
plot(as.data.frame(TrainXY.ls$TrainSubSec1$X.train[,6:10]),
     pch=20, col=TrainXY.ls$TrainSubSec1$Y.train)
plot(as.data.frame(TrainXY.ls$TrainSubSec1$X.train[,11:13]),
     pch=20, col=TrainXY.ls$TrainSubSec1$Y.train)

####baseline classifier

alg1_res <- alg1_baselineClass(
  TrainXls = lapply(TrainXY.ls, function(x){x$X.train}),
  TrainYls = lapply(TrainXY.ls, function(x){x$Y.train}),
  TestXls = lapply(TestXY.ls, function(x){x$X.test}),
  TestYls = lapply(TestXY.ls, function(x){x$Y.test}),
  K_forCrossV   =  10,
  svmGamma      =  .02,
  svmCost       =  .5,
  prnt2scr      =   T,
  X_cols2Keep   =  1:13,
  transX=F, sampleRed=F, doParalellSVM=F, datatyp="FC")


alg1_training.DF <- rbind(as.data.frame.list(lapply(alg1_res$results.al, function(sets){
  round(sets$train$overall*100,2)
})), as.data.frame.list(lapply(alg1_res$results.al, function(sets){
  round(sets$train$byClass*100,2)
}))); alg1_training.DF

alg1_res$baselineSVM

####Robust Mean and Cov (Alg2)
alg2_res <- alg2_rob_meanNCov(alg1_res$baselineSVM)
alg2_res



###Max Cross Cor (Alg3)

alg3_res <- alg3_shiftComp(
  task_list = lapply(TestXY.ls, function(x){x$X.test}),
  source_list    =  lapply(TrainXY.ls, function(x){x$X.train}),
  alg2_result    =  alg2_res,
  print2screen   =  F,
  ImpFeats = "",
  save2file      =  F,
  maximumLag = 0,
  CoreClassifier="LinSVM",
  datatyp="FC",
  useAbsCor = T,
  medianMediansBL = F)

alg3_res


###Bias Update (Alg4)

alg4_res <- alg4_BiasUpdate(task_list = lapply(TestXY.ls, function(x){x$X.test}),
                                      alg1_result = alg1_res,
                                      alg2_result = alg2_res,
                                      alg3_result = alg3_res,
                                      goodColumns = "", alg4MinFx = "gd",
                                      Marg = .3,
                                      save2file =F,
                                      ADM=F,
                                      useMedian = T,
                                      ZnormMappingBL=F,
                                      datatyp="FC",
                                      RCSmodeBL = F,
                                      CoreClassifier = "LinSVM")

###Hyperplane rotation update (Alg6)

alg6_res <- alg6_NormalVectorUpdate(task_list = lapply(TestXY.ls, function(x){x$X.test}),
                                    alg1_result = alg1_res,
                                    alg2_result = alg2_res,
                                    alg3_result = alg3_res,
                                    alg4_result = alg4_res,
                                    X_feat_cols = "",
                                    save2file = F,
                                    Marg = .3,
                                    ADM=F,
                                    datatyp="FC",
                                    RCSmodeBL = F,
                                    CoreClassifier = "LinSVM")


####Final Visualization of results
Viz_res <- FinalViz(TrainTestSet.ls = TestXY.ls,
                    alg1_result = alg1_res,
                    alg2_result = alg2_res,
                    alg3_result = alg3_res,
                    alg4_result = alg4_res,
                    alg6_result = alg6_res,
                    datatyp =  "FC",
                    ADM = F)



#Extra results dive

dtb.CSS <- as.data.table(Viz_res$comStats.sub)
dtb.CSS[, Mean:=mean(value, na.rm = T), by=list(variable, L2)]
dtb.CSS[, CI95:=quantile(value, c(.95), na.rm = T) , by=list(variable, L2)]
dtb.CSS[, CI05:=quantile(value, c(.05), na.rm = T) , by=list(variable, L2)]

dtb.CSS[, Meankfolds:=mean(value, na.rm = T), by=list(variable, L2)]
dtb.CSS[, MeankfoldsLCI:=quantile(value, c(0.05), na.rm = T), by=list(variable, L2)]
dtb.CSS[, MeankfoldsUCI:=quantile(value, c(0.95), na.rm = T), by=list(variable, L2)]

#classification results compared to baseline (non-adapated mean hyperplane)
gg6 <- ggplot(dtb.CSS, aes(x=factor(L2), y=(Meankfolds*100), fill=factor(L2))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=MeankfoldsLCI*100, ymax=MeankfoldsUCI*100),
                width=.2, position=position_dodge(.9)) +
  facet_wrap(~factor(variable), ncol=2) +
  theme_bw() +
  theme(plot.title = element_text( color="#666666", face="bold", size=25, hjust=0.5)) +
  theme(axis.title = element_text( color="#666666", face="bold", size=20)) + ylim(-5,105) +
  labs(title = paste("Classification by Hyperplanes\n", sep=""), y = "Mean (%) +/- 95CI", x = "statistic") +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = col_vector)
plot(gg6)






