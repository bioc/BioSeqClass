performance <- function(predictClass,factClass){
  ## predictClass - a vector indication the predict class label. 1 and -1 are used to express class1 and class2.
  ## factClass - a vector indication the fact class label. 1 and -1 are used to express class1 and class2.
  ##
  ## Prc: precision
  ## Sn: sensitivity
  ## Sp: specificity
  ## Acc: accuracy
  ## MCC: Matthews Correlation Coefficient
  ## PC: Performance Coefficient
  ## AUC: area under the ROC curve (receive operating characteristic curve)
  ##
  
  perform = list()
  tmp = predictClass[predictClass==factClass]
  perform$tp = sum(tmp=="+1")
  perform$tn = sum(tmp=="-1")
  perform$fp = sum(factClass=="-1")-perform$tn
  perform$fn = sum(factClass=="+1")-perform$tp 
  perform$prc = perform$tp/(perform$tp+perform$fp)
  perform$sn = perform$tp/(perform$tp+perform$fn)
  perform$sp = perform$tn/(perform$tn+perform$fp)
  perform$acc = (perform$tp+perform$tn)/(perform$tp+perform$tn+perform$fp+perform$fn)
  perform$mcc = (perform$tp*perform$tn-perform$fn*perform$fp)/((perform$tp+perform$fn)*(perform$tn+perform$fp))^0.5/((perform$tp+perform$fp)*(perform$tn+perform$fn))^0.5
  perform$pc = perform$tp/(perform$tp+perform$fn+perform$fp)
  unlist(perform)
}


