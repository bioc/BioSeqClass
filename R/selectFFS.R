selectFFS <- function(data,accCutoff,stop.n,
                   classifyMethod="knn",cv=10){
 if( missing(accCutoff) & missing(stop.n) ){
   stop("Parameter accCutoff and stop.n can not be missing together.")
 }
 columns_now = 1:(ncol(data)-1)
 name = colnames(data)[columns_now]
 
 tmp <- sapply(columns_now,function(x){ 
          classify(data[,c(x,ncol(data))],classifyMethod,cv)[["totalPerformance"]]
        }) 
 tmp1 <- tmp["acc",]
 features_now = columns_now[which(tmp1==max(tmp1))[1]]
 columns_now = setdiff(columns_now, features_now)
 acc_now = max(tmp1)
 result <- list()
 result[[1]] <- list(features=features_now, featureName=name[features_now], 
 performance=tmp[,which(tmp1==max(tmp1))[1]])
 
 for( i in 1:(ncol(data)-2) ){
   print(paste("Time",i,sep=": "))
   if( ! missing(stop.n) ){
     if(i==stop.n){
       break
     }
   }
   tmp <- sapply(columns_now,function(x){ 
     classify(data[,c(features_now,x,ncol(data))], classifyMethod,cv)[["totalPerformance"]]
   })
   tmp1 <- tmp["acc",]
   
   if( ! missing(accCutoff) ){
     if( abs(acc_now-max(tmp1))>accCutoff ){       
       features_now = c(features_now, columns_now[which(tmp1==max(tmp1))[1]] ) 
       columns_now = setdiff(columns_now, features_now)
       acc_now = max(tmp1)
       result[[length(result)+1]] <- list(features=features_now, 
         performance=tmp[,which(tmp1==max(tmp1))[1]])
     }else{
       break
     }
   }else{
     if( ! missing(stop.n) ){
       features_now = c(features_now, columns_now[which(tmp1==max(tmp1))[1]] )       
       columns_now = setdiff(columns_now, features_now)
       acc_now = max(tmp1)
       result[[length(result)+1]] <- list(features=features_now, featureName=name[features_now], 
         performance=tmp[,which(tmp1==max(tmp1))[1]])
     }
   }
 }
 result
}
