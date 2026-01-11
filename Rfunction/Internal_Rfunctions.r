
Tuning_setup_group_only = function(tauvec, subtuning, Lambda_limit, Lenlam, medianval, Nq){
   
   lambdavec = seq((min(Lambda_limit)),(max(Lambda_limit)),len=Lenlam)
   if(Nq == 1){
     ratios = 1
   }else{
     ratios = seq(1, 0, length.out=subtuning)
   }
   Start = 1
   for(i1 in 1:length(tauvec)){
     tau = tauvec[i1]
     for(ss in 1:length(lambdavec)){
       total = lambdavec[ss]
       subs1 = total*ratios
       subs2 = total*(1-ratios)
       if(Start == 1){
         TuningsMatrix = cbind(subs1,0,subs2*tau,tau)
         Start = 0
       }else{
         TuningsMatrix = rbind(TuningsMatrix,cbind(subs1,0,subs2*tau,tau))
       }
     }
   }
   test = apply(TuningsMatrix,1,paste0,collapse=":")
   uniq_test = unique(test)
   TuningsMatrix = TuningsMatrix[match(uniq_test,test),]
   return(TuningsMatrix)
}




Tuning_setup_group_func = function(lambdavec_func, lambdavec_func_limit_len, numfunc,  tauvec, subtuning, Lambda_limit, Lenlam, medianval, equal_lambdaf = 1, Nq){

  if(is.null(lambdavec_func)){
    lambdavec_func = seq(0, lambdavec_func_limit_len[1], length.out=lambdavec_func_limit_len[2])
  }
  if(equal_lambdaf==1){
    funcLambda_inputs = matrix(rep(lambdavec_func,each = numfunc), length(lambdavec_func), numfunc, byrow=T)
  }else{
    funcLambda0 = permutations(length(lambdavec_func), numfunc, repeats.allowed=T)
    funcLambda_inputs = matrix(lambdavec_func[c(funcLambda0)],nrow(funcLambda0),ncol(funcLambda0),byrow=F)
  }

  lambdavec = seq((min(Lambda_limit)),(max(Lambda_limit)),len=Lenlam)
  if(Nq == 1){
    ratios = 1
  }else{
    ratios = seq(1, 0, length.out=subtuning)
  }
  Start = 1
  for(i1 in 1:length(tauvec)){
    tau = tauvec[i1]
    for(ss in 1:length(lambdavec)){
       total = lambdavec[ss]
       subs1 = total*ratios
       subs2 = total*(1-ratios)
       for(subs in 1:length(subs1)){
         temp = cbind(subs1[subs],funcLambda_inputs, subs2[subs]*tau, tau)
         if(Start==1){
           AllTuningMatrix = temp
           Start = 0
         }else{
           AllTuningMatrix = rbind(AllTuningMatrix, temp)
         }
       }
     }
   }
   test = apply(AllTuningMatrix,1,paste0,collapse=":")
   uniq_test = unique(test)
   AllTuningMatrix = AllTuningMatrix[match(uniq_test,test),]
   return(AllTuningMatrix)
}

nonzero = function(xx){
  return(length(which(xx!=0)))
}

Cleaning = function(BetaMatrix, Numitervec, AllTuningMatrix){
   tuningvec = apply(AllTuningMatrix,1,paste0,collapse=":")
   uniqvec = unique(tuningvec)
   mat = match(uniqvec,tuningvec)
   Numitervec = Numitervec[mat]
   BetaMatrix = BetaMatrix[mat,]
   AllTuningMatrix = AllTuningMatrix[mat,]
   NumCounts = apply(BetaMatrix,1,nonzero)
   od = order(NumCounts)
   Numitervec = Numitervec[od]
   BetaMatrix = BetaMatrix[od,]
   AllTuningMatrix = AllTuningMatrix[od,]
   output = list()
   output[[1]] = Numitervec
   output[[2]] = BetaMatrix
   output[[3]] = AllTuningMatrix
   return(output)
}


