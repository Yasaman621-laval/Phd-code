Gen_plinkLDList = function(outnames_vec, Allkeepsnps, rcut = 0.03){
   plinkLDList = list()
   for(ii in 1:length(outnames_vec)){
     ld_outname = outnames_vec[ii]
     lddata = read.table(file=ld_outname,header=T,as.is=T)
     missR = c(NA,"NaN")
     wrm = which(lddata$R %in% missR)
     if(length(wrm)>0){
       lddata = lddata[-wrm,]
     }
     wkeep = which( lddata$SNP_A  %in% Allkeepsnps & lddata$SNP_B %in% Allkeepsnps)
     lddata = lddata[wkeep,]
     wkeep = which( abs(lddata$R) >= rcut)
     lddata = lddata[wkeep,]
     plinkLDList[[ii]] = lddata
     rm(lddata)
   }
   gc()
   return(plinkLDList)
}


Gen_plinkLDList3traits = function(outnames_vec, Allkeepsnps, rcut = 0.3) {
   plinkLDList = list()
   for(ii in 1:3){
     ld_outname = outnames_vec
     lddata = read.table(file=ld_outname,header=T,as.is=T)
     missR = c(NA,"NaN")
     wrm = which(lddata$R %in% missR)
     if(length(wrm)>0){
       lddata = lddata[-wrm,]
     }
     wkeep = which( lddata$SNP_A  %in% Allkeepsnps & lddata$SNP_B %in% Allkeepsnps)
     lddata = lddata[wkeep,]
     wkeep = which( abs(lddata$R) >= rcut)
     lddata = lddata[wkeep,]
     plinkLDList[[ii]] = lddata
     rm(lddata)
   }
   gc()
   return(plinkLDList)
}






es_alpha = function(xx,p){
        vec1 = xx[1:p]
        vec2 = xx[(p+1):length(xx)]
        c(mean(vec1==0 & vec2==0),mean(vec1==0 & vec2!=0),mean(vec1!=0 & vec2==0),mean(vec1!=0 & vec2!=0))
  }

es_alpha_three_pops = function(xx, p) {
  # Split the vector xx into three vectors, each corresponding to one population
  vec1 = xx[1:p]            # First population (1st group of SNPs)
  vec2 = xx[(p+1):(2*p)]    # Second population (2nd group of SNPs)
  vec3 = xx[(2*p+1):(3*p)]  # Third population (3rd group of SNPs)
  
  # Calculate the 8 components for three populations
  pi_0 = mean(vec1 == 0 & vec2 == 0 & vec3 == 0)        # SNPs that are 0 in all populations
  pi_1 = mean(vec1 == 0 & vec2 == 0 & vec3 != 0)        # SNPs that are 0 in first two pops, non-zero in third
  pi_2 = mean(vec1 == 0 & vec2 != 0 & vec3 == 0)        # SNPs that are 0 in first and third pops, non-zero in second
  pi_3 = mean(vec1 == 0 & vec2 != 0 & vec3 != 0)        # SNPs that are 0 in first pop, non-zero in second and third
  pi_12 = mean(vec1 != 0 & vec2 != 0 & vec3 == 0)       # SNPs that are non-zero in first and second pops, 0 in third
  pi_13 = mean(vec1 != 0 & vec2 == 0 & vec3 != 0)       # SNPs that are non-zero in first and third pops, 0 in second
  pi_23 = mean(vec1 == 0 & vec2 != 0 & vec3 != 0)       # SNPs that are non-zero in second and third pops, 0 in first
  pi_123 = mean(vec1 != 0 & vec2 != 0 & vec3 != 0)      # SNPs that are non-zero in all three populations

  # Return the 8 components as a vector
  return(c(pi_0, pi_1, pi_2, pi_3, pi_12, pi_13, pi_23, pi_123))
}




gen_sPostB = function(xx_bb_vec,P, num_alpha){
  uDMatrix = matrix(xx_bb_vec[1:(P*num_alpha)],P,num_alpha,byrow=T)
  sumD = apply(uDMatrix,1,sum)
  uPMatrix = uDMatrix/matrix(rep(sumD,each=num_alpha),P,num_alpha,byrow=T)
  
  wNA = unique(c(which(is.na(sumD)), which(sumD==0)))
  if(length(wNA)>0){
    uPMatrix[wNA,] = 0
  }
  
  bbMatrix = matrix(xx_bb_vec[(P*num_alpha+1):length(xx_bb_vec)],P,num_alpha,byrow=T)

  weightbb = apply(uPMatrix*bbMatrix,1,sum)
  return(weightbb)
}


gen_sPostB2 = function(xx_bb_vec,P, num_alpha){
  uDMatrix = matrix(xx_bb_vec[1:(P*num_alpha)],P,num_alpha,byrow=T)
  sumD = apply(uDMatrix,1,sum)
  uPMatrix = uDMatrix/matrix(rep(sumD,each=num_alpha),P,num_alpha,byrow=T)
  
  wNA = unique(c(which(is.na(sumD)), which(sumD==0)))
  if(length(wNA)>0){
    uPMatrix[wNA,] = 0
  }
  
  #bbMatrix = matrix(xx_bb_vec[(P*num_alpha+1):length(xx_bb_vec)],P,num_alpha,byrow=T)

  #weightbb = apply(uPMatrix*bbMatrix,1,sum)
  #return(weightbb)
  return(uPMatrix)
}





  gen_uprob_whmax = function(xx,P, num_alpha){
        uMatrix = matrix(xx,P,num_alpha,byrow=T)
        apply(uMatrix,1,which.max)
  }

  count_num = function(xx,num){
        sum(xx %in% num)
  }


  gen_Uden = function(xx,P, num_alpha){
    matrix(xx,P,num_alpha,byrow=T)
  }

#Gen_tuningMatrix_RealmixLOG_trans= function(maxL1,maxL2, minL, NumL = 50, subNumL = 50){


    #Lvec = seq(max(maxL1,maxL2),minL,length.out=(NumL+1))[1:NumL]
    #ratios = seq(1,0,length.out=subNumL)
  #  for(ss in 1:NumL){
    #  total = Lvec[ss]
   #   subs1 = total*ratios
      #subs2 = total*(1-ratios)
      #if(ss == 1){
       # TuningsMatrix = cbind(subs1,subs2)
      #}else{
      #  TuningsMatrix = rbind(TuningsMatrix,cbind(subs1,subs2))
    #  }
   # }
  #  return(TuningsMatrix)
#  }



Gen_tuningMatrix_RealmixLOG_trans= function(maxL, minL, NumL = 50, subNumL = 50){


    Lvec = seq(maxL,minL,length.out=(NumL+1))[1:NumL]
    ratios = seq(1,0,length.out=subNumL)
    for(ss in 1:NumL){
      total = Lvec[ss]
      subs1 = total*ratios
      subs2 = total*(1-ratios)
      if(ss == 1){
        TuningsMatrix = cbind(subs1,subs2)
      }else{
        TuningsMatrix = rbind(TuningsMatrix,cbind(subs1,subs2))
      }
    }
    return(TuningsMatrix)
  }



  
 Gen_tuningMatrix_RealmixLOG_trans3= function(maxL1,maxL2,maxL3, minL, NumL = 50, subNumL = 50){


    Lvec = seq(max(maxL1,maxL2,maxL3),minL,length.out=(NumL+1))[1:NumL]
    ratios = seq(1,0,length.out=subNumL)
    for(ss in 1:NumL){
      total = Lvec[ss]
      subs1 = total*ratios
      subs2 = total*(1-ratios)
      if(ss == 1){
        TuningsMatrix = cbind(subs1,subs2)
      }else{
        TuningsMatrix = rbind(TuningsMatrix,cbind(subs1,subs2))
      }
    }
    return(TuningsMatrix)
  } 
  


Gen_tuningMatrix_mixLOG = function(maxL, minL, NumL = 50, subNumL = 50){


    Lvec = seq(maxL,minL,length.out=(NumL+1))[1:NumL]
    ratios = seq(1,0,length.out=subNumL)
    for(ss in 1:NumL){
      total = Lvec[ss]
      subs1 = total*ratios
      subs2 = total*(1-ratios)
      if(ss == 1){
        TuningsMatrix = cbind(subs1,0,subs2,1)
      }else{
        TuningsMatrix = rbind(TuningsMatrix,cbind(subs1,0,subs2,1))
      }
    }
    return(TuningsMatrix)
  }


Gen_tuningMatrix_LOGfixweight = function(maxL, minL, NumL = 400){


    Lvec = seq(maxL,minL,length.out=NumL)
    TuningsMatrix = cbind(0,1,Lvec,1)
    return(TuningsMatrix)
  }



Gen_tuningMatrix_mixLOGfixweight = function(maxL, minL, NumL = 50, subNumL = 50){


    Lvec = seq(maxL,minL,length.out=(NumL+1))[2:(NumL+1)]

    for(ss in 1:NumL){
      diff = maxL - Lvec[ss]
      subs = seq(diff,0,length.out = subNumL)
      if(ss == 1){
        TuningsMatrix = cbind(Lvec[ss],0,subs,1)
      }else{
        TuningsMatrix = rbind(TuningsMatrix,cbind(Lvec[ss],0,subs,1))
      }
    }
    return(TuningsMatrix)
  }



Gen_tuningMatrix_LOG_trans = function(savefile0,savef0,perR = 0.8, NumL = 100){
    load(file=savefile0)
    
    InitTunnings = saveoutput[[7]]
    whmax = which.max(PreR2[1,])
    qPreR2 = quantile(PreR2[1,],perR,na.rm=T)

    wwcs = which(PreR2[1,] >= qPreR2)
    sww = min(wwcs)
    mww = max(wwcs)

    maxL = max(InitTunnings[c(sww,mww),1])
    minL = min(InitTunnings[c(sww,mww),1])

    Lvec = seq(maxL,minL,length.out=(NumL+1))[2:(NumL+1)]
    
    TuningsMatrix = cbind(0,Lvec)
    
    return(TuningsMatrix)
  }



Gen_tuningMatrix_DLasso_trans = function(savefile0, savef0, perR = 0.7, NumL = 100){
    load(file=savefile0)
    
    InitTunnings = saveoutput[[7]]
    whmax = which.max(PreR2[1,])
    qPreR2 = quantile(PreR2[1,],perR,na.rm=T)

    wwcs = which(PreR2[1,] >= qPreR2)
    sww = min(wwcs)
    mww = max(wwcs)

    maxL = max(InitTunnings[c(sww,mww),1])
    minL = min(InitTunnings[c(sww,mww),1])

    Lvec = seq(maxL,minL,length.out=NumL)
    
    return(Lvec)
  }


Gen_tuningMatrix_mixLOG_trans = function(savefile0,savef0,perR = 0.8, NumL = 50, subNumL = 50){
    load(file=savefile0)
    
    InitTunnings = saveoutput[[7]]
    whmax = which.max(PreR2[1,])
    qPreR2 = quantile(PreR2[1,],perR,na.rm=T)

    wwcs = which(PreR2[1,] >= qPreR2)
    sww = min(wwcs)
    mww = max(wwcs)

    maxL = max(InitTunnings[c(sww,mww),1])
    minL = min(InitTunnings[c(sww,mww),1])

    Lvec = seq(maxL,minL,length.out=(NumL+1))[2:(NumL+1)]
    subNumL = 50
    for(ss in 1:NumL){
      diff = maxL - Lvec[ss]
      subs = seq(0,diff,length.out = (ss+1))
      if(ss == 1){
        TuningsMatrix = cbind(Lvec[ss],subs)
      }else{
        TuningsMatrix = rbind(TuningsMatrix,cbind(Lvec[ss],subs))
      }
    }
    return(TuningsMatrix)
  }


Gen_tau = function(tauuse, vecs){
   if(tauuse <= 1){
     tau = quantile(vecs,tauuse)
   }else{
     tau = quantile(vecs,1)*tauuse
   }
  }

   Get_lambdasL = function(savefile10)
{
   output1 = Gen_One_BetaMatrix(savefile10, 1,1)
   nonzero1 = apply(output1[[1]],1,nonzero)
   minZR1 = min(nonzero1[nonzero1>0],na.rm=T)
   maxZR1 = max(nonzero1[nonzero1>0],na.rm=T)
   w1 = max(which(nonzero1==minZR1))
   w2 = min(which(nonzero1==maxZR1))

   total = output1[[2]][w1,1]
   lowv = output1[[2]][w2,1]
   rm(output1,nonzero1)
   return(c(total,lowv))
}



Gen_One_BetaMatrix_saveoutput = function(saveoutput, groupnum, gindex){

    iterindex = saveoutput[[1]]
    ncol_BetaMatrix = saveoutput[[2]]
    nrow_BetaMatrix = saveoutput[[3]]
    wpos = saveoutput[[4]]
    vals = saveoutput[[5]]
    AllBetaMatrix = matrix(0,nrow_BetaMatrix,ncol_BetaMatrix)
    AllBetaMatrix[wpos] = vals
    rm(wpos,vals)

    P = length(saveoutput[[6]])
    

    gindex1vec =seq(1,(groupnum*P),by=P)
    gindex2vec =seq(P,(groupnum*P),by=P)
    gindex1 = gindex1vec[gindex]
    gindex2 = gindex2vec[gindex]

    AllBetaMatrix = AllBetaMatrix[, gindex1:gindex2, drop=F]
    colnames(AllBetaMatrix) = saveoutput[[6]]
    
    tuningMatrix = saveoutput[[7]]
    rm(saveoutput)
    gc()
    output = list()
    output[[1]] = AllBetaMatrix
    output[[2]] = tuningMatrix
    output[[3]] = iterindex
    return(output)
}


Gen_One_BetaMatrix = function(savefile1, groupnum, gindex){

    load(file=savefile1)

    iterindex = saveoutput[[1]]
    ncol_BetaMatrix = saveoutput[[2]]
    nrow_BetaMatrix = saveoutput[[3]]
    wpos = saveoutput[[4]]
    vals = saveoutput[[5]]
    AllBetaMatrix = matrix(0,nrow_BetaMatrix,ncol_BetaMatrix)
    AllBetaMatrix[wpos] = vals
    rm(wpos,vals)

    P = length(saveoutput[[6]])
    

    gindex1vec =seq(1,(groupnum*P),by=P)
    gindex2vec =seq(P,(groupnum*P),by=P)
    gindex1 = gindex1vec[gindex]
    gindex2 = gindex2vec[gindex]

    AllBetaMatrix = AllBetaMatrix[, gindex1:gindex2, drop=F]
    colnames(AllBetaMatrix) = saveoutput[[6]]
    
    tuningMatrix = saveoutput[[7]]
    rm(saveoutput)
    gc()
    output = list()
    output[[1]] = AllBetaMatrix
    output[[2]] = tuningMatrix
    output[[3]] = iterindex
    return(output)
}



Gen_One_BetaMatrix_3 = function(savefile1, groupnum, gindex){

    load(file=savefile1)

    iterindex = saveoutput[[1]]
    ncol_BetaMatrix = saveoutput[[2]]
    nrow_BetaMatrix = saveoutput[[3]]
    wpos = saveoutput[[4]]
    vals = saveoutput[[5]]
    AllBetaMatrix = matrix(0,nrow_BetaMatrix,ncol_BetaMatrix)
    AllBetaMatrix[wpos] = vals
    rm(wpos,vals)

    P = length(saveoutput[[6]])
    

gindex1vec = seq(1, (groupnum * P), by = P)  # First group (population 1)
gindex2vec = seq(P + 1, (groupnum * P), by = P)  # Second group (population 2)
gindex3vec = seq(2 * P + 1, (groupnum * P), by = P)  # Third group (population 3)

# Select the indices for each population based on a specific group index (gindex)
gindex1 = gindex1vec[gindex]
gindex2 = gindex2vec[gindex]
gindex3 = gindex3vec[gindex]

    AllBetaMatrix = AllBetaMatrix[, gindex1:gindex3, drop=F]
    colnames(AllBetaMatrix) = saveoutput[[6]]
    
    tuningMatrix = saveoutput[[7]]
    rm(saveoutput)
    gc()
    output = list()
    output[[1]] = AllBetaMatrix
    output[[2]] = tuningMatrix
    output[[3]] = iterindex
    return(output)
}






Gen_outnames_mafnames = function(dirPlinkFormat, popvec, usedtrait){
    usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
    usedtraitIndex = as.numeric(usedtraitsvec)
    outnames_vec = c()
    mafnames_vec = c()
    for(ii in 1:length(usedtraitIndex)){
      iiIndex = usedtraitIndex[ii]
      pop = popvec[iiIndex]
      outname= paste0(dirPlinkFormat,pop,"chr",chr,"bedformat")
      ld_outname = paste0(outname,".ld")
      outnames_vec = c(outnames_vec,ld_outname)
      maffile =  paste0(dirPlinkFormat,pop,"chr",chr,"bedformat.frq")
      mafnames_vec = c(mafnames_vec,maffile)
   }
   output = list()
   output[[1]] = outnames_vec
   output[[2]] = mafnames_vec
   return(output)
 }

abs_diff = function(btemp){
    nq = length(btemp)
    temp = 0
    for(i in 1:(nq-1)){
      for(j in (i+1):nq){
        temp = temp + abs(btemp[i]-btemp[j])
      }
    }
    return(temp)
}
Get_saveoutput_trans = function(output, savefile1){
      saveoutput = list()
      saveoutput[[1]] = output$Numitervec
      saveoutput[[2]] = ncol(output$BetaMatrix)
      saveoutput[[3]] = nrow(output$BetaMatrix)
      saveoutput[[4]] = which(output$BetaMatrix!=0)
      wpos = which(output$BetaMatrix!=0)
      saveoutput[[5]] = output$BetaMatrix[wpos]
      saveoutput[[6]] = rownames(summaryZ)
      saveoutput[[7]] = output[[3]]
      saveoutput[[8]] = output$Dfq1
      save(saveoutput,file=savefile1)
      rm(output,saveoutput)
}
Get_summaryBetas = function(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr){

    GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
    if(!is.null(chr)){
      wk = which(GWASbeta$CHR==chr)
      GWASbeta = GWASbeta[wk,]
    }
    rownames(GWASbeta) = GWASbeta$SNP
    Allkeepsnps = GWASbeta$SNP
    od = order(GWASbeta[,paste0("p",mainindex)],decreasing=F)
    GWASbeta = GWASbeta[od,]
    namesBeta = paste0("b",usedtraitsvec)
    summaryBetas = GWASbeta[,match(namesBeta, colnames(GWASbeta)),drop=F]
    rownames(summaryBetas) = GWASbeta$SNP
    summaryBetas = check_summaryZ(summaryBetas)

    sumbetas = apply(abs(summaryBetas),1,sum)
    sumdiffbetas = apply(summaryBetas,1,abs_diff)
    
    sumbetas_after = abs(apply(summaryBetas,1,sum))


    namesBeta = paste0("SE",usedtraitsvec)
    SDvec = GWASbeta[,match(namesBeta, colnames(GWASbeta)),drop=F]
    mat = match(rownames(summaryBetas), rownames(SDvec))
    SDvec = SDvec[mat,,drop=F]

    rm(GWASbeta)


    output = list()
    output[[1]] = summaryBetas
    output[[2]] = Allkeepsnps
    output[[3]] = sumbetas
    output[[4]] = sumdiffbetas
    output[[5]] = SDvec
    output[[6]] = sumbetas_after
    rm(summaryBetas,Allkeepsnps, SDvec)
    return(output)
 }


#Get_summaryZ = function(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr){
    #GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)

#if (!"A2" %in% names(GWASbeta)) {
 # GWASbeta <- GWASbeta %>%
  #  mutate(
      # Extract second allele (A2) from pattern like "1:752721_A_G"
   #   A2 = str_extract(SNP, "(?<=_[A-Za-z])_[A-Za-z]+$") %>%
     #      str_remove("^_")
    #) %>%
   # select(
   #   CHR, SNP, A1, A2,
    #  Zobs1, b1, SE1, p1,
    #  Zobs2, b2, SE2, p2,
    #  Zobs3, b3, SE3, p3
   #$ )
  
#} else {
  #GWASbeta <- GWASbeta %>%
   # select(
    #  CHR, SNP, A1, A2,
     # Zobs1, b1, SE1, p1,
     # Zobs2, b2, SE2, p2,
     # Zobs3, b3, SE3, p3
   # )
#}

    #if(!is.null(chr)){
    #  wk = which(GWASbeta$CHR==chr)
    #  GWASbeta = GWASbeta[wk,]
   # }
   # rownames(GWASbeta) = GWASbeta$SNP
   # Allkeepsnps = GWASbeta$SNP
  #  od = order(GWASbeta[,paste0("p",mainindex)],decreasing=F)
  #  GWASbeta = GWASbeta[od,]
   # namesBeta = paste0("Zobs",usedtraitsvec)
   # summaryZ = GWASbeta[,match(namesBeta, colnames(GWASbeta)),drop=F]
  #  rownames(summaryZ) = GWASbeta$SNP
  #  summaryZ = check_summaryZ(summaryZ)
   # rm(GWASbeta)
   # suboutput = Get_sum_diff_betas(Nvec, summaryZ)
   
   # sumbetas = suboutput[[1]]
   # sumdiffbetas = suboutput[[2]]

   # output = list()
   # output[[1]] = summaryZ
   # output[[2]] = Allkeepsnps
   # output[[3]] = sumbetas
   # output[[4]] = sumdiffbetas
   # output[[5]] = suboutput[[3]] #abs(summaryBetas)
  # output[[6]] = suboutput[[4]] #summaryBetas
    #output[[7]] = suboutput[[5]] #sumbetas_after
    #rm(summaryZ,Allkeepsnps, suboutput)
   # return(output)
#}

Get_summaryZ = function(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr){
    GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
    if(!is.null(chr)){
      wk = which(GWASbeta$CHR==chr)
      GWASbeta = GWASbeta[wk,]
    }
    rownames(GWASbeta) = GWASbeta$SNP
    Allkeepsnps = GWASbeta$SNP
    od = order(GWASbeta[,paste0("p",mainindex)],decreasing=F)
    GWASbeta = GWASbeta[od,]
    namesBeta = paste0("Zobs",usedtraitsvec)
    summaryZ = GWASbeta[,match(namesBeta, colnames(GWASbeta)),drop=F]
    rownames(summaryZ) = GWASbeta$SNP
    summaryZ = check_summaryZ(summaryZ)
    rm(GWASbeta)
    suboutput = Get_sum_diff_betas(Nvec, summaryZ)
   
    sumbetas = suboutput[[1]]
    sumdiffbetas = suboutput[[2]]

    output = list()
    output[[1]] = summaryZ
    output[[2]] = Allkeepsnps
    output[[3]] = sumbetas
    output[[4]] = sumdiffbetas
    output[[5]] = suboutput[[3]] #abs(summaryBetas)
    output[[6]] = suboutput[[4]] #summaryBetas
    output[[7]] = suboutput[[5]] #sumbetas_after
    rm(summaryZ,Allkeepsnps, suboutput)
    return(output)
}



Get_summary_3 = function(GWASbetafile, Nvec, chr){
    GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
    GWASbeta$CHR <- sapply(strsplit(GWASbeta$SNP, ":"), `[`, 1)

    if(!is.null(chr)){
      wk = which(GWASbeta$CHR==chr)
      GWASbeta = GWASbeta[wk,]
    }
    rownames(GWASbeta) = GWASbeta$SNP
    Allkeepsnps = GWASbeta$SNP
   od = order(GWASbeta[["p1"]], GWASbeta[["p2"]], decreasing = FALSE)

    GWASbeta = GWASbeta[od,]
    namesBeta = paste0("Zobs",c(1,2,3))
    summaryZ = GWASbeta[,match(namesBeta, colnames(GWASbeta)),drop=F]
    rownames(summaryZ) = GWASbeta$SNP
    summaryZ = check_summaryZ(summaryZ)
    rm(GWASbeta)
    suboutput = Get_sum_diff_betas(Nvec, summaryZ)
   
    sumbetas = suboutput[[1]]
    sumdiffbetas = suboutput[[2]]

    output = list()
    output[[1]] = summaryZ
    output[[2]] = Allkeepsnps
    output[[3]] = sumbetas
    output[[4]] = sumdiffbetas
    output[[5]] = suboutput[[3]] #abs(summaryBetas)
    output[[6]] = suboutput[[4]] #summaryBetas
    output[[7]] = suboutput[[5]] #sumbetas_after
    rm(summaryZ,Allkeepsnps, suboutput)
    return(output)
}


Get_summaryBetas_only = function(Nvec, summaryZ){
  Nq = length(Nvec)
  summaryBetas = matrix(0, nrow(summaryZ), Nq)
  rownames(summaryBetas) = rownames(summaryZ)
  SDvec = matrix(0, nrow(summaryZ), Nq)
  for(ii in 1:Nq){
    summaryBetas[,ii] = summaryZ[,ii]/sqrt(Nvec[ii])
    SDvec[,ii] = 1/sqrt(Nvec[ii])
  }
  return(summaryBetas)
}

Get_sum_diff_betas = function(Nvec, summaryZ){
    Nq = length(Nvec)
    summaryBetas = matrix(0, nrow(summaryZ), Nq)
    rownames(summaryBetas) = rownames(summaryZ)
    SDvec = matrix(0, nrow(summaryZ), Nq)

    for(ii in 1:Nq){
      summaryBetas[,ii] = summaryZ[,ii]/sqrt(Nvec[ii])
      SDvec[,ii] = 1/sqrt(Nvec[ii])
    }
    sumbetas = apply(abs(summaryBetas),1,sum)
    sumdiffbetas = apply(summaryBetas,1,abs_diff)

    sumbetas_after = abs(apply(summaryBetas,1,sum))
    output = list()
    output[[1]] = sumbetas
    output[[2]] = sumdiffbetas
    output[[3]] = abs(summaryBetas)
    output[[4]] = summaryBetas
    output[[5]] = sumbetas_after
    return(output)
  }

check_summaryZ = function(summaryZ){
    for(ii in 1:ncol(summaryZ)){
      wna = which(is.na(summaryZ[,ii]) | summaryZ[,ii]=="Inf" | summaryZ[,ii]=="NaN")
      if(length(wna)>0){
        summaryZ[wna,ii] = 0
      }
    }
    return(summaryZ)
}

Gen_ThetaMatrix = function(savefile1){

    load(file=savefile1)

    iterindex = saveoutput[[1]]
    nrow_ThetaMatrix = saveoutput[[9]]
    ncol_ThetaMatrix = saveoutput[[10]]
    wpos = saveoutput[[11]]
    vals = saveoutput[[12]]

    ThetaMatrix = matrix(0,nrow_ThetaMatrix,ncol_ThetaMatrix)
    ThetaMatrix[wpos] = vals
    colnames(ThetaMatrix) = saveoutput[[6]]
    rm(wpos,vals)

    P = length(saveoutput[[6]])
    
    wk = which(iterindex > 0)
    ThetaMatrix = ThetaMatrix[wk,]
    tuningMatrix = saveoutput[[7]]
    tuningMatrix = tuningMatrix[wk,,drop=F]
    rm(saveoutput)
    gc()
    output = list()
    output[[1]] = ThetaMatrix
    output[[2]] = tuningMatrix
    return(output)
}
Gen_simpleTrans_input = function(usedtrait1, TrainingNsam, Zscale, dirSimuData, gwasfilename, mafnames_vec1,low0_single, top0_single, chr){
     usedtraitsvec = unlist(strsplit(as.character(usedtrait1),split=","))
     usedtraitIndex = as.numeric(usedtraitsvec)
     Nvec = TrainingNsam[usedtraitIndex]
     GWASbetafile = paste0(dirSimuData,gwasfilename)
     GWASbeta = read.table(file=paste0(dirSimuData,gwasfilename), sep = "\t", header=T, as.is=T)
     if(!is.null(chr)){
      wk = which(GWASbeta$CHR==chr)
      GWASbeta = GWASbeta[wk,]
     }
     rownames(GWASbeta) = GWASbeta$SNP
     Allkeepsnps = GWASbeta$SNP
     usedvec = c(paste0("b",usedtraitIndex),paste0("SE",usedtraitIndex))
     inputmatrix = GWASbeta[,match(usedvec,colnames(GWASbeta) )]
     Meta = t(apply(inputmatrix,1,meta.F,2))
     maxZ = max(abs(Meta[,3]),na.rm=T)
     if(Zscale==1){
       summaryZ = matrix(Meta[,3],nrow(Meta),1)
       rownames(summaryZ) = GWASbeta$SNP
       summaryZ = check_summaryZ(summaryZ)
       NamesSNP = rownames(summaryZ)
       #suboutput = Get_sum_diff_betas(Nvec, summaryZ)
       #sumbetas = suboutput[[1]]
       #sumdiffbetas = suboutput[[2]]
       summaryBetas =  Get_summaryBetas_only(sum(Nvec), summaryZ)
       #sumbetas_after = suboutput[[5]]
       SSvec = SDvec = NULL
     }else{
       summaryZoutput = Get_summaryBetas(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr)
       SDvec = summaryZoutput[[5]]
       mat = match(rownames(SDvec),Allkeepsnps)
       SDvec = SDvec[mat,,drop=F]
       summaryZ = NULL
       summaryBetas = matrix(Meta[,1],nrow(Meta),1)
       rownames(summaryBetas) = GWASbeta$SNP
       SSvec = matrix(0,nrow(summaryBetas),ncol(summaryBetas))
       rownames(SSvec) = rownames(summaryBetas)
       NamesSNP = rownames(summaryBetas)

       SSvec0 = matrix(0,nrow(summaryBetas),2)
       rownames(SSvec0) = rownames(summaryBetas)

       for(ii in 1:length(mafnames_vec1)){
         MAFdata0 = read.table(file=mafnames_vec1[ii],as.is=T,header=T)
         mat = match(rownames(SSvec),MAFdata0$SNP)
         SSvec0[,ii] = 2*MAFdata0$MAF[mat]*(1-MAFdata0$MAF[mat])*Nvec[ii]
         rm(MAFdata0)
       }
       SSvec[,1] = apply(SSvec0,1,sum)
       rm(SSvec0)
     }
     rm(GWASbeta)
     gc()
     if(is.null(top0_single)){
       top0_single = maxZ
     }
     if(Zscale==0){
       topvs = top0_single*median(Meta[,2],na.rm=T)*median(c(SSvec))
       lowvs = low0_single*min(Meta[,2],na.rm=T)*min(SSvec)
     }else{
       topvs = top0_single
       lowvs = low0_single
     }
     rm(Meta)
     output = list()
     output[[1]] = summaryZ
     output[[2]] = summaryBetas
     output[[3]] = SSvec
     output[[4]] = SDvec
     output[[5]] = topvs
     output[[6]] = lowvs
     output[[7]] = Nvec
     output[[8]] = Allkeepsnps
     return(output)
}

Get_inputs = function(Zscale, usedtrait, TrainingNsam, dirSimuData,outnames_vec, RupperVal, gwasfilename, mafnames_vec, penalty, low1, top1){
     sumbetas = NULL
     sumdiffbetas = NULL
     usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
     usedtraitIndex = as.numeric(usedtraitsvec)
     Nvec = TrainingNsam[usedtraitIndex]
     GWASbetafile = paste0(dirSimuData,gwasfilename)
     if(Zscale==1){
      summaryZoutput = Get_summaryZ(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr)
      summaryZ = summaryZoutput[[1]]
      Allkeepsnps = summaryZoutput[[2]]
      sumbetas = summaryZoutput[[3]]
      sumdiffbetas = summaryZoutput[[4]]
      #absbeta = summaryZoutput[[5]]
      summaryBetas = summaryZoutput[[6]]
      sumbetas_after = summaryZoutput[[7]]
      SSvec = NULL
      namesSNP = rownames(summaryZ)
      if(!is.null(RupperVal)){
        RupperVal = max(1/sqrt(Nvec))*RupperVal
      }
      topv = top1
      lowv = low1
    }else{
      summaryZoutput = Get_summaryBetas(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr)
      summaryBetas = summaryZoutput[[1]]
      Allkeepsnps = summaryZoutput[[2]]
      sumbetas = summaryZoutput[[3]]
      sumdiffbetas = summaryZoutput[[4]]
      SDvec = summaryZoutput[[5]]
      sumbetas_after = summaryZoutput[[6]]
      if(!is.null(RupperVal)){
        RupperVal = max(unlist(SDvec))*RupperVal
      }
      summaryZ = NULL
      SSvec = matrix(0,nrow(summaryBetas),ncol(summaryBetas))
      rownames(SSvec) = rownames(summaryBetas)
      namesSNP = rownames(summaryBetas)

      for(ii in 1:length(mafnames_vec)){
        MAFdata0 = read.table(file=mafnames_vec[ii],as.is=T,header=T)
        mat = match(rownames(SSvec),MAFdata0$SNP)
        SSvec[,ii] = 2*MAFdata0$MAF[mat]*(1-MAFdata0$MAF[mat])*Nvec[ii]
        rm(MAFdata0)
        }
        topv = top1*median(c(unlist(SDvec[,1])))*median(c(SSvec[,1]))
        lowv = low1*min(c(unlist(SDvec[,1])))*min(SSvec[,1])
      }
    output = list()
    output[[1]] = summaryZ
    output[[2]] = summaryBetas
    output[[3]] = SSvec
    output[[4]] = namesSNP
    output[[5]] = RupperVal
    output[[6]] = sumbetas
    output[[7]] = sumdiffbetas
    output[[8]] = topv
    output[[9]] = lowv
    output[[10]] = Nvec
    output[[11]] = sumbetas_after
    return(output)
  }

Gen_tuningMatrix = function(penalty, lowv , topv, totaltuning1, totaltuning2, tau, totalv){
    if(penalty=="ratiow" || penalty=="mixLOG_after"|| penalty=="mixLOG_afterS"){tau0 = tau}
    if(penalty=="ratiowS"){tau0 = 1}
    if(penalty=="mixLOG_afterS"){
      tuningMatrix = cbind(0,seq((lowv),(topv),len=totaltuning2)*tau0)
    }else{
      lambda0vec = seq((lowv),(topv),len=totaltuning1)
      for(ll in 1:length(lambda0vec)){
        lam = lambda0vec[ll]
        diff = totalv - lam
        if(ll==1){
          tuningMatrix = cbind(rep(lam,totaltuning2),seq(0,diff,length.out=totaltuning2)*tau0)
        }else{
          tuningMatrix = rbind(tuningMatrix,cbind(rep(lam,totaltuning2),seq(0,diff,length.out=totaltuning2)*tau0))
        }
      }
    }
    
    return(tuningMatrix)
}

Get_saveoutput_trans = function(output, savefile1, namesSNP){
      saveoutput = list()
      saveoutput[[1]] = output$Numitervec
      saveoutput[[2]] = ncol(output$BetaMatrix)
      saveoutput[[3]] = nrow(output$BetaMatrix)
      saveoutput[[4]] = which(output$BetaMatrix!=0)
      wpos = which(output$BetaMatrix!=0)
      saveoutput[[5]] = output$BetaMatrix[wpos]
      saveoutput[[6]] = namesSNP
      saveoutput[[7]] = output[[3]]
      saveoutput[[8]] = output$Dfq1
      save(saveoutput,file=savefile1)
      rm(output,saveoutput)
}

Gen_One_BetaMatrix_noteskip = function(savefile1, groupnum, gindex){
  load(file=savefile1)

  iterindex = saveoutput[[1]]
  ncol_BetaMatrix = saveoutput[[2]]
  nrow_BetaMatrix = saveoutput[[3]]
  wpos = saveoutput[[4]]
  vals = saveoutput[[5]]
  AllBetaMatrix = matrix(0,nrow_BetaMatrix,ncol_BetaMatrix)
  AllBetaMatrix[wpos] = vals
  rm(wpos,vals)

  P = length(saveoutput[[6]])

  gindex1vec = seq(1, groupnum*P, by=P)
  gindex2vec = seq(P, groupnum*P, by=P)
  gindex1 = gindex1vec[gindex]
  gindex2 = gindex2vec[gindex]

  AllBetaMatrix = AllBetaMatrix[, gindex1:gindex2, drop=F]
  colnames(AllBetaMatrix) = saveoutput[[6]]

  tuningMatrix = saveoutput[[7]]

  wpos2 = saveoutput[[8]]
  noteskipMatrix = matrix(0,nrow_BetaMatrix,ncol_BetaMatrix)
  vals2 = saveoutput[[9]]
  noteskipMatrix[wpos2] = vals2

  #colnames(noteskipMatrix) = colnames(AllBetaMatrix)

  rm(saveoutput)
  gc()
  output = list()
  output[[1]] = AllBetaMatrix
  output[[2]] = tuningMatrix
  output[[3]] = iterindex
  output[[4]] = noteskipMatrix
  return(output)
}

   filter_Beta = function(Beta, tunings, iterinfo, noteskipMatrix){
  wk = which(iterinfo > 0)
  Beta = Beta[wk,]
  tunings = tunings[wk,]
  iterinfo = iterinfo[wk]
  noteskipMatrix=noteskipMatrix[wk,]
  
  nonzero_test = apply(Beta, 1, nonzero)
  wrm = which(nonzero_test==0)
  if(length(wrm)>0){
    Beta = Beta[-wrm,]
    tunings = tunings[-wrm,]
    iterinfo = iterinfo[-wrm]
    noteskipMatrix = noteskipMatrix[-wrm,]
  }
  
  filteroutput = list()
  filteroutput[[1]] = Beta
  filteroutput[[2]] = tunings
  filteroutput[[3]] = iterinfo
  filteroutput[[4]] =  noteskipMatrix
  return(filteroutput)
} 
