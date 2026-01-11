
  meta.F <- function(inputs, Index){
    b.est = inputs[1:Index]
    se = inputs[(Index+1):length(inputs)]
    #returns inverse-variance weighted meta-analysis estimate, SE and P-value.
    b.F = sum(b.est / se^2) / sum(1 / se^2)
    se.F = 1 / sqrt(sum(1 / se^2))
    Z.F = b.F/se.F
    return(c(b.F, se.F, Z.F))
  }

nullRsquare = function(score,pheno){
  N = length(score)
  sampleIndex = sample(1:N,N,replace=T)
  return(cor(score,pheno[sampleIndex])^2)
}





getOriScale = function(Beta, AF){
	    Beta/sqrt(2*AF*(1-AF))
	  }
nonzero = function(xx){
      length(which(xx!=0))
}
scaleBetas = function(inputBetas, MAFdata){
	  wkeep = which(names(inputBetas) %in% MAFdata$SNP)
	  inputBetas = inputBetas[wkeep]
	  
	  mat = match(names(inputBetas), MAFdata$SNP)
	  MAFdata = MAFdata[mat,]
	  
	  scaledBeta  = inputBetas*sqrt(2*MAFdata$AF*(1-MAFdata$AF))
	  
	  wzero = which(is.na(scaledBeta))
	  if(length(wzero)>0){
	    scaledBetas[wzero] = 0
	  }
	  return(scaledBeta)
	}


Generate.Data.inLD <- function(plink.file.discovery, plink.file.validation, plink.freq.file, plink.r.file,
 hsq = c(0.2, 0.3), r2.range = r2.range, n.causal.snp, RatioSNPinLD = c(1/4, 1/3, 1/2),numsim=20){

    r2.rangeName = paste(r2.range,collapse="-")
	n.causal.snpName = paste(n.causal.snp,collapse="-")
	hsq.snpName = paste(hsq,collapse="-")
    
	#system(paste("cp ",plink.file.validation, ".??? ", sub.folder, "/", sep=""));
    #setwd(sub.folder);
    #x = read.table(file= plink.freq.file, skip=1, as.is=T)[,c(2,3,5)];
    MAF = read.table(file= plink.freq.file, skip=1, as.is=T)[,c(1,2,3,5)];	 
	n.indep.snps = nrow(MAF) 
	allsnpnames = MAF[,2]
	
    lddata = read.table(file=plink.r.file, header=T,as.is=T)
	lddata = lddata[which(lddata$R>min(sqrt(r2.range)) & lddata$R <max(sqrt(r2.range))),]
	snpused = unique(lddata$SNP_A)
    nsnp = length(snpused)

	numLDpairs = ceiling(sum(n.causal.snp)*RatioSNPinLD/2)
    
	for(nn in 1:length(RatioSNPinLD)){
	  pair = numLDpairs[nn]
	  findex = sample(c(1:nsnp),pair, replace=F)
	  namesnp1 = snpused[findex]
	  wrm = which(lddata$SNP_B %in% namesnp1)
	  
	  if(length(wrm)>0){
	    tmplddata = lddata[-wrm,]
	  }else{
	    tmplddata = lddata
	  }
	  
	  ssnps = unique(tmplddata$SNP_B)
	  sindex = sample(c(1:length(ssnps)),pair, replace=F)
	  
	  selectsnps = c(namesnp1,ssnps[sindex])
	  if(length(unique(selectsnps))!=2*pair){stop("error in selecting snps")}
	  rm(namesnp1,ssnps,tmplddata)
	  
	  remainsnps = allsnpnames[-which(allsnpnames %in% selectsnps)]
	  rindex = sample(c(1:length(remainsnps)),(sum(n.causal.snp)-length(selectsnps)),replace=F)
	  remainselect = remainsnps[rindex]
	  
	  
	  # simulate MAF-adjusted beta for causal SNPs
      beta = rep(0, n.indep.snps);
	  betatrue = c()
	  for(k in 1:length(hsq)){
     	  betatrue = c(rnorm(n.causal.snp[k])*sqrt(hsq[k]/n.causal.snp[k]), betatrue)
	  }
      beta[which(allsnpnames %in% c(selectsnps,remainselect))] = betatrue 
	  
	  savefolder = paste0(dir0,dataname,"_RatioSNPinLD",as.character(signif(RatioSNPinLD[nn],2)),"-r2.range",r2.rangeName,
	  "-n.causal.snp",n.causal.snpName,"-hsq.snp",hsq.snpName,"/")	 
      system(paste("mkdir", savefolder))

	  # Recover orginal beta coefficients by adjusting maf.
      score = cbind(MAF[,1:3],beta/sqrt(2*MAF[,4]*(1-MAF[,4])))[beta!=0,];
      write.table(score,file=paste0(savefolder,"score"),row.names=F,col.names=F,quote=F)
      
	  setwd(savefolder)
      #command = paste("plink --bfile ../", plink.file.discovery, " --noweb --score score --out temp", sep="")
	  command = paste("plink --bfile ", plink.file.discovery, " --noweb --score score --out temp", sep="")
      system(command);
    
	
	  for(sim in 1:numsim){
	    y = read.table("temp.profile", as.is=T, skip=1)[, c(1,2,6)];
        y[,3] = y[,3]/sqrt(var(y[,3])) * sqrt(sum(hsq)) + rnorm(nrow(y))*sqrt(1-sum(hsq));
        write.table(y, file=paste0("discovery_sim",sim,".pheno"), row.names=F, col.names=F, quote=F);
      }
      #command = paste("plink --bfile ../", plink.file.validation, " --noweb --score score --out temp.validation", sep="");
	  command = paste("plink --bfile ", plink.file.validation, " --noweb --score score --out temp.validation", sep="");
      system(command);
	  
      for(sim in 1:numsim){
	    y = read.table("temp.validation.profile", as.is=T, skip=1)[, c(1,2,6)];
        y[,3] = y[,3]/sqrt(var(y[,3])) * sqrt(sum(hsq)) + rnorm(nrow(y))*sqrt(1-sum(hsq));
        write.table(y, file=paste0("validation_sim",sim,".pheno"), row.names=F, col.names=F, quote=F);
      }
      #system(paste("/data/chent5/script/Match.Two.File.py ", plink.file.validation, ".fam 1 validation.pheno 1 out", sep=""));
      #system(paste("awk '{print $1,$2,$3,$4,$5,$9}' out > ", plink.file.validation, ".fam", sep=""));
	}  
}

	
#input file: SNP (clumped), reference allele, beta (uncorrected for Winnner's curse), P-value

Create.Score.File.For.Plink.LASSO <- function(input.file, p.Threshold.file = "/data/jianxins/GWAS.Risk.Prediction/Shared.Files/polygenic.p.thresholds"){

    # Need to specify path
    p.Threshold = read.table(p.Threshold.file)[,3];
    temp = read.table(input.file, as.is=T);
    temp = temp[temp[,4]<0.95,]

    temp = cbind(temp, -qnorm(temp[,4]/2), abs(temp[,3]/qnorm(temp[,4]/2)));
    n.threshold = length(p.Threshold);

    for(k in 1:n.threshold){
       my.temp = temp[temp[,4]<=p.Threshold[k],];
       z0 =  -qnorm(p.Threshold[k]/2);
       lambda = z0*my.temp[,6];
       beta = (abs(my.temp[,3])-lambda)*sign(my.temp[,3]);
       write.table(cbind(my.temp[,1:2],beta),paste("LASSO.Score.",k,sep=""),row.names=F,col.names=F,quote=F);
    }
}


# plink.file is validation file.
Polygenic.Score <- function(plink.file, score.file.pre, n.threshold){

    command = paste("plink --bfile ", plink.file, " --noweb --make-bed --out temp1 > 1", sep="");
    system(command);

    for(k in 1:n.threshold){
        system(paste("wc -l ", score.file.pre, ".Score.", k, " > n.rows", sep=""));
        tt = read.table("n.rows")[1,1];
        if(tt>0){
            command = paste("plink --bfile temp1 --noweb --extract ", score.file.pre, ".Score.", k, " --make-bed --out temp2 > 1", sep="");
            system(command);
            command = paste("plink --bfile temp2 --noweb --score ", score.file.pre, ".Score.", k, " --out Single.", k, " > 1", sep="");
            print(k);
        }
    }

    system("rm *.log *.nopred *.nosex");
}
	
	sumx2fun = function(xx){
      sum(xx^2)
    }
	
	outputvector = function(tXX, testXX, trainOutput, pcutvec, tauvec, lambdavec){
  
	  outputMatrix = as.data.frame(matrix(,1,6))
	  names(outputMatrix) = c("naiveI","naiveJ","lassoI","AlassoI","lassoJ","AlassoJ")
	
      sumsquare = apply(tXX,2,sumx2fun)
  
      D = matrix(0,ncol(tXX),ncol(tXX))
      diag(D) = sumsquare
  
      inv_tXX = solve(t(tXX) %*% tXX )
      indBeta = matrix(trainOutput[,1],ncol(tXX),1)
      joinBeta = inv_tXX %*% D %*% indBeta 
   
      mediansd = median(trainOutput[,2])
  
      #native estimator at individual
      olsBetas = trainOutput[,1]
  
      Naive_preY_ind = testXX %*% olsBetas
    outputMatrix$naiveI[1] = cor(Naive_preY_ind,testY)^2
   
	
    #native estimator at joint
    Naive_preY_joint = testXX %*% joinBeta
    outputMatrix$naiveJ[1] = cor(Naive_preY_joint,testY)^2
    
   
    #lasso estimator at individual with bias
    preYlist = list()
    SSEvec = corvec = c()
    
	for(p1 in 1:length(pcutvec)){
      Zcut = qnorm(pcutvec[p1],lower.tail = FALSE)  
      lambda = Zcut*mediansd 
      olsBetas = trainOutput[,1]
      Betahat = sign(olsBetas)*(abs(olsBetas)-lambda)	
	  wzero = which(abs(olsBetas) <= lambda)
	  Betahat[wzero] = 0
	  corvec[p1] = cor((testXX %*% Betahat),testY)^2
    }
	max(corvec,na.rm=T)
    outputMatrix$lassoI[1] = max(corvec,na.rm=T)
  
    #adaptive lasso estimator at individual
    preYlist = list()
    SSEvec = c()
	corMatrix = matrix(,length(lambdavec),length(tauvec))
	for(t1 in 1:length(tauvec)){
      tau0 = tauvec[t1]
      
	  for(l1 in 1:length(lambdavec)){
	    lambda = lambdavec[l1]
		olsBetas = trainOutput[,1]
	    indlambdas = lambda/(abs(olsBetas))^tau0
	    
		Betahat = c()
	    
		for(jj in 1:length(olsBetas)){
	      if(abs(olsBetas[jj]) <= indlambdas[jj]){
	        Betahat[jj] = 0
	      }else{
	        Betahat[jj] = sign(olsBetas[jj])*(abs(olsBetas[jj])-indlambdas[jj])	
	      } 
	    }
	    corMatrix[l1,t1] = cor((testXX %*% Betahat),testY)^2
      }
	}  
	max(corMatrix,na.rm=T)
    outputMatrix$AlassoI[1] = max(corMatrix,na.rm=T)
  
    #lasso estimator at joint with bias
    preYlist = list()
    SSEvec = corvec = c()
	
    for(p1 in 1:length(pcutvec)){
      Zcut = qnorm(pcutvec[p1],lower.tail = FALSE)  
      lambda = Zcut*mediansd 
	  Betahat = sign(joinBeta)*(abs(joinBeta)-lambda)
      wzero = which(abs(joinBeta)<=lambda)
	  Betahat[wzero] = 0	
	  corvec[p1] = cor((testXX %*% Betahat),testY)^2
    }
	max(corvec,na.rm=T)
    outputMatrix$lassoJ[1] = max(corvec,na.rm=T)
  
    #adaptive lasso estimator at joint
 
    preYlist = list()
    SSEvec = c()
	corMatrix = matrix(,length(lambdavec),length(tauvec))
	
	for(t1 in 1:length(tauvec)){
      tau0 = tauvec[t1]
      
	  for(l1 in 1:length(lambdavec)){
    
	    lambda = lambdavec[l1]
	    indlambdas = lambda/(abs(joinBeta))^tau0
    
	    Betahat = c()
	    
		for(jj in 1:length(joinBeta)){
	      if(abs(joinBeta[jj])<=indlambdas[jj]){
	        Betahat[jj] = 0
	      }else{
	        Betahat[jj] = sign(joinBeta[jj])*(abs(joinBeta[jj])-indlambdas[jj])	
	      } 
	    }
	    corMatrix[l1,t1] = cor((testXX %*% Betahat),testY)^2
      }
	}  
	max(corMatrix,na.rm=T)
    outputMatrix$AlassoJ[1] = max(corMatrix,na.rm=T)
	return(outputMatrix)
   } 
  
	
	
	imputeMiss = function(xx){
		    wmiss = which(is.na(xx))
		    if(length(wmiss)==length(xx)){
			  meanval = 0
			}else{meanval = mean(xx,na.rm=T)}
			
			if(length(wmiss)>0){
			  xx[wmiss] = meanval
			}
		    return(xx)
    }
	
	imputeMissZero = function(xx){
      wmiss = which(is.na(xx))
	  if(length(wmiss)>0){
	    xx[wmiss] = 0
	  }
	  return(xx)
    }
	
	
	df_count = function(xx){
	  length(which(xx!=0))
	}
	
	Betascut_func = function(pcutvec, summaryBetasInfo){
	  betascut = c()
	  
	  for(p1 in 1:length(pcutvec)){
	    pcut = pcutvec[p1]
		wk = which(summaryBetasInfo$pvalue <= pcut)
		if(length(wk)==0){next}
		betascut = c(betascut,min(abs(summaryBetasInfo$betas[wk]),na.rm=T))
	  }
	  return(betascut)
	}
	
	
	Lasso_estimator = function(betascut, summaryBetas){
	  BetaMatrix = matrix(,length(betascut), length(summaryBetas))
	  
	  for(p1 in 1:length(betascut)){  
        lambda = betascut[p1] 
        Betahat = sign(summaryBetas)*(abs(summaryBetas)-lambda)	
	    wzero = which(abs(summaryBetas) <= lambda)
	    Betahat[wzero] = 0
		BetaMatrix[p1,] = Betahat 
	  }
	  return(BetaMatrix)
	}
	
    adaptiveLasso_estimator = function(betascut, summaryBetas, tauvec, lambdavec){
	  all_lambdavec = rep(lambdavec, times = length(tauvec))
	  all_tauvec = rep(tauvec, each = length(lambdavec))
	  
	  BetaMatrix = matrix(,length(all_lambdavec),length(summaryBetas))
	  
	  for(t1 in 1:length(all_lambdavec)){
        tau0 = tauvec[t1]
	    lambda = lambdavec[t1]
	    indlambdas = lambda/(abs(summaryBetas))^tau0
	    
		Betahat = c()	    
		for(jj in 1:length(summaryBetas)){
	      if(abs(olsBetas[jj]) <= indlambdas[jj]){
	        Betahat[jj] = 0
	      }else{
	        Betahat[jj] = sign(summaryBetas[jj])*(abs(summaryBetas[jj])-indlambdas[jj])	
	      } 
	    }
	    BetaMatrix[t1,] = Betahat
      }
	  return(BetaMatrix)
	}  
	
   
	
	Pruning_SNPs = function(LDdata, summaryBetasInfo, pruncut = 0.1){	
	  od = order(summaryBetasInfo$pvalue,decreasing=F)
	  summaryBetasInfo = summaryBetasInfo[od,,drop=F]
	  mat1 = match(summaryBetasInfo$SNP, colnames(LDdata) )
	  LDdata = LDdata[mat1,mat1]
      pruningLDdata = screeningLD(LDdata0 = LDdata, cor2cut = pruncut)
	  rm(pruningLDdata, LDdata)
	  return(summaryBetasInfo$SNP)
	}
	
	
	Pruning_SNPs_by_pcut = function(LDdata, betascut, summaryBetasInfo, pruncut = 0.1){
	
	  pruningSNPlist = list()	
	  for(p1 in 1:length(betascut)){
	    bcut = betascut[p1]
	    wkp = which(abs(summaryBetasInfo$betas) >= bcut)
	    if(length(wkp)==0){next}
	    pruning_summaryBetasInfo = summaryBetasInfo[wkp,,drop=F]
	    	
	    od = order(pruning_summaryBetasInfo$pvalue,decreasing=F)
	    pruning_summaryBetasInfo = pruning_summaryBetasInfo[od,,drop=F]
		keepsnps = pruning_summaryBetasInfo$SNP
		mat = match(keepsnps,colnames(LDdata))
		pruningLDdata = LDdata[mat,mat]
	
		if(length(wkp)>1){
		  pruningLDdata = screeningLD(LDdata0 = pruningLDdata, cor2cut = pruncut)
	    }
		pruningSNPlist[[p1]] = pruning_summaryBetasInfo$SNP
	  }
	  return(pruningSNPlist)
	}
	
	
	
	Pruning_SNPs_by_pcut_rawld = function(ld_outname, betascut, summaryBetasInfo, pruncut = 0.1){
	
	  pruningSNPlist = list()	
	  for(p1 in 1:length(betascut)){
	    bcut = betascut[p1]
	    wkp = which(abs(summaryBetasInfo$betas) >= bcut)
	    if(length(wkp)==0){next}
	    pruning_summaryBetasInfo = summaryBetasInfo[wkp,,drop=F]
	    	
	    od = order(pruning_summaryBetasInfo$pvalue,decreasing=F)
	    pruning_summaryBetasInfo = pruning_summaryBetasInfo[od,,drop=F]
		keepsnps = pruning_summaryBetasInfo$SNP
		if(p1==1){
	      lddata = read.table(file=ld_outname,header=T,as.is=T)
	      dosenames = keepsnps
		  LDmatrix = matrix(0,length(dosenames),length(dosenames))
	      colnames(LDmatrix) = rownames(LDmatrix) = dosenames
	      diag(LDmatrix) = 1
		
	      for(ii in 1:nrow(lddata)){
		    snp1 = lddata$SNP_A[ii]
		    snp2 = lddata$SNP_B[ii]
		    wd1 = which(dosenames==snp1)
		    wd2 = which(dosenames==snp2)
		    LDmatrix[wd1,wd2] = LDmatrix[wd2,wd1] = lddata$R[ii]
	      }	
		}
		LDdata = LDmatrix
		mat = match(keepsnps,colnames(LDdata))
		pruningLDdata = LDdata[mat,mat]
	
		if(length(wkp)>1){
		  pruningLDdata = screeningLD(LDdata0 = pruningLDdata, cor2cut = pruncut)
	    }
		pruningSNPlist[[p1]] = pruning_summaryBetasInfo$SNP
	  }
	  return(pruningSNPlist)
	}
	
	
	pruningLasso_estimator = function(betascut, summaryBetas, pruningSNP){
	  BetaMatrix = matrix(0,length(betascut), length(summaryBetas))
	  for(p2 in 1:length(betascut)){ 
          lambda = betascut[p2]   
		  Betahat  = summaryBetas
		  wN = which(!names(summaryBetas) %in% pruningSNP)
		  Betahat[wN] = 0
		  wzero = which(abs(Betahat)  <= lambda )  
          Betahat = sign(Betahat)*(abs(Betahat)-lambda)	  
	      Betahat[wzero] = 0
		  BetaMatrix[p2,] = Betahat  
	  } 
	  return(BetaMatrix)
	}
	
	pruningLasso_estimator_by_pcut = function(betascut, summaryBetas, pruningSNPlist){
	  BetaMatrix = matrix(0,length(betascut), length(summaryBetas))
	  for(p2 in 1:length(betascut)){
          pruningSNP = pruningSNPlist[[p2]]	  
          lambda = betascut[p2]   
		  Betahat  = summaryBetas
		  wN = which(!names(summaryBetas) %in% pruningSNP)
		  Betahat[wN] = 0
		  wzero = which(abs(Betahat)  <= lambda )  
          Betahat = sign(Betahat)*(abs(Betahat)-lambda)	  
	      Betahat[wzero] = 0
		  BetaMatrix[p2,] = Betahat  
	  } 
	  return(BetaMatrix)
	}
	
	
	
	
	mcounts = function(xx,mcut=-0.9){
      length(which(xx<mcut))
	}
	
	
	screeningNegativeLD_plink_ld = function(ld_outname, LDdata0, mcut = -0.7){
	  lddata = read.table(file=ld_outname,header=T,as.is=T)
	  
	  ncolnum = ncol(LDdata0)
	  jj = 1
	  while(jj <=ncolnum){
	    wrm = which(LDdata0[jj,]<mcut)
		if(length(wrm)>0){
		  LDdata0 = LDdata0[-wrm,-wrm]
	    }
        jj = jj + 1
        ncolnum = ncol(LDdata0)		  
	  }
	  return(LDdata0)
	}
	
	screeningNegativeLD = function(LDdata0, mcut = -0.7){
	  ncolnum = ncol(LDdata0)
	  jj = 1
	  while(jj <=ncolnum){
	    wrm = which(LDdata0[jj,]<mcut)
		if(length(wrm)>0){
		  LDdata0 = LDdata0[-wrm,-wrm]
	    }
        jj = jj + 1
        ncolnum = ncol(LDdata0)		  
	  }
	  return(LDdata0)
	}
	
	screeningLD = function(LDdata0, cor2cut = 0.1){
	  diag(LDdata0) = NA
	  ncolnum = ncol(LDdata0)
	  jj = 1
	  while(jj <=ncolnum){
	    wrm = which(LDdata0[jj,]^2>cor2cut)
		if(length(wrm)>0){
		  LDdata0 = LDdata0[-wrm,-wrm, drop=F]
	    }
        jj = jj + 1
        ncolnum = ncol(LDdata0)		  
	  }
	  diag(LDdata0) = 1
	  return(LDdata0)
	}
	
	
	
	iterSNP_Rversion1 = function(summaryBetas, lambda, tau, corXX, nter = 1000){
	
	  tmpbetas = tmpbetas0 = rep(0,ncol(corXX)) 
	  names(tmpbetas) = names(tmpbetas0) = colnames(corXX)
	  for(iter in 1:nter){
        for(jj in 1:length(tmpbetas)){
          tmp = (summaryBetas[jj] - sum(corXX[jj,-jj]*tmpbetas[-jj]))
		  
		  if(iter==1){
		    threshold = lambda/tau
		  }else{
		    threshold = lambda/(abs(tmpbetas0[jj]) + tau)
		  }
		  if(abs(tmp)>threshold){
		    tmp = sign(tmp)*(abs(tmp)-threshold)
		  }else{
		    tmp = 0
		  }
		  if(tmp>100){stop("tmp inf")}
	      tmpbetas[jj] = tmp
	    }
        if(max(abs(tmpbetas-tmpbetas0))<0.0001){break}
        tmpbetas0 = tmpbetas
      }
	  output = list()
	  output[[1]] = iter
	  output[[2]] = tmpbetas
	  return(output)
	}
	
	iterSNP_Rversion = function(summaryBetas, lambda, tau, corXX, nter = 1000){
	
	  tmpbetas = tmpbetas0 = rep(0,ncol(corXX)) 
	  names(tmpbetas) = names(tmpbetas0) = colnames(corXX)
	  for(iter in 1:nter){
        for(jj in 1:length(tmpbetas)){
          tmp = (summaryBetas[jj] - sum(corXX[jj,-jj]*tmpbetas0[-jj]))
		  
		  if(iter==1){
		    threshold = lambda/tau
		  }else{
		    threshold = lambda/(abs(tmpbetas0[jj]) + tau)
		  }
		  if(abs(tmp)>threshold){
		    tmp = sign(tmp)*(abs(tmp)-threshold)
		  }else{
		    tmp = 0
		  }
		  if(tmp>100){stop("tmp inf")}
	      tmpbetas[jj] = tmp
	    }
        if(max(abs(tmpbetas-tmpbetas0))<0.0001){break}
        tmpbetas0 = tmpbetas
      }
	  output = list()
	  output[[1]] = iter
	  output[[2]] = tmpbetas
	  return(output)
	}
	
	
	
