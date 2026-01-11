Create_sigma2K_allAlpha_List = function(K,sigma2Kvec,rhoMat){
    SharedPattern = permutations(n=2,r=K,v=c(0,1),repeats.allowed=T)
    numAlpha = nrow(SharedPattern)

    sigma2K_allAlpha_List = list()

    for(I in 2:numAlpha){
      Smat = matrix(0,K,K)
      tempindex = SharedPattern[I,]

      for(i in 1:ncol(SharedPattern)){
        for(j in i:ncol(SharedPattern)){
          if(i==j){
            if(tempindex[i]==1){
              Smat[i,j] = sigma2Kvec[i]
            }else{
              Smat[i,j] = 0
            }
          }else{
            if(tempindex[i]==1 & tempindex[j]==1){
              Smat[i,j] = Smat[j,i] = rhoMat[i,j]*sqrt(sigma2Kvec[i])*sqrt(sigma2Kvec[j])
            }else{
              Smat[i,j] = Smat[j,i] = 0
            }
          }
        }
      }
      sigma2K_allAlpha_List[[I]] = Smat
      rm(Smat)
    }
    return(sigma2K_allAlpha_List)
  }



create_names_sort = function(xx){
    xx = sort(xx)
    paste0(xx,collapse="-")
}
create_names = function(xx){
    paste0(xx,collapse="-")
}
check_order = function(xx){
  all(xx[,1]<xx[,2])
}


Gen_All_lddata_simpleTrans = function(outnames_vec, Allkeepsnps, rcut = 0.03){
    lddata1 = read.table(file=outnames_vec[1],header=T,as.is=T)
    lddata2 = read.table(file=outnames_vec[2],header=T,as.is=T)

    missR = c(NA,"NaN")
    wrm = which(lddata1$R %in% missR)
    if(length(wrm)>0){
      lddata1 = lddata1[-wrm,]
    }
    wkeep = which( abs(lddata1$R) >= rcut)
    lddata1 = lddata1[wkeep,]
    
    wkeep = which(lddata1$SNP_A  %in% Allkeepsnps & lddata1$SNP_B %in% Allkeepsnps)
    lddata1 = lddata1[wkeep,]

    wrm = which(lddata2$R %in% missR)
    if(length(wrm)>0){
      lddata2 = lddata2[-wrm,]
    }
    wkeep = which( abs(lddata2$R) >= rcut)
    lddata2 = lddata2[wkeep,]
    wkeep = which(lddata2$SNP_A  %in% Allkeepsnps & lddata2$SNP_B %in% Allkeepsnps)
    lddata2 = lddata2[wkeep,]

    if(!check_order(lddata1[,c(2,5)])){print("error name1")}
    if(!check_order(lddata2[,c(2,5)])){print("error name2")}

    name1 = apply(lddata1[,c(2,5)],1,create_names)
    name2 = apply(lddata2[,c(2,5)],1,create_names)

    lddata1$name = name1
    lddata2$name = name2

    All_lddata= merge(lddata1,lddata2,by="name",all=T)
    wkX = which(!is.na(All_lddata$BP_A.x))
    wkY = which(is.na(All_lddata$BP_A.x) & !is.na(All_lddata$BP_A.y))

    Xnames = c("CHR_A.x","BP_A.x","SNP_A.x",
    "CHR_B.x", "BP_B.x","SNP_B.x", "R.x","R.y")
    Ynames = c("CHR_A.y","BP_A.y","SNP_A.y",
    "CHR_B.y", "BP_B.y","SNP_B.y", "R.x","R.y")

    Finalnames = c("CHR_A","BP_A","SNP_A",
    "CHR_B", "BP_B","SNP_B", "R.1","R.2")

    if(length(wkX)!=0 & length(wkY)!=0){
      part1 = All_lddata[wkX,Xnames]
      part2 = All_lddata[wkY,Ynames]
      names(part1) = Finalnames
      names(part2) = Finalnames
      All_lddata_final = rbind(part1,part2)
      names(All_lddata_final) = Finalnames
    }else{
      All_lddata_final = All_lddata[wkX,Xnames]
      names(All_lddata_final) = Finalnames
    }
    rm(All_lddata)
    gc()
    for(ii in 1:ncol(All_lddata_final)){
       wNA = which(is.na(All_lddata_final[,ii]))
       if(length(wNA)>0){
         All_lddata_final[wNA,ii] = 0
       }
    }
    names = apply(All_lddata_final[,c(2,5)],1,create_names)
    if(!all(lddata1$name %in% names)){print("error in lddata1")}
    if(!all(lddata2$name %in% names)){print("error in lddata2")}
    return(All_lddata_final)
}



Gen_All_lddata_simpleTrans3 <- function(outnames_vec, Allkeepsnps, rcut = 0.03) {
    
    # ---- Load 3 LD files ----
    ld1 <- read.table(outnames_vec[1], header = TRUE, as.is = TRUE)
    ld2 <- read.table(outnames_vec[2], header = TRUE, as.is = TRUE)
    ld3 <- read.table(outnames_vec[3], header = TRUE, as.is = TRUE)

    # ---- Clean function ----
    clean_ld <- function(ld) {
        missR <- c(NA, "NaN")
        ld <- ld[!(ld$R %in% missR), ]
        ld <- ld[abs(ld$R) >= rcut, ]
        ld <- ld[ld$SNP_A %in% Allkeepsnps & ld$SNP_B %in% Allkeepsnps, ]
        ld
    }

    ld1 <- clean_ld(ld1)
    ld2 <- clean_ld(ld2)
    ld3 <- clean_ld(ld3)

    # ---- Name checks ----
    if (!check_order(ld1[, c(2, 5)])) print("error name1")
    if (!check_order(ld2[, c(2, 5)])) print("error name2")
    if (!check_order(ld3[, c(2, 5)])) print("error name3")

    # ---- Create merge key ----
    name1 <- apply(ld1[, c(2, 5)], 1, create_names)
    name2 <- apply(ld2[, c(2, 5)], 1, create_names)
    name3 <- apply(ld3[, c(2, 5)], 1, create_names)

    ld1$name <- name1
    ld2$name <- name2
    ld3$name <- name3

    # ---- Rename R ? R.1, R.2, R.3 ----
    names(ld1)[names(ld1) == "R"] <- "R.1"
    names(ld2)[names(ld2) == "R"] <- "R.2"
    names(ld3)[names(ld3) == "R"] <- "R.3"

    # ---- Merge 3 LD datasets by "name" ----
    All_lddata <- merge(ld1, ld2, by = "name", all = TRUE)
    All_lddata <- merge(All_lddata, ld3, by = "name", all = TRUE)

    # ---- Identify rows from each dataset ----
    wk1 <- which(!is.na(All_lddata$BP_A.x))                                     # pop1
    wk2 <- which(is.na(All_lddata$BP_A.x) & !is.na(All_lddata$BP_A.y))          # pop2
    wk3 <- which(is.na(All_lddata$BP_A.x) & is.na(All_lddata$BP_A.y) & 
                 !is.na(All_lddata$BP_A))                                       # pop3

    # ---- Column names for each part ----
    part1_cols <- c("CHR_A.x","BP_A.x","SNP_A.x",
                    "CHR_B.x","BP_B.x","SNP_B.x","R.1","R.2","R.3")

    part2_cols <- c("CHR_A.y","BP_A.y","SNP_A.y",
                    "CHR_B.y","BP_B.y","SNP_B.y","R.1","R.2","R.3")

    part3_cols <- c("CHR_A","BP_A","SNP_A",
                    "CHR_B","BP_B","SNP_B","R.1","R.2","R.3")

    Finalnames <- c("CHR_A","BP_A","SNP_A",
                    "CHR_B","BP_B","SNP_B",
                    "R.1","R.2","R.3")

    # ---- Build final dataset ----
    part1 <- All_lddata[wk1, part1_cols]
    part2 <- All_lddata[wk2, part2_cols]
    part3 <- All_lddata[wk3, part3_cols]

    names(part1) <- Finalnames
    names(part2) <- Finalnames
    names(part3) <- Finalnames

    All_lddata_final <- rbind(part1, part2, part3)
    names(All_lddata_final) <- Finalnames

    # ---- Replace missing with zero ----
    for (ii in 1:ncol(All_lddata_final)) {
        wNA <- which(is.na(All_lddata_final[, ii]))
        if (length(wNA) > 0) {
            All_lddata_final[wNA, ii] <- 0
        }
    }

    # ---- Final consistency checks ----
    names_check <- apply(All_lddata_final[, c(2, 5)], 1, create_names)
    if (!all(ld1$name %in% names_check)) print("error in lddata1")
    if (!all(ld2$name %in% names_check)) print("error in lddata2")
    if (!all(ld3$name %in% names_check)) print("error in lddata3")

    return(All_lddata_final)
}



Gen_All_lddata_simpleTrans3traits <- function(outname, Allkeepsnps, rcut = 0.03) {

    # ---- Load 1-pop LD ----
    ld <- read.table(file = outname, header = TRUE, as.is = TRUE)

    # ---- Clean missing / invalid R ----
    missR <- c(NA, "NaN")
    ld <- ld[!(ld$R %in% missR), ]

    # ---- Apply |R| cutoff ----
    ld <- ld[abs(ld$R) >= rcut, ]

    # ---- Keep only SNP pairs in Allkeepsnps ----
    ld <- ld[ld$SNP_A %in% Allkeepsnps & ld$SNP_B %in% Allkeepsnps, ]

    # ---- Ordering check (same structure as 2-pop function) ----
    if (!check_order(ld[, c("BP_A", "BP_B")])) {
        print("warning: LD ordering check failed")
    }

    # ---- Create R.1, R.2, R.3 (all identical for 1-pop case) ----
    ld$R.1 <- ld$R
    ld$R.2 <- ld$R
    ld$R.3 <- ld$R

    # ---- Final output columns (same as 3-pop version) ----
    final <- ld[, c("CHR_A", "BP_A", "SNP_A",
                    "CHR_B", "BP_B", "SNP_B",
                    "R.1", "R.2", "R.3")]

    # ---- Replace NA with 0 ----
    final[is.na(final)] <- 0

    return(final)
}






GenPreR2_Chrs_saveScore = function(AllBetaMatrix_list, savef2_list, refAllelefile_list, output_sub_folder_predi, PV, popuseY, dirSimuData, sim, savename, Scoreoutput,savecore){

    if(PV==""){
      plinkver = 1
    }
    if(PV==2){
      plinkver = 2
    }


    ytestvalfile = paste0(dirSimuData,"test_validation.pheno_", popuseY, "_sim", sim)
    tvfile =  paste0(dirSimuData,"Test_ValiDid_", popuseY, "_sim", sim, ".txt")
    testfile = paste0(dirSimuData,"TestDid_", popuseY, "_sim", sim, ".txt")
    valfile = paste0(dirSimuData,"ValiDid_", popuseY, "_sim", sim, ".txt")


    
    for(ll in 1:length(AllBetaMatrix_list)){
      #-- get phenotype
      phenoOri = read.table(file=ytestvalfile,as.is=T)
      rownames(phenoOri) = phenoOri[,1]
      
      AllBetaMatrix = AllBetaMatrix_list[[ll]]
      savef2 = savef2_list[[ll]]
      refAllelefile = refAllelefile_list[[ll]]
      nonzero_test = apply(AllBetaMatrix,2,nonzero)
      wcol = which(nonzero_test>0)
      AllBetaMatrix = AllBetaMatrix[,wcol,drop=F]

      OriBeta = t(AllBetaMatrix)
      rm(AllBetaMatrix)
      AllSNPs = rownames(OriBeta)
      if(length(AllSNPs)==0){
        PreR2 = NULL
        Finish = 1
        save(PreR2,Finish, file=savef2)
      }else{
        savefile = paste0(output_sub_folder_predi, "allchrs.AllSNPs",ll,".txt")
        write.table(AllSNPs,file=savefile,row.names=F,col.names=F,quote=F);

        refAllele = read.table(file=refAllelefile, as.is=T, header=F)
        mat = match(AllSNPs, refAllele[,1])
        refallele = refAllele[mat,]
        dataname = paste0("/lustre06/project/6005709/yatah3/simulation/SimuGenotype/",popuseY,"AllChrs_bedformat")
        command = paste0("plink",PV," --bfile ",dataname," --keep ",tvfile," --allow-no-sex --extract ", savefile,  " --threads 1 --make-bed --out ",paste0(output_sub_folder_predi,"temp",ll," > 1"));
        system(command);

        #if(file.exists(savef2)){
          #load(file=savef2)
        #}else{
          PreR2 = matrix(NA, 2, ncol(OriBeta))
          numbetasvec = rep(NA,ncol(OriBeta))
          Finish = 0
        #}
        testid = read.table(file=testfile,as.is=T)
        valid = read.table(file=valfile,as.is=T)

        testk = 1:ncol(OriBeta)
        for(iik in 1:length(testk)){
          k = testk[iik]
          if(k==ncol(OriBeta)){Finish=1}
          if(savename!="TRUE"){
            if(length(testk) > 25){
              if(iik %in% seq(1,length(testk),by=5)){
                print(apply(PreR2,1,max,na.rm=T))
                print(table(numbetasvec))
                save(PreR2, numbetasvec, Finish,file=savef2)
              }
            }
          }
          if(!is.na(PreR2[2,k])){next}
          oriBetas = OriBeta[,k]
          if(any(is.na(oriBetas))){next}
          wnonzero = which(oriBetas!=0)
          numbetasvec[k] = length(wnonzero)

          if(length(wnonzero)==0){next}
          scoreinfo = cbind(refallele,oriBetas)[wnonzero,,drop=F]
          savefile2 = paste0(output_sub_folder_predi, "Score.",ll,"_",k)
          write.table(scoreinfo,file=savefile2,row.names=F,col.names=F,quote=F);
          rm(scoreinfo)
          outfile = paste0(Scoreoutput,savecore)
          command = paste0("plink",PV,"  --allow-no-sex --bfile ", paste0(output_sub_folder_predi,"temp",ll), " --score ", savefile2, " --threads 1 --out ",outfile, "Single.",ll,"_", k,".profile > 1", sep="");
          system(command);

          if(plinkver==1){
            preYplink0 = read.table(file=paste0(outfile,"Single.",ll,"_", k,".profile.profile"), as.is=T,header=T)
            file.remove(paste0(outfile,"Single.",ll,"_", k,".profile.profile"))
          }
          if(plinkver==2){
            preYplink0 = read.table(file=paste0(outfile,"Single.",ll,"_", k,".profile.sscore"), as.is=T)
            file.remove(paste0(outfile,"Single.",ll,"_", k,".profile.sscore"))
            names(preYplink0)[6] = "SCORE"
          }

          wk = which(rownames(phenoOri) %in% preYplink0[,1])
          pheno = phenoOri[wk,]
          mat = match(rownames(pheno),preYplink0[,1])
          preYplink0 = preYplink0[mat,]

          mat1 = match(testid[,1],preYplink0[,1])
          mat2 = match(valid[,1],preYplink0[,1])

          PreR2[1,k] = cor(pheno[mat1,phenoIndex],preYplink0$SCORE[mat1])^2
          PreR2[2,k] = cor(pheno[mat2,phenoIndex],preYplink0$SCORE[mat2])^2
          rm(pheno,preYplink0)
        }
        rm(OriBeta)
        gc()
        save(PreR2, numbetasvec, Finish, file=savef2)
        
        rm(PreR2, numbetasvec, Finish)
        
      }
      rm(AllBetaMatrix)
    }
}



TransPEN_MixLog <- function(
  Zscale = 1,
  usedtrait,                 # <-- matches your call
  TrainingNsam,
  dirSimuData,               # <-- matches your call (if you use it later)
  warmStart = 1,
  penalty = "mixLOG",
  NumIter = 1000,
  outnames_vec,
  mainindex = 1,
  RupperVal = NULL,
  gwasfilename,
  mafnames_vec = NULL,
  singleStart = 0,
  output_sub_folder,         # <-- matches your call
  AbsTauvec = NULL,
  sim,
  savefile0,
  savef0,
  chr = NULL,
  WeightN,
  NumL = 10,
  subNumL = 10,
  low0_single = 0,
  top0_single = 0,
  rcut = 0.3                 # use same default here as helper
){
  
    if(penalty == "RealmixLOG"){
      tuningMatrix0 = Gen_tuningMatrix_RealmixLOG_trans(top0_single, low0_single, NumL = NumL, subNumL = subNumL)
    }
    if(penalty == "mixLOG"){
      tuningMatrix0 = Gen_tuningMatrix_mixLOG_trans(savefile0, savef0, perR = 0.8, NumL = 50, subNumL = 50)
    }
    if(penalty == "Log"){
      tuningMatrix0 = Gen_tuningMatrix_LOG_trans(savefile0,savef0, perR = 0.6, NumL = 100)
    }

    penalty1 = "mixLOG"
    usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
    usedtraitIndex = as.numeric(usedtraitsvec)
    Nvec = TrainingNsam[usedtraitIndex]
    GWASbetafile = gwasfilename

    summaryZoutput = Get_summaryZ(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr)
    summaryZ = summaryZoutput[[1]]
    Allkeepsnps = summaryZoutput[[2]]
    plinkLDList = Gen_plinkLDList(outnames_vec, Allkeepsnps,rcut=0.03)
    ordersequse=1
    for(tt in 1:ncol(AbsTauvec)){
      tauuse = AbsTauvec[sim,tt]
      savefile1 = paste0(output_sub_folder, penalty,"chr",chr,"usedtrait",usedtrait,"warmStart",warmStart,"sim",sim,"Zscale",Zscale,"tauuse",tauuse,"singleStart",singleStart,".RData")
     # if(file.exists(savefile1)){next}
      tuningMatrix = tuningMatrix0
      tuningMatrix[,2] = tuningMatrix[,2]*tauuse
      tau = tauuse
      output = transPEN(summaryZ=summaryZ, plinkLD = plinkLDList, NumIter=NumIter,
        breaking = 0,  Nvec = Nvec, penalty=penalty1,
        tuningMatrix= tuningMatrix, dfMax = nrow(summaryZ), warmStart = warmStart, outputAll=1, orderseq = ordersequse, tau = tau, singleStart = singleStart, weightN = WeightN)

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
      rm(output,saveoutput,tuningMatrix)
      gc()
   }
   rm(plinkLDList)
   gc()

}

TransPEN_MixLog3 = function(Zscale = 1,
  usedtrait,                 # <-- matches your call
  TrainingNsam,
  dirSimuData,               # <-- matches your call (if you use it later)
  warmStart = 1,
  penalty = "mixLOG",
  NumIter = 1000,
  outnames_vec,
  mainindex = 1,
  RupperVal = NULL,
  gwasfilename,
  mafnames_vec = NULL,
  singleStart = 0,
  output_sub_folder,         # <-- matches your call
  AbsTauvec = NULL,
  sim,
  savefile0,
  savef0,
  chr = NULL,
  WeightN = 0,
  NumL = 10,
  subNumL = 10,
  low0_single = 0,
   top0_single1 = 0,
    top0_single2 = 0,
    top0_single3 = 0, 
    rcut=0.03){
  
    if(penalty == "RealmixLOG"){
      tuningMatrix0 = Gen_tuningMatrix_RealmixLOG_trans3(top0_single1,top0_single2, top0_single3 , low0_single, NumL = NumL, subNumL = subNumL)
    }
    if(penalty == "mixLOG"){
      tuningMatrix0 = Gen_tuningMatrix_mixLOG_trans(savefile0, savef0, perR = 0.8, NumL = 50, subNumL = 50)
    }
    if(penalty == "Log"){
      tuningMatrix0 = Gen_tuningMatrix_LOG_trans(savefile0,savef0, perR = 0.6, NumL = 100)
    }

    penalty1 = "mixLOG"
    usedtraitIndex = c(1,2,3)
    Nvec = TrainingNsam

    summaryZoutput = Get_summary_3(gwasfilename, Nvec, chr)
    summaryZ = summaryZoutput[[1]]
    Allkeepsnps = summaryZoutput[[2]]
    
    
    length(Allkeepsnps)

    
    plinkLDList = Gen_plinkLDList(outnames_vec, Allkeepsnps,rcut=0.03)

    ordersequse=1
    
    for(tt in 1:ncol(AbsTauvec)){
      tauuse = AbsTauvec[sim,tt]



 savefile1 = paste0(output_sub_folder, penalty,"chr",chr,"_3pop","warmStart",warmStart,"sim",sim,"Zscale",Zscale,"tauuse",tauuse,"singleStart",singleStart,".RData")
 
      #if(file.exists(savefile1)){next}
      tuningMatrix = tuningMatrix0
      tuningMatrix[,2] = tuningMatrix[,2]*tauuse
      tau = tauuse
      
      output = transPEN(summaryZ=summaryZ, plinkLD = plinkLDList, NumIter=NumIter,
        breaking = 0,  Nvec = Nvec, penalty=penalty1,
        tuningMatrix= tuningMatrix, dfMax = nrow(summaryZ), warmStart = warmStart, outputAll=1, orderseq = ordersequse, tau = tau, singleStart = singleStart, weightN = WeightN)

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
      rm(output,saveoutput,tuningMatrix)
      gc()
   }
   rm(plinkLDList)
   gc()

}


TransPEN_MixLog3traits = function(Zscale = 1,
  usedtrait,                 # <-- matches your call
  TrainingNsam,
  dirSimuData,               # <-- matches your call (if you use it later)
  warmStart = 1,
  penalty = "mixLOG",
  NumIter = 1000,
  outnames_vec,
  mainindex = 1,
  RupperVal = NULL,
  gwasfilename,
  mafnames_vec = NULL,
  singleStart = 0,
  output_sub_folder,         # <-- matches your call
  AbsTauvec = NULL,
  sim,
  savefile0,
  savef0,
  chr = NULL,
  WeightN = 0,
  NumL = 10,
  subNumL = 10,
  low0_single = 0,
   top0_single1 = 0,
    top0_single2 = 0,
    top0_single3 = 0, 
    rcut=0.03){
  
    if(penalty == "RealmixLOG"){
      tuningMatrix0 = Gen_tuningMatrix_RealmixLOG_trans3(top0_single1,top0_single2, top0_single3 , low0_single, NumL = NumL, subNumL = subNumL)
    }
    if(penalty == "mixLOG"){
      tuningMatrix0 = Gen_tuningMatrix_mixLOG_trans(savefile0, savef0, perR = 0.8, NumL = 50, subNumL = 50)
    }
    if(penalty == "Log"){
      tuningMatrix0 = Gen_tuningMatrix_LOG_trans(savefile0,savef0, perR = 0.6, NumL = 100)
    }

    penalty1 = "mixLOG"
    usedtraitIndex = c(1,2,3)
    Nvec = rep(TrainingNsam,3)

    summaryZoutput = Get_summary_3(gwasfilename, Nvec, chr)
    summaryZ = summaryZoutput[[1]]
    Allkeepsnps = summaryZoutput[[2]]
    
    
    length(Allkeepsnps)


  plinkLDList <- Gen_plinkLDList3traits(outnames_vec, Allkeepsnps,rcut=0.03)


    
 
    ordersequse=1
    
    for(tt in 1:ncol(AbsTauvec)){
      tauuse = AbsTauvec[sim,tt]



 savefile1 = paste0(output_sub_folder, penalty,"chr",chr,"_3traits","warmStart",warmStart,"sim",sim,"Zscale",Zscale,"tauuse",tauuse,"singleStart",singleStart,".RData")
 
      #if(file.exists(savefile1)){next}
      tuningMatrix = tuningMatrix0
      tuningMatrix[,2] = tuningMatrix[,2]*tauuse
      tau = tauuse
      
      output = transPEN(summaryZ=summaryZ, plinkLD = plinkLDList, NumIter=NumIter,
        breaking = 0,  Nvec = Nvec, penalty=penalty1,
        tuningMatrix= tuningMatrix, dfMax = nrow(summaryZ), warmStart = warmStart, outputAll=1, orderseq = ordersequse, tau = tau, singleStart = singleStart, weightN = WeightN)

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
      rm(output,saveoutput,tuningMatrix)
      gc()
   }
   rm(plinkLDList)
   gc()

}





    
simpleTrans = function(Zscale = 1, popvec, usedtrait, TrainingNsam, dirSimuData, savefile1, low0_single, top0_single, totaltuning_single, NumIter= 1000, outnames_vec, gwasfilename = "GWASbetaStandard.txt",mafnames_vec, mainindex = 1, chr = NULL){
  usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
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
  All_lddata = Gen_All_lddata_simpleTrans(outnames_vec, Allkeepsnps)
  if(Zscale==1){
    summaryZ = matrix(Meta[,3],nrow(Meta),1)
    rownames(summaryZ) = GWASbeta$SNP
    summaryZ = check_summaryZ(summaryZ)
    NamesSNP = rownames(summaryZ)
    summaryBetas = SSvec = SDvec = NULL
  }else{
    summaryZoutput = Get_summaryBetas(GWASbetafile, mainindex, Nvec, chr)
    SDvec = summaryZoutput[[5]]
    mat = match(rownames(SDvec),Allkeepsnps)
    SDvec = SDvec[mat,]
    summaryZ = NULL
    summaryBetas = matrix(Meta[,1],nrow(Meta),1)
    rownames(summaryBetas) = GWASbeta$SNP
    SSvec = matrix(0,nrow(summaryBetas),ncol(summaryBetas))
    rownames(SSvec) = rownames(summaryBetas)
    NamesSNP = rownames(summaryBetas)

    SSvec0 = matrix(0,nrow(summaryBetas),2)
    rownames(SSvec0) = rownames(summaryBetas)

    for(ii in 1:length(mafnames_vec)){
      MAFdata0 = read.table(file=mafnames_vec[ii],as.is=T,header=T)
      mat = match(rownames(SSvec),MAFdata0$SNP)
      SSvec0[,ii] = 2*MAFdata0$MAF[mat]*(1-MAFdata0$MAF[mat])*Nvec[ii]
      rm(MAFdata0)
    }
    SSvec[,1] = apply(SSvec0,1,sum)
    rm(SSvec0)
  }
  if(is.null(top0_single)){
    print(paste0("maxZ=",maxZ))
    top0_single = maxZ
  }
  if(Zscale==0){
    topv = top0_single*median(Meta[,2],na.rm=T)*median(c(SSvec))
    lowv = low0_single*min(Meta[,2],na.rm=T)*min(SSvec)
  }else{
    topv = top0_single
    lowv = low0_single
  }

  lambdavec = seq(topv,lowv,len=totaltuning_single)
  tuningMatrix = cbind(lambdavec,1,0,1)
  rm(GWASbeta)

  output = simpletransPEN(summaryZ=summaryZ, summaryBetas = summaryBetas, SSvec = SSvec, SDvec = SDvec, plinkLD= All_lddata, NumIter=NumIter, breaking = 0,  Nvec = Nvec,Zscale = Zscale,
        tuningMatrix= tuningMatrix, dfMax = length(NamesSNP), outputAll=1)
  names(output)
  saveoutput = list()
  saveoutput[[1]] = output$Numitervec
  saveoutput[[2]] = ncol(output$BetaMatrix)
  saveoutput[[3]] = nrow(output$BetaMatrix)
  saveoutput[[4]] = which(output$BetaMatrix!=0)
  wpos = which(output$BetaMatrix!=0)
  saveoutput[[5]] = output$BetaMatrix[wpos]
  saveoutput[[6]] = NamesSNP
  saveoutput[[7]] = output[[3]]

  rm(output,wpos,summaryBetas,summaryZ, All_lddata)
  save(saveoutput,file=savefile1)
  rm(saveoutput)
  gc()
}




singleTrait = function(Zscale = 1, usedtrait, TrainingNsam, dirPlinkFormat1, savefile1, low0_single, top0_single, totaltuning_single, NumIter= 1000, ld_outname, warmStart = 1, mainindex = 1, RupperVal = NULL, gwasfilename, notedSkip = 0, mafnames_vec = NULL , savefile0 = NULL, totaltuning1 = 0, totaltuning2 = 0, singleStart = 0, penalty= "mixLOG" , causalSNPs=NULL, chr = chr, fixWeight = 0, NumL= 20, subNumL = 20, timesT = 1, rcut = 0.6, alpha = 1){

  penalty1 = penalty
  
  if(penalty == "RealmixLOG"){
    penalty1 = "mixLOG"
  }
  
  Nvec =c(10000,10000)

   summaryZoutput = Get_summaryZ(gwasfilename, mainindex, Nvec, usedtraitsvec,chr)
    summaryZ = summaryZoutput[[1]]
    Allkeepsnps = summaryZoutput[[2]]
    NamesSNP = rownames(summaryZ)
      Nq = length(Nvec)
      
    summaryBetas = matrix(0, nrow(summaryZ), Nq)
    SDvec = matrix(0, nrow(summaryZ), Nq)

    for(ii in 1:Nq){
      summaryBetas[,ii] = summaryZ[,ii]/sqrt(Nvec[ii])
      SDvec[,ii] = 1/sqrt(Nvec[ii])
    }
    rownames(summaryBetas) = rownames(summaryZ)
    SSvec = matrix(0,1,1)
    
  plinkLDList = Gen_plinkLDList(ld_outname, Allkeepsnps)


  if(is.null(top0_single)){
    #print(paste0("maxZ=",maxZ))
    top0_single = maxZ
  }
  if(Zscale==0){
    topv = top0_single*median(c(unlist(SDvec)))*median(c(SSvec))
    lowv = low0_single*min(c(unlist(SDvec)))*min(SSvec)
  }else{
    topv = top0_single
    lowv = low0_single
  }

  addweights = 0
  
  if(fixWeight==0){
    usefixWeight = 0
  }
  
    Lambda_limit = c(0.5,0.9) 
   subtuning = 10
    Lenlam = 10
    taufactor = c(1/25, 1, 10)
   
 medianval = median(apply(abs(summaryBetas),1,sum),na.rm=T)
 tauvec = sort(medianval*taufactor)
    Lambda_limit = quantile(abs(summaryZ[,1]), Lambda_limit)

tuningMatrix = Tuning_setup_group_only(tauvec, subtuning, Lambda_limit, Lenlam, medianval, Nq)



  output = gsPEN(summaryZ=summaryZ, summaryBetas = summaryBetas, SSvec = SSvec, SDvec = SDvec, ,  Nvec = Nvec , plinkLD= plinkLDList, NumIter=NumIter,penalty = penalty1,
        breaking = 0, RupperVal = RupperVal, Zscale = Zscale,
        tuningMatrix= tuningMatrix, dfMax = length(NamesSNP), outputAll=1, warmStart = warmStart, notedSkip = 0,  addweights =  addweights, singleStart = 1, fixWeight = usefixWeight)
  names(output)
  saveoutput = list()
  saveoutput[[1]] = output$Numitervec
  saveoutput[[2]] = ncol(output$BetaMatrix)
  saveoutput[[3]] = nrow(output$BetaMatrix)
  saveoutput[[4]] = which(output$BetaMatrix!=0)
  wpos = which(output$BetaMatrix!=0)
  saveoutput[[5]] = output$BetaMatrix[wpos]
  saveoutput[[6]] = NamesSNP
  saveoutput[[7]] = output[[3]]
  saveoutput[[10]] = apply(output$BetaMatrix,1,nonzero)

  if(notedSkip==1){
    saveoutput[[8]] = which(output$NotedSkipMatrix!=0)
    wpos2 = which(output$NotedSkipMatrix!=0)
    saveoutput[[9]] = output$NotedSkipMatrix[wpos2]
  }

  rm(output,wpos,summaryZ)
  save(saveoutput,file=savefile1)
  rm(saveoutput,plinkLDList)
  gc()
}


Gen_Prediction_CV_Chrs = function(savefile1List,savef2, groupnum, output_sub_folder_predi, refAllelefile, popuse, phenoIndex,gindex, plinkver, cindex, popuseY, dirSimuDataSet, sim, savename,penalty){
    if(is.null(popuseY)){
      popuseY = popuse
    }
    if(plinkver==1){
      PV = ""
    }
    if(plinkver==2){
      PV = "2"
    }

    for(chr in 1:22){
      savefile1 = savefile1List[[chr]]
      load(file=savefile1)
      iterindex = saveoutput[[cindex]]
      ncol_BetaMatrix = saveoutput[[2]]
      nrow_BetaMatrix = saveoutput[[3]]
      wpos = saveoutput[[4]]
      vals = saveoutput[[5]]
      BetaMatrix = matrix(0,nrow_BetaMatrix,ncol_BetaMatrix)
      BetaMatrix[wpos] = vals
      rm(wpos,vals)
      P = length(saveoutput[[6]])
      gindex1vec =seq(1,(groupnum*P),by=P)
      gindex2vec =seq(P,(groupnum*P),by=P)
      gindex1 = gindex1vec[gindex]
      gindex2 = gindex2vec[gindex]
      BetaMatrix = BetaMatrix[, gindex1:gindex2, drop=F]
      colnames(BetaMatrix) = saveoutput[[6]]

      if(chr==1){
        Alliter = iterindex
        AllBetaMatrix = BetaMatrix
      }else{
        Alliter = rbind(Alliter, iterindex)
        AllBetaMatrix = cbind(AllBetaMatrix, BetaMatrix)
      }
      rm(saveoutput,BetaMatrix)
    }
    gc()
    test = apply(Alliter,2,min)
    for(jj in 1:length(test)){
      if(test[jj] <= 0){
        AllBetaMatrix[jj,] = 0
      }
    }
    #print(table(apply(AllBetaMatrix,1,nonzero)))
    GenPreR2_Chrs(AllBetaMatrix, savef2, refAllelefile, output_sub_folder_predi, PV, popuseY, dirSimuDataSet, sim, savename)
   
}

GenPreR2_Chrs = function(AllBetaMatrix, savef2, refAllelefile, output_sub_folder_predi, PV, popuseY, dirSimuDataSet, sim, savename){

    if(PV=="1"){
      plinkver = 1
    }
    if(PV==2){
      plinkver = 2
    }
    
    ytestvalfile = paste0(dirSimuData,"test_validation.pheno_",popuseY,"sim",sim)
    tvfile =  paste0(dirSimuData,"Test_ValiDid",popuseY,"sim",sim,".txt")
    testfile = paste0(dirSimuData,"TestDid",popuseY,"sim",sim,".txt")
    valfile = paste0(dirSimuData,"ValiDid",popuseY,"sim",sim,".txt")


    #-- get phenotype
    phenoOri = read.table(file=ytestvalfile,as.is=T)
    rownames(phenoOri) = phenoOri[,1]
    nonzero_test = apply(AllBetaMatrix,2,nonzero)
    wcol = which(nonzero_test>0)
    AllBetaMatrix = AllBetaMatrix[,wcol,drop=F]

    OriBeta = t(AllBetaMatrix)
    rm(AllBetaMatrix)
    AllSNPs = rownames(OriBeta)
    if(length(AllSNPs)==0){
      PreR2 = NULL
      Finish = 1
      save(PreR2,Finish, file=savef2)
    }else{

    savefile = paste0(output_sub_folder_predi, "allchrs.AllSNPs.txt")
    write.table(AllSNPs,file=savefile,row.names=F,col.names=F,quote=F);

    refAllele = read.table(file=refAllelefile, as.is=T, header=F)
    mat = match(AllSNPs, refAllele[,1])


    refallele = refAllele[mat,]
    dataname = paste0(dirSimuDataSet,popuseY,"AllChrs_bedformat")
    command = paste0("plink",PV," --bfile ",dataname," --keep ",tvfile," --allow-no-sex --extract ", savefile,  " --threads 1 --make-bed --out ",paste0(output_sub_folder_predi,"temp > 1"));
    system(command);

    if(file.exists(savef2)){
      load(file=savef2)
    }else{
      PreR2 = matrix(NA, 2, ncol(OriBeta))
      numbetasvec = rep(NA,ncol(OriBeta))
      Finish = 0
    }
    testid = read.table(file=testfile,as.is=T)
    valid = read.table(file=valfile,as.is=T)
    
    
    
    #if(savename=="single"){
    #  int = 10
    #  testk = c(seq(1,ncol(OriBeta),by = int),ncol(OriBeta))
    #}else{
    #  testk = 1:ncol(OriBeta)
    #}
    testk = 1:ncol(OriBeta)
    for(iik in 1:length(testk)){
      k = testk[iik]
      if(k==ncol(OriBeta)){Finish=1}
      if(savename!="TRUE"){
      if(iik %in% seq(25,length(testk),by=25)){
         print(apply(PreR2,1,max,na.rm=T))
         save(PreR2, numbetasvec, Finish,file=savef2)
      }
      }
      if(!is.na(PreR2[2,k])){next}
      oriBetas = OriBeta[,k]
      if(any(is.na(oriBetas))){next}
      wnonzero = which(oriBetas!=0)
      numbetasvec[k] = length(wnonzero)

      if(length(wnonzero)==0){next}
      scoreinfo = cbind(refallele,oriBetas)[wnonzero,,drop=F]
      savefile2 = paste0(output_sub_folder_predi_fem, "Score.fem",k)
      write.table(scoreinfo,file=savefile2,row.names=F,col.names=F,quote=F);
      rm(scoreinfo)
      outfile_fem = paste0(output_sub_folder_predi_fem, "iter")
      command = paste0("plink",PV,"  --allow-no-sex --bfile ", paste0(output_sub_folder_predi_fem,"temp"), " --score ", savefile2,
        " --threads 1 --out ",outfile_fem, "Single.", k,".profile > 1", sep="");
      system(command);

      if(plinkver==1){
       preYplink0 = read.table(file=paste0(outfile_fem,"Single.", k,"fem.profile.profile"), as.is=T,header=T)
       file.remove(paste0(outfile_fem,"Single.", k,"fem.profile.profile"))
      }
      if(plinkver==2){
       preYplink0 = read.table(file=paste0(outfile_fem,"Single.", k,"fem.profile.sscore"), as.is=T)
       file.remove(paste0(outfile_fem,"Single.", k,"fem.profile.sscore"))
       names(preYplink0)[6] = "SCORE"
      }

      wk = which(rownames(phenoOri) %in% preYplink0[,1])
      pheno = phenoOri[wk,]
      mat = match(rownames(pheno),preYplink0[,1])
      preYplink0 = preYplink0[mat,]

      mat1 = match(testid[,1],preYplink0[,1])
      mat2 = match(valid[,1],preYplink0[,1])

      PreR2[1,k] = cor(pheno[mat1,phenoIndex],preYplink0$SCORE[mat1])^2
      PreR2[2,k] = cor(pheno[mat2,phenoIndex],preYplink0$SCORE[mat2])^2
      rm(pheno,preYplink0)
    }
    rm(OriBeta)
    gc()
    save(PreR2, numbetasvec, Finish, file=savef2)
    }
}


transETpen = function(Zscale , usedtrait, dirPlinkFormat1, sim, savefile1, low0_single, top0_single, NumIter= 1000, ld_outname, warmStart = 1, penalty , mainindex = 1, RupperVal = NULL, gwasfilename , notedSkip = 0, mafnames_vec = NULL , savefile0 = NULL, totaltuning1 = 0, totaltuning2 = 0, causalSNPs=NULL, chr , rcut = 0.6, alpha = 1,NumL= 20, subNumL = 20){

    
 usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
    usedtraitIndex = as.numeric(usedtraitsvec)
    Nvec = c(10000,10000)
   

    summaryZoutput = Get_summaryZ(gwasfilename, mainindex, Nvec, usedtraitsvec,chr)
    summaryZ = summaryZoutput[[1]]
    Allkeepsnps = summaryZoutput[[2]]
     summaryBetas = SSvec = NULL
    NamesSNP = rownames(summaryZ)
      Nq = length(Nvec)
      
    summaryBetas = matrix(0, nrow(summaryZ), Nq)
    SDvec = matrix(0, nrow(summaryZ), Nq)

    for(ii in 1:Nq){
      summaryBetas[,ii] = summaryZ[,ii]/sqrt(Nvec[ii])
      SDvec[,ii] = 1/sqrt(Nvec[ii])
    }
    rownames(summaryBetas) = rownames(summaryZ)
    SSvec = matrix(0,1,1)
    
  plinkLDList = Gen_plinkLDList(ld_outname, Allkeepsnps)
# Rename the column R to R1
names(plinkLDList)[names(plinkLDList) == "R"] <- "R1"

# Add a new column R2 that is a duplicate of R1
plinkLDList$R2 <- plinkLDList$R1


  if(is.null(top0_single)){
    print(paste0("maxZ=",maxZ))
    top0_single = maxZ
  }

    topv = top0_single
   
    Lambda_limit = c(0.02,0.9) 
   subtuning = 35
    Lenlam = 10
    taufactor = c(1/25, 1, 10)
   
 medianval = median(apply(abs(summaryBetas),1,sum),na.rm=T)
 tauvec = sort(medianval*taufactor)
    Lambda_limit = quantile(abs(summaryZ[,1]), Lambda_limit)

tuningMatrix = Tuning_setup_group_only(tauvec, subtuning, Lambda_limit, Lenlam, medianval, Nq)


   # Run the transET_PEN function with the necessary parameters
  output = transET_PEN(summaryZ= summaryZ, summaryBetas = summaryBetas, SSvec = SSvec, Nvec = Nvec, plinkLD = plinkLDList , NumIter = 500, breaking = 0, Zscale = 1,tuningMatrix = tuningMatrix, penalty="RealmixLOG",outputAll = 1, StopRule = 0 , outBefore=1)
   
   # Initialize saveoutput list and populate it with the required data
   saveoutput = list()
   saveoutput[[1]] = output$Numitervec
   saveoutput[[2]] = ncol(output$BetaMatrix)
   saveoutput[[3]] = nrow(output$BetaMatrix)
   saveoutput[[4]] = which(output$BetaMatrix != 0)
   wpos = which(output$BetaMatrix != 0)
   saveoutput[[5]] = output$BetaMatrix[wpos]
   saveoutput[[6]] = NamesSNP
   saveoutput[[7]] = output[[3]]
   saveoutput[[8]] = nrow(output$ThetaMatrix)
   saveoutput[[9]] = ncol(output$ThetaMatrix)
   saveoutput[[10]] = which(output$ThetaMatrix != 0)
   wpos1 = which(output$ThetaMatrix != 0)
   saveoutput[[11]] = output$ThetaMatrix[wpos1]


   # Save the output to the specified file
   save(saveoutput, file = savefile1)

   # Remove output and saveoutput to free up memory
   rm(output, saveoutput)
  rm(plinkLDList)
  gc()
}
 

