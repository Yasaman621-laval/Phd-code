  rm(list=ls())

  #module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
  #module load intel/2016.4 r/3.5.0
  popvec = c("AFR","EAS","EUR")

  dirOutputSave = "/home/yatah3/projects/def-sduchesn/yatah3/TransEthnic/SimuGenotype/"
  dirSimuDataSet = paste0(dirOutputSave,"simuSet/SimuData_2/")

  rmallele = c("AT","TA","CG","GC")

  for(chr in 1:22){
    print(chr)
    dirSimuDataSetChr = paste0(dirSimuDataSet,"Chr",chr,"/")
    dirPlinkFormat = paste0(dirSimuDataSetChr,"PlinkFormat/")
    sfile = paste0(dirPlinkFormat ,"keepSNPs_chr",chr,".txt")
    if(file.exists(sfile)){next}

    Start  = 1
    for(iiIndex in 1:length(popvec)){
      pop = popvec[iiIndex]
      plinkfile = paste0(dirPlinkFormat,pop,"chr",chr,"bedformat")
      bim = read.table(file=paste0(plinkfile,".bim"),as.is=T)
      bim$nv5 = nchar(bim$V5)#determine the length of the alleles associated with each SNP.
      bim$nv6 = nchar(bim$V6)#determine the length of the alleles associated with each SNP.
      wk = which(bim$nv5==1 & bim$nv6==1)#nv5 and nv6 are equal to 1 correspond to SNPs with single-nucleotide alleles.
      bim = bim[wk,]
      allles = apply(bim[,5:6],1,paste0,collapse="")# This line extracts the alleles of the SNPs from columns 5 and 6 of the bim data frame and concatenates them into a single string for each SNP
      wrm = which(allles %in% rmallele)
      if(length(wrm)>0){
        bim = bim[-wrm,]
      }

      if(Start==1){
        AllSNPs = bim[,2]
        Start = 0
      }else{
        AllSNPs = intersect(AllSNPs,bim[,2])
      }
      rm(bim)
      gc()
    }
    AllSNPs = unique(AllSNPs)
    print(length(AllSNPs))
    write.table(AllSNPs, file=sfile,quote = F, row.names = F,col.names = F)
    rm(AllSNPs)
 }

  totalsfile = paste0(dirSimuDataSet,"SNPs_allChrs.txt")
  if(!file.exists(totalsfile)){
    AllSNPs = c()
    for(chr in 1:22){
      dirSimuDataSetChr = paste0(dirSimuDataSet,"Chr",chr,"/")
      dirPlinkFormat = paste0(dirSimuDataSetChr,"PlinkFormat/")
      sfile = paste0(dirPlinkFormat ,"keepSNPs_chr",chr,".txt")
      AllSNPs = c(AllSNPs, read.table(file=sfile,as.is=T)[,1])
    }
    write.table(AllSNPs, file=totalsfile,quote = F, row.names = F,col.names = F)#collects all  filtered SNPs from different chromosomes and consolidates them into a single file for further analysis or reference.
    print(length(AllSNPs))
  }



  CreateMergefile = 1
  if(CreateMergefile == 1){
    for(iiIndex in 1:length(popvec)){
      pop = popvec[iiIndex]
      savefile = paste0(dirSimuDataSet,pop,"bedformatlist.txt")
      FilesNames = c()
      for(chr in 1:22){
        dirSimuDataSetChr = paste0(dirSimuDataSet,"Chr",chr,"/")
        dirPlinkFormat = paste0(dirSimuDataSetChr,"PlinkFormat/")
        outfile = paste0(dirPlinkFormat,pop,"chr",chr,"bedformat")#represents the Plink-format file for that particular population and chromosome (.bim ,.bed and .fam).
        FilesNames = c(FilesNames,outfile) 
      }
      write.table(FilesNames,file=savefile,quote=F,row.names =F, col.names = F)
      bfile = paste0(dirSimuDataSet,pop,"AllChrs_bedformat")#merged dataset that contains data from all chromosomes for a specific population.
      command4 = paste("plink --allow-no-sex --merge-list",savefile, "--extract", totalsfile,"--make-bed --out",bfile)# merged plink format dataset that contains data from all the chromosomes for each population based on only filtered SNPs (totalsfile).The "--extract" option in Plink is used to include specific SNPs in the analysis, not exclude them.
      system(command4)
      test = read.table(file=paste0(bfile,".bim"),as.is=T)
      print(dim(test))

      fam = read.table(file=paste0(bfile,".fam"),as.is=T)
      print(dim(test))

    }
  }






