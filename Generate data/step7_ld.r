    rm(list=ls())
    library(MASS)
    library(mvtnorm)
    library(gtools)
    testset = 1
    #module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
    #module load intel/2016.4 r/3.5.0

    popvec = c("AFR","EAS","EUR")

    dirOutputSave = "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/"
    dirSimuData = paste0(dirOutputSave,"sim_hsq0.3_rho0.2_train5000-5000-50000/")


    totalsfile = usedsnpfile = paste0(dirSimuData,"SNPs_allChrs.txt")

    refAllelefile = paste0(dirSimuData,"refereceAlleleSNPs.txt")
    if(!file.exists(refAllelefile)){
      pop = popvec[3]
      plinkfile = paste0(dirSimuData,pop,"AllChrs_bedformat.bim")
      bim = read.table(file=plinkfile,as.is=T)
      refSNP = bim[,c(2,5)]#the SNP ID is in the second column and the reference allele is in the fifth column. The resulting reference SNP information is stored in the refSNP variable.
      rm(bim)
      write.table(refSNP,file=refAllelefile,quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
    }


    LD_Index = 1
    numsnp = 400
    if(LD_Index==1){
      for(iiIndex in 3){
        pop = popvec[iiIndex]
        for(chr in 1:22){
        print(chr)
        dirSimuDataSetChr = paste0(dirSimuData,"Chr",chr,"/")
        dirPlinkFormat = paste0(dirSimuDataSetChr,"PlinkFormat/")
        plinkfile = paste0(dirPlinkFormat,pop,"chr",chr,"bedformat")
        outname = paste0(dirPlinkFormat,pop,"chr",chr,"bedformat_numsnp",numsnp)
        
        lineR = paste0("plink --noweb --allow-no-sex --bfile ",plinkfile ," --reference-allele ",refAllelefile ," --extract ",usedsnpfile," --r --ld-window ",numsnp," --out ", outname)#--noweb : It ensures that Plink operates solely with the data and files available locally on your system not online resources.
        system(lineR)

      }
    }
  }
