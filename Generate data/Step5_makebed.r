  rm(list=ls())
  iiIndex = 3
  #module load nixpkgs/16.09 plink
  #module load intel/2016.4 r/3.5.0
  popvec = c("AFR","EAS","EUR")

  dirOutputSave = "/home/yatah3/projects/def-sduchesn/yatah3/TransEthnic/SimuGenotype/"

  dirSimuDataSet = paste0(dirOutputSave,"simuSet/")

  Nvec = c("40000","40000","100000")
  Tobedfile = 1
  if(Tobedfile == 1){
    #for(iiIndex in 1:3){
    for(chr in 1:22){
      print(chr)
      dirSimuDataSetChr = paste0(dirSimuDataSet,"Chr",chr,"/")
      dirPlinkFormat = paste0(dirSimuDataSetChr,"PlinkFormat/")
      system(paste0("mkdir -p ",dirPlinkFormat))
      hapname = paste0("ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")
      pop = popvec[iiIndex]
      outputfile = paste0(dirSimuDataSetChr,pop,hapname,"allset")
      genfile = paste0(outputfile,".controls.gen")#This file typically contains the genotype data for the control individuals in a study. It represents the SNPs and their corresponding genotypes for each control individual.
      samfile = paste0(outputfile,".controls.sample")#This file contains information about the samples or individuals included in the control group. It provides details such as sample IDs, family IDs (if applicable), phenotype information, and other sample-level annotations. 


      plinkfile = paste0(dirPlinkFormat,pop,"chr",chr,"bedformat")
      if(file.exists(paste0(plinkfile,".fam"))){
        test = read.table(file=paste0(plinkfile,".fam"),as.is=T)
        if(nrow(test)==Nvec[iiIndex]){
        rm(test)
        next

        }
        rm(test)
      }
        print(paste0(pop,"chr",chr,"bedformat"))
        command3 = paste("plink2 --oxford-single-chr",chr,"--gen", genfile,"ref-unknown --sample",samfile, "--maf 0.01 --hwe 1.0e-03 --geno 0.2 --mind 0.2 --make-bed --out",plinkfile)#The Oxford single-chromosome format is in  format .fam, .bim and .bed that based on  to which chromosome the data corresponds to. 
 #      The ref-unknown option in the code ref-unknown --sample is used to indicate that reference allele information is unknown. This option is typically used when working with genotype data for which the reference alleles are not known or are not specified in the dataset.  
        system(command3)
      
     #}
  }
}
#}
