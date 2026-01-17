rm(list=ls())
iiIndex = 2
chrS = 4
chrE = 22

#chrS = 10
#chrE = 17
#iiIndex = 3
#chrS = 18
#---subset haplotype files
#module load nixpkgs/16.09 intel/2016.4 r/3.5.0
# module load shapeit/2.r904



#--using IMPUTE haps

dirOutput = "/home/yatah3/projects/def-sduchesn/yatah3/randomsam/TransEthnic/"
dirOutputSave = "/home/yatah3/projects/def-thchlava/TransEthnic/SimuGenotype/"

dir0 = "/home/yatah3/projects/def-sduchesn/yatah3/randomsam/Research/"
dirGenome = paste0(dir0,"Data1000Genomes/")
dirIMPUTE = paste0(dir0,"Data1000Genomes/IMPUTE/1000GP_Phase3/")
dirGunZip = paste0(dir0,"Data1000Genomes/IMPUTE/1000GP_Phase3/gunzip/")
#dirhapsample = paste0(dir0,"Data1000Genomes/hapsample/")

dirhapmap = "/home/yatah3/projects/def-sduchesn/yatah3/randomsam/PRSdata/phase_3/"

  popvec = c("AFR","EAS","EUR")

  for(ii in iiIndex){
    for(chr in c(chrS:chrE)){
    keepposfile = paste0(dirhapmap,"filtered_hapmap3_chr",chr,"pos.txt")

    print(paste0("chr=",chr))
    filterhap3file = paste0(dirOutputSave,"filtered_hapmap3_chr",chr,".txt")
    hapname = paste0("ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")
    hapfile = paste0(filterhap3file,hapname)
    outputfile = paste0(dirOutputSave,popvec[ii],hapname)
    indsfile = paste0(dirOutput,popvec[ii],"1000GP_Phase3.inds")
    command1 = paste("shapeit -convert --input-haps", hapfile, "--output-haps",outputfile,"--include-ind",indsfile,"--include-snp",keepposfile, "--thread 12")
    system(command1)
 }
}
