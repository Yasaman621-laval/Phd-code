
#---subset haplotype files
# module load nixpkgs/16.09 intel/2016.4 r/3.5.0

# module load shapeit/2.r904


#--using IMPUTE haps

dirOutputSave = "/home/yatah3/projects/def-sduchesn/yatah3/TransEthnic/"

dir0 = "/home/yatah3/projects/def-thchlava/2023_0227/"
dirIMPUTE = paste0(dir0,"Data1000Genomes/IMPUTE/1000GP_Phase3/")
dirhapsample = paste0(dir0,"Data1000Genomes/hapsample/")
popvec = c("AFR","EAS","EUR")


Chrvec = c(1:22)
Chrvec
  chr = Chrvec
for(iiIndex in c(1:3)){
  for(chr in 1:22){
  legendfile = paste0(dirOutputSave,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypesFirst6cols.txt")

  hapname = paste0("ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")
  hapfile = paste0(dirOutputSave,popvec[iiIndex],hapname,".haps")

  if(!file.exists(legendfile)){
    command = paste0("awk \'{ print $1,$2,$3,$4,$5,$6   }\' ",hapfile," > ",legendfile)
    system(command)
  }
  }
}

for(ii in 1:length(popvec)){
  for(chr in Chrvec){
    print(paste0("chr",chr))
    hapname = paste0("ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")
    hapfile = paste0(dirOutputSave,popvec[ii],hapname,".haps")
    haps= read.table(file=hapfile,header=F,as.is=T,nrow=6)
    colprint = paste0(c(6:ncol(haps)),collapse=",$")

    hapfile2 = paste0(dirOutputSave ,"notitle",popvec[ii],"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.haps")

    command = paste0("awk \'{ print $",colprint,"   }\' ",hapfile," > ",hapfile2)
    system(command)
  }
}
