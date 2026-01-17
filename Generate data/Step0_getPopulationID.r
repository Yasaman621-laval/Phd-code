
  dirOutput = "/home/yatah3/projects/def-sduchesn/yatah3/TransEthnic/"

#---subset haplotype files
# module load shapeit/2.r904
#--using IMPUTE haps


  dir0 = "/home/yatah3/projects/def-thchlava/2023_0227/"
  dirIMPUTE = paste0(dir0,"Data1000Genomes/IMPUTE/1000GP_Phase3/")
  dirhapsample = paste0(dir0,"Data1000Genomes/hapsample/")
  dirSave = paste0(dirhapsample,"SubsetHaps/")
  dirSubsetLegendTags = paste0(dir0,"Data1000Genomes/IMPUTE/SubsetLegendTags/")

  dirOutputSave = "/home/yatah3/projects/def-thchlava/TransEthnic/SimuGenotype/"
  system(paste0("mkdir -p ",dirOutputSave))


  samsIMPUTEfile = paste0(dirIMPUTE,"1000GP_Phase3.sample")
  samsIMPUTE = read.table(file=samsIMPUTEfile,as.is=T,header=T)


  popvec = c("AFR","AMR","EAS","EUR","SAS")

  for(ii in 1:length(popvec)){
    subsams = samsIMPUTE$ID[which(samsIMPUTE$GROUP==popvec[ii])]
    indsfile = paste0(dirOutputSave,popvec[ii],"1000GP_Phase3.inds")
    ids = cbind(0,subsams)
    write.table(ids,file=indsfile,quote=F,row.names =F, col.names = F)
  }

AFR AMR EAS EUR SAS
661 347 504 503 489
