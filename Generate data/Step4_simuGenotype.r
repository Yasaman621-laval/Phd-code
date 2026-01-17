
iiIndex = 3


chrS = 22
chrE = 1
if(iiIndex==3){

#chrsvec = c(2,22,21)
#chrsvec = c(3,20,19)
chrsvec = c(1:22)
}else{
chrsvec = c(chrS:chrE)
}

chrsvec

#---subset haplotype files
#module load nixpkgs/16.09 intel/2016.4 r/3.5.0
# module load shapeit/2.r904
# module load hapgen2


#--using IMPUTE haps

  dirOutputSave = "/home/yatah3/projects/def-sduchesn/yatah3/TransEthnic/SimuGenotype/"


  dir0 = "/home/yatah3/projects/def-thchlava/2023_0227/"

  dirSimuData = paste0(dir0,"Data1000Genomes/SimuGenotype/")


  popvec = c("AFR","EAS","EUR")
  Nvec = c("40000","40000","100000")
  NvecN = c(40000,40000,100000)

#for(iiIndex in 1:3){
  paste0(iiIndex)
  pop = popvec[iiIndex]
  Ninput = Nvec[iiIndex]

for(chr in chrsvec){


  dirSimuDataSet = paste0(dirOutputSave,"simuSet/")
  system(paste("mkdir -p",dirSimuDataSet))


  dirTag = paste0(dirSimuDataSet,"Tag/")
  system(paste("mkdir -p",dirTag))

  dirSimuDataSetChr = paste0(dirSimuDataSet,"Chr",chr,"/")
  system(paste("mkdir -p",dirSimuDataSetChr))

  hapname = paste0("ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")

  legendfile = paste0(dirOutputSave,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypesFirst6cols.txt")
  legend = read.table(file=legendfile,as.is=T,header=F)
  legend = legend[order(legend[,3]),]
  dim(legend)
  legend = legend[,2:5]
  colnames(legend) = c("rs", "position", "X0", "X1")

  sublegendfile = paste0(dirTag,"ALL.chr",chr,".leg")
  if(!file.exists(sublegendfile)){
    write.table(legend,file=sublegendfile,quote=F,row.names =F, col.names = T)
  }


  tagsfile = paste0(dirTag,"Tag_chr",chr,".tags")
  if(!file.exists(tagsfile)){
    tages = legend[,2,drop=F]
    write.table(tages,file=tagsfile,quote=T,row.names =F, col.names = F)
  }

  tag = read.table(file=tagsfile,as.is=T)[,1]
  outputfile = paste0(dirSimuDataSetChr,pop,hapname,"allset")

  samplefile = paste0(outputfile,".controls.sample")
  if(file.exists(samplefile)){

  sampleID = read.table(file=samplefile,as.is=T,header=T)


  if(nrow(sampleID) >= NvecN[iiIndex]){next}

  print(samplefile)
  }
  #stop("")

  mapfile = paste0(dir0,"Data1000Genomes/map/genetic_map_chr",chr,"_combined_b37.20140701.txt")
  hapfile2 = paste0(dirOutputSave,"notitle",pop,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.haps")

  alld = paste(tag[1:10], "1 1 1")
  alld = paste0(alld,collapse=" ")

  command2 = paste("hapgen2 -m",mapfile, "-l", sublegendfile, "-h", hapfile2, "-o",outputfile,"-dl ",alld," -n ",Ninput," 0")
  system(command2)
}

