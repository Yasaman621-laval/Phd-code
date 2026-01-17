  rm(list=ls())
  #module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
  #module load intel/2016.4 r/3.5.0

  dirOutput ="/home/yatah3/projects/def-sduchesn/yatah3/TransEthnic/"
  dirBed = "/home/yatah3/projects/def-thchlava/2023_0227/Data1000Genomes/Bed1000Genome/"

  dirOutputSave = "/home/yatah3/projects/def-thchlava/TransEthnic/SimuGenotype/"
  system(paste0("mkdir -p ",dirOutputSave))


  chr=22
  mafcut = 0.01
  missingcut = 0.02
  hwecut = 1e-6

  dirhapmap = "/home/yatah3/projects/def-sduchesn/yatah3/PRSdata/phase_3/"# ".map" extension is often used for files that contain information about genetic markers, such as SNPs (Single Nucleotide Polymorphisms), and their physical positions on a reference genome. These files typically include columns specifying the chromosome number, SNP identifier, genetic position, and physical position
  system(paste0("mkdir -p ",dirhapmap ))

  sfiles = list.files(dirhapmap,pattern="qc.poly.recode.map")

  hap3file = paste0(dirhapmap,"hapmap3.txt")
  if(!file.exists(hap3file)){
    for(s in 1:length(sfiles)){
      map = read.table(file=paste0(dirhapmap,sfiles[s]),as.is=T)
      if(s==1){
        AllSNPs = map[,2]
      }else{
        AllSNPs = intersect(AllSNPs,map[,2])#find the common elements between AllSNPs and map[,2]
      }
      rm(map)
    }
    AllSNPs = unique(AllSNPs)
    write.table(AllSNPs, file=hap3file,quote = F, row.names = F,col.names = F)
  }


  popvec = c("AFR","EAS","EUR")

  for(chr in 1:22){
    for(ii in 1:length(popvec)){

      indsfile = paste0(dirOutputSave,popvec[ii],"1000GP_Phase3.inds")
      bedfile = paste0(dirBed,"ALL.chr",chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")#related to genotypes derived from the "phase 3 based  on haplotype phasing using the Shapeit2 software
      outname = paste0(dirOutputSave,"ALL.chr",chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes",popvec[ii])#contains the original haplotype data. This file typically contains information about the filtered SNPs based on  missingcut, mafcut, and hwecut and the corresponding haplotypes for a set of individuals.
      if(file.exists(paste0(outname,".fam"))){next}

      lineR = paste0("plink --allow-no-sex --bfile ",bedfile, " --keep ", indsfile," --geno ",missingcut," --maf ",mafcut," --hwe ",hwecut," --extract ",hap3file," --make-bed --out ", outname)#use PLINK software to performs various operations such as data filtering, quality control, and extraction of specific SNPs for each population in .bim format,
      system(lineR)
    }
  }

    for(chr in 1:22){
       AllSNPs = c()
       filterhap3file = paste0(dirOutputSave,"filtered_hapmap3_chr",chr,".txt")

       if(file.exists(filterhap3file)){next}
       Start  = 1
       for(ii in 1:length(popvec)){
         outname = paste0(dirOutputSave,"ALL.chr",chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes",popvec[ii])
         bim = read.table(file=paste0(outname,".bim"),as.is=T)
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
      write.table(AllSNPs, file=filterhap3file,quote = F, row.names = F,col.names = F)# extract filtered SNPS from  outname file which is in .bim format and is related to genotypes derived from the "phase 3 based  on haplotype phasing for each population and chr. These SNPs are stored in filterhap3file
    }

    allpose = c()
    for(chr in 1:22){
       posfile = paste0(dirOutputSave,"filtered_hapmap3_chr",chr,"pos.txt")
       if(file.exists(posfile)){next}
       filterhap3file = paste0(dirhapmap,"filtered_hapmap3_chr",chr,".txt")
       AllSNPs = read.table(file=filterhap3file,as.is=T)

       bedfile = paste0(dirBed,"ALL.chr",chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes")

       bim = read.table(file=paste0(bedfile,".bim"),as.is=T)

       pos = bim[which(bim[,2] %in% unlist(AllSNPs)),4]#This line extracts a subset of positions from the .bim file. It uses the which() function to identify the rows in the .bim file where the second column (bim[,2]) matches any of the SNPs listed in the AllSNPs vector. The unlist() function is used to convert the AllSNPs vector into a regular vector before performing the matching. Finally, the resulting subset of positions is assigned to the pos variable, which likely represents genomic positions corresponding to the selected SNPs.


       write.table(pos,file=posfile,quote = F, row.names = F,col.names = F)

       print(length(pos))
       allpose = c(allpose,pos)

       rm(bim, AllSNPs,pos)
    }


    length(allpose)




