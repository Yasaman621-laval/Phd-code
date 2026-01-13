SC_DIR="/lustre06/project/6005709/yatah3/simulation/SimuGenotype/sim_hsq0.5_rho0.4_train10000-10000-80000"
LD_DIR="${SC_DIR}/ld_mixer"

# AFR
cat ${LD_DIR}/AFR_chr*.prune_maf0p05_rand2M_r2p8.snps \
  > ${LD_DIR}/AFR_allchr.prune_maf0p05_rand2M_r2p8.snps

# EAS
cat ${LD_DIR}/EAS_chr*.prune_maf0p05_rand2M_r2p8.snps \
  > ${LD_DIR}/EAS_allchr.prune_maf0p05_rand2M_r2p8.snps

# EUR
cat ${LD_DIR}/EUR_chr*.prune_maf0p05_rand2M_r2p8.snps \
  > ${LD_DIR}/EUR_allchr.prune_maf0p05_rand2M_r2p8.snps

# quick sanity check
wc -l ${LD_DIR}/*_allchr.prune_maf0p05_rand2M_r2p8.snps
head -n 3 ${LD_DIR}/AFR_allchr.prune_maf0p05_rand2M_r2p8.snps
