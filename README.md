# Phd-code  
**Cross-population & multi-trait Polygenic Risk Score pipelines**  

---

## Overview

This repository contains all code used in my PhD thesis for developing and evaluating **penalized polygenic risk score (PRS) methods using GWAS summary statistics**, with a focus on:

- Cross-population PRS transferability  
- Multi-trait information sharing  
- Robust modeling of heterogeneous genetic architectures  

The code directly implements the methodological workflow described in:

- **Chapter 1**: Simulation design, GWAS construction, and data processing  
- **Chapter 2**: MAPSPEN methodology, nuisance parameter estimation, and empirical evaluation  

Both **simulation studies** and **real-data analyses** are included, with **PRS-CSx** used as a baseline comparator throughout.

---

## High-level workflow (thesis → code)

The thesis workflow is implemented in the following order:

1. GWAS input preparation and QC  
2. Estimation of genetic architecture / nuisance parameters  
3. Penalized PRS modeling (MAPSPEN / SSPEN)  
4. Baseline PRS-CSx construction  
5. Predictive performance evaluation (AUC / R²)  
6. Method comparison under different sharing regimes  

This structure is mirrored exactly in the repository.

---

## Repository structure
```
Phd-code/
│
├── Generate data/ # Simulation data generation & GWAS standardization
├── Project1/
│ ├── Simulation/ # Simulation experiments (Project 1)
│ └── Real data/ # Real-data pipeline (Project 1)
├── Project2/
│ ├── Simulation/ # Simulation experiments (Project 2)
│ └── Real data/ # Real-data pipeline (Project 2)
├── PRS-CSX/ # PRS-CSx baseline implementation
├── run_mixer/ # MiXeR-style nuisance parameter estimation
├── Rfunction/ # Shared R helper functions
├── SummaryAUC/ # Predictive performance evaluation
├── SummaryLasso/ # Penalized regression R package
└── README.md
```

---

## Folder descriptions
## 1. `Generate data/`  
**(Chapter 1: Simulation setup & GWAS construction)**

This folder contains scripts used to construct **simulation inputs** and **standardized GWAS-like summary statistics**.

Main tasks:
- Phenotype simulation under predefined heritability and correlation structures  
- Multi-population and multi-trait GWAS generation  
- GWAS harmonization and standardization  
- LD preparation  

Typical scripts:
- `step6A_sim*_y*.r` – phenotype simulation  
- `step7_ld*.r` – LD block preparation  
- `step8B_combine_GWASstandard*.r` – GWAS harmonization  

These scripts implement the simulation framework described in **Chapter 1**, including sparse causal SNPs and mixture-based effect sharing.

---

## 2. `Project1/` and `Project2/`

Both projects follow the same internal logic but represent **different experimental settings**.

---

### 2.1 `Simulation/`  
**(Chapter 1 & Chapter 2: controlled experiments)**

These folders contain step-wise simulation pipelines.  
File names indicate execution order.

Typical workflow:
run_step1_.r # Effect transformation / mixture-logit modeling
run_step2_.r # Effect-size density estimation
step3_.r # Parameter estimation
step4_.r # PRS construction
step5_*.r # Method comparison and visualization


Examples:
- `run_step1_Trans_mixLog_3pop.r`
- `densityU_3traits.r`
- `boxplot_pi_3pops.r`

These scripts are used to:
- Estimate mixture proportions (π)  
- Recover true sharing structures  
- Compare MAPSPEN, SSPEN, and PRS-CSx under different architectures  

---

### 2.2 `Real data/`  
**(Chapter 1 & Chapter 2: real-data analysis)**

These folders implement the **real-data PRS pipeline**, using GWAS summary statistics (e.g., T2D).

Example scripts:
- `run_step1_Trans_mixLog.r`
- `run_step2_densityU.r`
- `step3_AUC_computation-2.r`
- `step5_MAPSPen_vs_pRS_CSX.r`

Pipeline logic:
1. Harmonize and QC GWAS summary statistics  
2. Estimate effect-size distributions  
3. Fit MAPSPEN / SSPEN  
4. Construct PRS  
5. Compare predictive performance against PRS-CSx  

This corresponds directly to the real-data analyses in **Chapter 2**.

---

## 3. `run_mixer/`  
**(Chapter 2: nuisance parameter estimation)**

This folder contains **MiXeR-style estimation code** used to infer:

- Number of causal variants  
- Cross-population / cross-trait overlap  
- Genetic architecture parameters  

Contents include:
- Shell scripts for SNP and LD preparation  
  - `mixer_ld_*.sh`
  - `mixer_snps_*.sh`
- R wrappers
  - `run_step1_runMixer_crossPop.r`
  - `run_step2_runMixer_summarystatistics*.r`

These outputs are used as **inputs to MAPSPEN**, not as final PRS models.

---

## 4. `PRS-CSX/`  
**(Baseline comparison method)**

This folder contains:
- Official PRS-CSx Python implementation (`PRScsx/`)
- Helper scripts for:
  - SNP matching  
  - LD block creation  
  - Running PRS-CSx with and without meta-analysis  

PRS-CSx results are used exclusively for benchmarking in **Chapter 2**.

---

## 5. `Rfunction/`  
**(Shared utilities)**

Reusable R functions shared across all pipelines:

- `AllPRS_Rfunctions.r` – PRS construction helpers  
- `Internal_Rfunctions.r` – internal utilities  
- `Iterative_Rfunctions.r` – iterative optimization  
- `PRS_utility.r` – common PRS calculations  
- `PlinkLD_transform.R` – LD processing  

Most scripts in `Project1/` and `Project2/` source functions from this folder.

---

## 6. `SummaryAUC/`  
**(Predictive performance evaluation)**

This folder computes:
- AUC (and standard error) across p-value thresholds  
- Performance plots and summary tables  

Used for real-data experiments and corresponds to evaluation figures in the thesis.

---

## 7. `SummaryLasso/`  
**(Penalized regression implementation)**

This folder is an **R package** implementing summary-statistics-based penalized regression:

- `R/` – R interface  
- `src/` – C++ backend  
- `man/` – documentation  

MAPSPEN / SSPEN rely on this package internally.

---

## Recommended execution order

### Simulation experiments
1. `Generate data/`
2. `run_mixer/` (optional, for architecture estimation)
3. `Project*/Simulation/` (step1 → step5)
4. `SummaryAUC/`

### Real-data analysis
1. GWAS QC & harmonization (`Project*/Real data/`)
2. `run_mixer/`
3. MAPSPEN / SSPEN fitting
4. `PRS-CSX/` baseline
5. `SummaryAUC/`
6. `step5_MAPSPen_vs_pRS_CSX.r`

---

## Notes

- Script names encode execution order (`step1`, `step2`, …).
- Simulation and real-data pipelines share logic but use different inputs.
- All figures and tables in **Chapters 1–2** of the thesis can be reproduced using this repository.

---
