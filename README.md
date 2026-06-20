# PhD Code: Cross-Population Genetic Architecture and Polygenic Prediction

This repository contains the code used for simulation studies and real-data analyses in a PhD project on **cross-population genetic architecture estimation** and **summary-statistics-based polygenic prediction**. The code is built around GWAS summary statistics, linkage disequilibrium (LD) information, sparse effect-size estimation, posterior sharing estimation, and comparison with external PRS methods.

The repository supports two connected goals:

1. **Genetic architecture estimation**: estimating how SNP effects are distributed across population-specific, pairwise shared, and fully shared components.
2. **Polygenic prediction evaluation**: estimating SNP effects from summary statistics and evaluating prediction performance, including comparison with PRS-CSx.

The main analyses are implemented in `R`, with supporting shell scripts, PLINK commands, Python utilities, and a local `SummaryLasso` R package.

---

## Repository Status

This is a research-code repository. The scripts reproduce the thesis analyses when the required input data, reference panels, and directory paths are available. Most large data files are not included because they are external GWAS or genotype reference files.

Several scripts were written for a high-performance computing environment and contain absolute paths to Compute Canada / cluster directories. Before running the code, update all path variables to match your local or cluster environment.

---

## Main Workflow

The full analysis pipeline has the following structure:

1. Prepare simulated genotype data and GWAS-like summary statistics.
2. Harmonize and quality-control real GWAS summary statistics.
3. Match GWAS variants to ancestry-specific reference panels.
4. Compute LD and create clumped SNP sets.
5. Estimate sparse SNP effects using summary-statistics-based penalized regression.
6. Estimate cross-population sharing proportions using ES-alpha, posterior soft counts, and calibrated posterior methods.
7. Evaluate prediction performance using AUC or related metrics.
8. Compare the proposed method with PRS-CSx.

---

## Repository Structure

```text
Phd-code/
├── Generate data/              # Simulation genotype, phenotype, LD, and GWAS preparation
├── Prepration real_data/        # Real-data GWAS harmonization, QC, clumping, and LD preparation
├── Project1/                   # Genetic architecture estimation
│   ├── Simulation/              # Simulation analyses for sharing-proportion estimation
│   └── Real data/               # Real-data architecture estimation for T2D
├── Project2/                   # Polygenic prediction and performance evaluation
│   ├── Simulation/              # Simulation-based prediction evaluation
│   └── Real data/               # Real-data prediction and AUC evaluation
├── Rfunction/                  # Shared helper functions used across analysis scripts
├── SummaryLasso/               # Local R package for summary-statistics-based penalized regression
├── SummaryAUC/                 # AUC and prediction-performance utilities
├── PRS-CSX/                    # PRS-CSx baseline comparison workflow
├── run_mixer/                  # MiXeR-style LD and summary-statistic utilities
└── README.md
```

---

## Directory Descriptions

### `Generate data/`

This folder contains scripts for constructing simulation datasets from reference genotype data and generating GWAS-like summary statistics.

Main tasks include:

- Selecting SNPs from reference panels.
- Simulating genotypes and phenotypes.
- Generating SNP effects under different heritability and correlation settings.
- Preparing LD information.
- Combining and standardizing GWAS summary statistics.

Important scripts include:

```text
Step0_getPopulationID.r
Step1_selectSNPs.r
Step2_subsetHaps.r
Step4_simuGenotype.r
Step5_makebed.r
Step6_clean_bed .r
step6A_simu_y.r
step6B_simu_y.r
step6A_simu_y_3traits.r
step6B_simu_y_3traits.r
step7B_generate_Beta.r
step7B_generate_Beta_3traits.r
step7_ld.r
step7_ld_3traits.r
step8B_combine_GWASstandard.r
step8B_combine_GWASstandard_3traits.r
```

These scripts are mainly used before the project-level simulation analyses.

---

### `Prepration real_data/`

This folder contains the real-data preparation pipeline for type 2 diabetes GWAS summary statistics.

Main tasks include:

- Harmonizing AFR, EAS, and EUR GWAS summary statistics.
- Converting odds ratios to log odds ratios where needed.
- Reconstructing Z-scores and standard errors.
- Merging SNPs across populations.
- Removing ambiguous and low-quality SNPs.
- Matching GWAS variants to ancestry-specific 1000 Genomes reference panels.
- Performing LD clumping and preparing LD files.

Important scripts include:

```text
Step1_Harmonized.r
Step2_(step0_to_step8_ld).r
```

Expected real-data inputs include ancestry-specific GWAS summary statistics and PLINK-format 1000 Genomes reference panels for AFR, EAS, and EUR.

---

### `Project1/`

`Project1` focuses on **genetic architecture estimation**. The goal is to estimate the vector of sharing proportions across three dimensions, such as:

```text
pi0, pi1, pi2, pi3, pi12, pi13, pi23, pi123
```

For the cross-population T2D application, the population order is:

```text
AFR, EAS, EUR
```

Therefore, the sharing components correspond to null, population-specific, pairwise shared, and fully shared SNP-effect patterns.

#### `Project1/Simulation/`

This folder contains simulation scripts for estimating and evaluating genetic sharing proportions under known truth.

Important scripts include:

```text
Estimate_pi_3pop
Estimate_pi_3traits.r
optim_3pop.r
Optim_3traits.r
es-alphamatrix_sim_3pops.r
es-alphamatrix_sim_3traits.r
3pop_opt_real_data.r
```

Main outputs include:

- ES-alpha estimates.
- Posterior soft-count estimates.
- Optimized or calibrated estimates of the sharing vector.
- RMSE-based comparison against the true simulation values.

#### `Project1/Real data/`

This folder applies the architecture-estimation pipeline to real T2D GWAS summary statistics.

Important scripts include:

```text
run_step1_Trans_mixLog_3pops_project1
run_step2_densityU_3pops_project1.r
run_step3_Estimate_pi_3pops_project1.r
run_step3_OPTimized_Estimate_pi_3pops_project1.r
```

Main tasks include:

- Fitting sparse effect-size models.
- Computing posterior component support.
- Estimating sharing proportions from posterior soft counts.
- Computing ES-alpha as a comparison.
- Applying simulation-calibrated correction parameters to obtain the final real-data architecture estimate.

---

### `Project2/`

`Project2` focuses on **polygenic prediction and performance evaluation**. It uses fitted SNP effects to evaluate prediction performance across population settings and compares the proposed method with PRS-CSx.

#### `Project2/Simulation/`

This folder contains simulation-based prediction analyses.

Important scripts include:

```text
run_step1_Trans_mixLog_3pops_project1.r
run_step2_densityU.r
step3_preR2.r
step4_results_R2_TH.r
mapspen_PRS_CSX_Plot.r
```

Main tasks include:

- Fitting the proposed summary-statistics model.
- Generating prediction weights.
- Computing prediction-performance metrics.
- Comparing MAPSPEN-style estimates with PRS-CSx outputs.

#### `Project2/Real data/`

This folder contains the real-data prediction and performance evaluation workflow.

Important scripts include:

```text
run_step1_Trans_mixLog.r
run_step2_densityU.r
step3_AUC_computation-2.r
step4.r
step5.r
step5_MAPSpen_vs_pRS_CSX.r
```

Main outputs include:

- AUC summaries.
- Best-performing tuning-parameter results.
- Comparison plots for the proposed method and PRS-CSx.

---

### `Rfunction/`

This folder contains shared helper functions used by the project scripts.

Important files include:

```text
AllPRS_Rfunctions.r
Internal_Rfunctions.r
Iterative_Rfunctions.r
PRS_utility.r
PlinkLD_transform.R
```

These functions support:

- Summary-statistics penalized regression.
- Iterative optimization.
- LD transformation.
- PRS construction.
- Common utility operations.

---

### `SummaryLasso/`

This folder contains a local R package for summary-statistics-based penalized regression.

The package includes:

```text
R/        # R functions
src/      # Compiled C source code
man/      # Documentation
```

Install locally before running scripts that load `SummaryLasso`:

```bash
R CMD INSTALL SummaryLasso
```

---

### `SummaryAUC/`

This folder contains utilities for computing AUC and related prediction-performance summaries from PRS outputs.

It includes:

- `auc.R`
- C source and compiled objects for adjusted correlation calculations.
- Example model files and output files.

If recompilation is required, rebuild the shared object file for your own system rather than relying on the included compiled file.

---

### `PRS-CSX/`

This folder contains scripts used to run PRS-CSx as a baseline comparison method.

Main tasks include:

- Creating PRS-CSx-compatible SNP information files.
- Creating LD blocks.
- Running PRS-CSx with and without meta-analysis.
- Summarizing PRS-CSx prediction results.

The folder includes the PRS-CSx Python implementation and related R wrapper scripts.

---

### `run_mixer/`

This folder contains MiXeR-style utilities for LD and summary-statistic preparation.

Main tasks include:

- Creating `.snps` files.
- Creating LD files by chromosome and ancestry.
- Running summary-statistic processing steps for cross-population and multi-trait analyses.

Several scripts are designed for SLURM array jobs.

---

## Software Requirements

The main analysis environment used:

```text
R version 4.3.1
PLINK version 1.9b_6.21-x86_64
```

Some legacy comments in the scripts refer to older module versions. The thesis analyses were run using the versions above unless otherwise specified.

Required R packages include:

```r
install.packages(c(
  "data.table",
  "dplyr",
  "stringr",
  "readr",
  "ggplot2",
  "tidyr",
  "tibble",
  "purrr",
  "MASS",
  "mvtnorm",
  "gtools",
  "foreach",
  "iterators",
  "doParallel",
  "jsonlite",
  "parallel",
  "openxlsx"
))
```

The local `SummaryLasso` package should also be installed:

```bash
R CMD INSTALL SummaryLasso
```

Additional external tools used in parts of the pipeline include:

```text
PLINK / PLINK2
HAPGEN2
SHAPEIT
Python 3
PRS-CSx
MiXeR
```

---

## Data Requirements

Large input data files are not included in this repository. The user must provide or generate:

- GWAS summary statistics for each population or trait.
- 1000 Genomes or other ancestry-matched reference panels in PLINK format.
- Simulated genotype files, if running the full simulation pipeline.
- Train/test summary-statistic files for prediction analyses.
- Any external PRS-CSx or MiXeR reference data required by those tools.

For real-data T2D analyses, the expected populations are:

```text
AFR = African ancestry
EAS = East Asian ancestry
EUR = European ancestry
```

The real-data preparation scripts assume harmonized summary statistics with fields such as:

```text
rsid
CHR
POS
Effect_Allele
Other_Allele
beta_POP
se_POP
pval_POP
Z_POP
N_POP
EAF_POP
```

where `POP` is one of `AFR`, `EAS`, or `EUR`.

---

## Running the Code

Because the scripts were designed for an HPC environment, there is no single command that runs the entire project from start to finish. The scripts should be run step by step after updating input and output paths.

### 1. Install the local package

```bash
R CMD INSTALL SummaryLasso
```

### 2. Update paths

Before running any script, edit path variables such as:

```r
dirBase
out_dir
plink_dir
bfile_EUR
bfile_AFR
bfile_EAS
functionsfolder
output_sub_folder
```

Several scripts contain absolute paths such as `/lustre...` or `/home/yatah3/...`; these must be changed for a new environment.

### 3. Real-data preparation example

```bash
Rscript "Prepration real_data/Step1_Harmonized.r"
Rscript "Prepration real_data/Step2_(step0_to_step8_ld).r"
```

### 4. Architecture estimation example

```bash
Rscript "Project1/Real data/run_step1_Trans_mixLog_3pops_project1"
Rscript "Project1/Real data/run_step2_densityU_3pops_project1.r"
Rscript "Project1/Real data/run_step3_Estimate_pi_3pops_project1.r"
Rscript "Project1/Real data/run_step3_OPTimized_Estimate_pi_3pops_project1.r"
```

### 5. Prediction evaluation example

```bash
Rscript "Project2/Real data/run_step1_Trans_mixLog.r"
Rscript "Project2/Real data/run_step2_densityU.r"
Rscript "Project2/Real data/step3_AUC_computation-2.r" <job_index>
Rscript "Project2/Real data/step5_MAPSpen_vs_pRS_CSX.r"
```

Scripts that require `<job_index>` are intended for SLURM array execution and expect an `input.txt` file defining scenario or parameter combinations.

---

## Main Outputs

Depending on the workflow, the scripts generate:

- Harmonized GWAS summary statistics.
- Merged cross-population SNP files.
- LD-clumped SNP sets.
- Population-specific LD files.
- Sparse fitted SNP-effect estimates.
- Posterior component-support files.
- ES-alpha sharing estimates.
- Posterior soft-count sharing estimates.
- Final calibrated sharing estimates.
- Prediction-performance summaries such as AUC.
- Comparison plots for MAPSPEN-style results and PRS-CSx.

---

## Important Notes

- Raw GWAS and genotype data are not included.
- Some folders contain third-party tools or comparison methods, such as PRS-CSx and SummaryAUC. These retain their original licenses and should be cited separately where appropriate.
- Compiled files, such as `.so` files, may need to be rebuilt on a different system.
- Directory names containing spaces, such as `Generate data/` and `Prepration real_data/`, should be quoted when used in shell commands.
- The repository is designed for reproducibility of thesis analyses after paths and input data are configured.

---

## Citation

If you use this code, please cite the related PhD thesis and the relevant external methods used in the workflow, including PRS-CSx, PLINK, MiXeR, and SummaryLasso where applicable.

---

## Contact

For questions about this repository, please contact:

```text
Yasaman Tahernezhad
GitHub: Yasaman621-laval
```
