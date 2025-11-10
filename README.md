# Genomic Prediction in Quinoa: Comparing GBLUP and Machine Learning Approaches

This repository contains code and data for genomic prediction analysis in quinoa (*Chenopodium quinoa*), comparing traditional genomic best linear unbiased prediction (GBLUP) with a machine learning method (LightGBM).

## Overview

This project implements two complementary approaches for genomic prediction:

1. **GBLUP** - Traditional quantitative genetics approach using genomic relationship matrices
2. **LightGBM** - Machine learning approach using gradient boosting with principal components

Both methods were evaluated using group k-fold cross-validation for new accessions across multiple traits and environments in Australian and Pakistani quinoa trials.

## Repository Structure

```
.
├── README.md
├── data/
│   ├── AUS_phenotypes_raw.csv                    # Raw phenotypes (Australia)
│   ├── PAK_phenotypes_raw.csv                    # Raw phenotypes (Pakistan)
│   ├── AUSPAK_phenotypes_means_BLUEs.csv         # Processed phenotypes (BLUEs)
│   ├── quinoa_551_AUSPAK_PCs.csv                 # Principal components
│   ├── kinship_matrix_VanRaden_auspak_maxmissing20.RData  # Genomic relationship matrix
│   ├── GBLUP_results.csv                         # GBLUP cross-validation results
│   └── LightGBM_results.pkl                      # LightGBM cross-validation results
│
├── notebooks/
│   ├── SNPfiltering_for_AUSPAK_samples.sh       # VCF filtering pipeline
│   ├── Kinship_PCA.ipynb                        # Kinship matrix & PCA calculation
│   ├── BLUEs.ipynb                              # Phenotype processing (BLUEs)
│   ├── GBLUP.ipynb                              # GBLUP genomic prediction
│   ├── LightGBM.ipynb                           # Machine learning prediction
│   └── Visualizations.ipynb                     # Results visualization
```

## Requirements

### R Environment

```r
# Install required R packages
install.packages(c("dplyr", "Matrix"))
install.packages("asreml")        # Requires license
install.packages("AGHmatrix")
install.packages("ASRgenomics")
install.packages("data.table")
install.packages("SNPRelate")
```

### Python Environment

```python
# Create conda environment
conda create -n genomic_pred python=3.10
conda activate genomic_pred

# Install required packages
pip install pandas numpy scikit-learn lightgbm matplotlib seaborn
```

### External Tools

- **bcftools** (v1.10+) - For VCF file manipulation
- **GATK** (optional) - If starting from raw sequencing data, see https://github.com/queenoa/Quinoa-SNP-calling

## Data Description

### Phenotypic Data

**Traits analyzed:**
- **DTF** - Days to flowering (days)
- **DTH** - Days to harvest maturity (days)
- **PtHt** - Plant height (cm)
- **PcleLng** - Panicle length (cm)
- **SdLen** - Seed length (mm)
- **TGW** - Thousand grain weight (g)
- **SdW_z** - Seed weight per plant, z-transformed

**Environments:**
- Australia: 3 years (2017, 2018, 2019)
- Pakistan: 3 years (2019-20, 2020-21, 2021-22)

### Genotypic Data

- **Samples:** 551 quinoa accessions total
- **SNPs:** Filtered for quality 
- **Chromosomes:** 18 nuclear chromosomes (9 homeologous pairs: Cq1A-Cq9A, Cq1B-Cq9B)

## Workflow

### 1. SNP Filtering and Quality Control

Filter VCF file for genomic prediction analysis:

```bash
bash SNPfiltering_for_AUSPAK_samples.sh
```

**Quality filters applied:**
- Nuclear chromosomes only (excludes organellar genomes)
- Biallelic SNPs only
- MAF ≥ 0.01
- Missing data < 20% per variant
- Mean depth 5-30× (excludes low coverage and likely paralogs)
- Genotype calls only (GT field retained)

### 2. Kinship Matrix and PCA Calculation

Generate genomic relationship matrix and principal components:

```r
# Run Kinship_PCA.ipynb

# Creates:
# - kinship_matrix_VanRaden_auspak_maxmissing20.RData
# - quinoa_551_AUSPAK_PCs.csv
```

**Methods:**
- Kinship matrix: VanRaden method (AGHmatrix package)
- PCA: All principal components calculated (SNPRelate package)

### 3. Phenotype Processing - BLUEs Estimation

Calculate Best Linear Unbiased Estimates (BLUEs) for each trait:

```r
# Run BLUEs.ipynb
```

**Model specification:**

```r
# Fixed effects: accession + year
# Random effects: year:trial_replicate (PAK only)
# Residual: units
```

The analysis:
- Accounts for year and replicate effects
- Estimates generalized heritability
- Provides diagnostic plots for model assessment
- Filters predictions to observed year:accession combinations only

**Outputs:**
- `AUSPAK_phenotypes_means_BLUEs.csv` - BLUEs for all traits
- Heritability estimates
- Diagnostic plots (residuals, Q-Q plots)

### 4. GBLUP Genomic Prediction

Perform genomic prediction using GBLUP with ASReml-R:

```r
# Run GBLUP.ipynb
```

**Model specification:**

```r
# Fixed effects: location
# Random effects: vm(sample.id, Ginv_sparse) + location:year
# Residual: units
```

**Cross-validation strategy:**
- Group k-fold CV (k=5)
- Genotypes assigned to folds (prevents data leakage)
- 15 iterations with different random seeds
- Location-specific evaluation
- z-score transformation by environment

**Evaluation metrics:**
- Pearson correlation (prediction accuracy)
- NDCG@10 (ranking ability for top 10 genotypes within a fold, ~ 20% selection intensity)

### 5. Machine Learning Genomic Prediction

Implement gradient boosting approach using LightGBM:

```python
# Run LightGBM.ipynb
```

**Feature set:**
- All principal components (capturing population structure)
- One-hot encoded location and environment variables
- No year variable (not informative across locations)

**Model parameters:**
```python
model_params = {
    'max_depth': 3,
    'learning_rate': 0.05,
    'n_estimators': 500
}
```

**Cross-validation strategy:**
- Group k-fold CV (k=5, matching GBLUP)
- Genotype-based folding
- 15 iterations with different random seeds
- Location-specific evaluation
- z-score transformation by environment

**Evaluation metrics:**
- Pearson correlation
- NDCG@10

### 6. Visualization and Comparison

Generate plots comparing methods:

```python
# Run Visualizations.ipynb
```

## Key Methodological Details

### Z-score Transformation

Both GBLUP and LightGBM include optional z-score transformation by environment to account for environmental effects. Recommended, since we are not accounting for environmental factors and want to focus on relative performance of accessions in each environment:

```
z_ijl = (y_ijl - μ_l) / σ_l
```

where:
- `y_ijl` = phenotype of genotype i in replicate j at environment l
- `μ_l` = mean of environment l
- `σ_l` = standard deviation of environment l

### Matrix Bending (GBLUP only)

The genomic relationship matrix is "bent" to ensure positive definiteness:

1. Eigenvalue decomposition
2. Adjustment of eigenvalues below threshold (default: 1e-6)
3. Matrix reconstruction
4. Inversion to sparse format for computational efficiency

### Evaluation Strategy

**Group k-fold cross-validation:**
- Prevents data leakage by grouping all observations from the same genotype
- Maintains temporal structure (all years included in training/testing)
- Location-specific evaluation prevents artifical inflation of correlation values due to different population means

**Ranking metrics:**
- NDCG@10 evaluates ability to identify top 10 performers within each fold, ~ 20 % selection intensity
- Critical for breeding programs focused on selecting superior genotypes
- Accounts for trait directionality (e.g., lower DTF is better)

## Results Structure

### GBLUP Results (`GBLUP_results.csv`)

Columns:
- `iteration` - CV iteration number
- `fold` - Fold number
- `trait` - Trait name
- `location` - Testing location
- `pearson` - Pearson correlation
- `ndcg_at_10` - NDCG@10 score
- `seed` - Random seed used
- `n_test_genotypes` - Number of test genotypes

### LightGBM Results (`LightGBM_results.pkl`)

Same structure as GBLUP results (pandas DataFrame)

## Citation

If you use this code or data, please cite:



## License



## Contact

clara.stanschewski@kaust.edu.sa



---



## Version History

- **v1.0** (2025-11-10) - Initial release

