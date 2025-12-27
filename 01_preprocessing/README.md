# 🧬 Genotyping Quality Control Pipeline (Axiom PMDA)

This repository contains a reproducible preprocessing and quality control (QC) pipeline for genome-wide genotyping data generated using the **Axiom™ Precision Medicine Diversity Array (PMDA)**.

The pipeline is designed to obtain a high-quality set of genetic variants suitable for **statistical association analyses** in neurodegenerative disease cohorts.

---

## 📥 Input Data

- Genotyping data exported from **Axiom Analysis Suite**
- File format: tab-delimited `.txt`
- Initial number of variants: ~800,000
- The first 5 rows contain metadata and are skipped during import (`skiprows=5`)

---

## 🔎 Quality Control Workflow

### 1. Technical Quality Filtering

Variants are filtered according to Axiom technical QC metrics:

- `BestandRecommended == 1`
- `ConversionType == "PolyHighResolution"`

The following variants are excluded:
- Monomorphic variants
- Variants with low-quality or missing conversion types

**Output:**  
A set of high-quality polymorphic variants.

---

### 2. Genetic Quality Filtering

Standard genetic QC filters are applied:

- **Minor Allele Frequency (MAF ≥ 0.01)**
- **Hardy–Weinberg Equilibrium (HWE p-value ≥ 0.05)**

Variants with missing MAF or HWE values are excluded prior to filtering.

---

### 3. Functional Subset Selection

Variants are optionally restricted to biologically relevant functional modules (`super_module`), including:

- GWAS Grid  
- Immunity, Inflammation, and HLA  
- Published GWAS Hits  
- Pharmacogenetics / ADME  
- Wellness and Lifestyle  

This step allows focusing the analysis on variants of potential biological and clinical relevance.

---

### 4. Genotype Conversion (A/B → Real Alleles)

Axiom genotypes are converted from A/B encoding to real allelic representation using the annotation fields `Allele_A`, `Allele_B`, and `Ref_Allele`.

Output format examples: `AA`, `AG`, `GG`.

This representation enables biological interpretation of genotypes.

---

### 5. Missingness Assessment and Filtering

- Missing genotype rates are calculated per variant
- Variants with **>5% missing genotypes** are excluded
- Residual missing values are retained

**No genotype imputation is performed.**

---

### 6. Numerical Genotype Encoding (0 / 1 / 2)

Genotypes are encoded using an additive genetic model:

| Genotype | Encoding |
|--------|----------|
| Homozygous reference | 0 |
| Heterozygous | 1 |
| Homozygous alternative | 2 |
| Missing | NaN |

This representation is intended for downstream **statistical association analyses**.

---

## 📊 Statistical Analysis

This pipeline prepares data exclusively for **statistical testing** (e.g. Fisher’s exact test, chi-square test, logistic regression).

- Missing genotypes are handled using a **complete-case approach on a per-variant basis**
- No machine learning models are applied
- No genotype imputation is performed

---

## 💾 Output Files

| File | Description |
|----|------------|
| `*_stats_QC.txt` | Variants passing all QC filters |
| `*_QC_allelic.txt` | Genotypes in real allelic format |
| `*_QC_encoded.pkl` | Numerical genotype encoding (0/1/2) |

Files are also saved in `.pkl` format for efficient loading in **Pandas**.

---

## 📦 Summary

Starting from ~800,000 variants, this pipeline produces:
- A high-quality set of polymorphic variants
- Genotype matrices suitable for robust statistical analyses
- A transparent and reproducible QC workflow aligned with standard GWAS practices

---

## 📌 Notes

This repository focuses on **data preprocessing and quality control**.  
Downstream statistical analyses are performed in separate notebooks.


