# Statistical Association Analysis

This folder contains the statistical analyses performed on quality-controlled genotyping data
generated using the Axiom™ Precision Medicine Diversity Array (PMDA).

All input datasets used in this folder are obtained from the preprocessing pipeline
(`01_preprocessing`) and are assumed to have already passed technical and genetic quality control.

---

## Overview

The analyses implemented in this folder focus on **genome-wide allelic association testing**
between cases and controls using **Fisher’s exact test**.

The aim is to identify genetic variants associated with disease status in a robust and
conservative statistical framework.

---

## Analysis Workflow

### 1. Input Data

- Genotype matrices encoded as 0/1/2 (additive model)
- Phenotype labels (case/control)
- Missing genotypes are retained and handled on a per-variant basis during testing

No genotype imputation is performed at this stage.

---

### 2. Fisher’s Exact Test

- Allelic association testing is performed for each variant using Fisher’s exact test
- A 2×2 contingency table is constructed based on allele counts in cases and controls
- Variants with insufficient information (e.g. monomorphic variants or zero minor allele count)
  are excluded from testing

---

### 3. Multiple Testing Correction

- P-values are adjusted for multiple testing using the **Benjamini–Hochberg False Discovery Rate (FDR)**
- A conservative **Bonferroni correction** is optionally applied post hoc to Fisher-significant variants
  for additional stringency

---

### 4. Post-Association Quality Checks

For variants passing statistical significance thresholds, additional checks are performed:

- Hardy–Weinberg equilibrium (HWE) tested separately in cases and controls
- Removal of variants showing significant deviation from HWE in controls

---

## Output Files

The main outputs generated in this folder include:

- Tables containing genome-wide Fisher test results
- Lists of statistically significant variants after multiple testing correction
- Filtered result tables suitable for downstream biological interpretation

---

## Notes

- This folder is dedicated exclusively to **statistical association testing**
- Biological annotation, pathway analysis, linkage disequilibrium analyses, and data visualization
  are performed in separate folders
- No machine learning models are implemented in this stage

`

