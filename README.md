# BXD Genetic Association Analysis

## Overview

This project investigates the genetic basis of phenotypic variation in BXD recombinant inbred mouse strains.

The goal is to analyze the association between genetic variation (SNP markers) and a selected phenotype using statistical models and compare their significance.

This study applies methods commonly used in **statistical genetics** and **genome-wide association studies (GWAS)**.

---

# Dataset

The analysis uses two datasets:

- **Phenotypes file**: quantitative trait measurements across strains
- **Genotypes file**: SNP markers across the genome

Genotype encoding:
- `H` → Heterozygous
- `U` → Unknown
- Other values represent homozygous states

---

# Methodology

A single phenotype and a SNP with at least one heterozygous strain were selected.

Three statistical models were applied:

---

## 1. Linear Regression (Ignoring Heterozygotes)

A simple linear regression model where heterozygous markers are excluded.

This model estimates the association between SNP genotype and phenotype.

---

## 2. Linear Regression (Including Heterozygotes)

A regression model where heterozygous genotypes are explicitly included as a third category.

This allows evaluation of potential dosage effects.

---

## 3. ANOVA Model (Ignoring Heterozygotes)

A one-way ANOVA model comparing phenotype means across genotype groups.

This model tests whether genotype groups significantly differ.

---

# Genome-Wide Association Scan

A full genome-wide scan was performed using linear regression (ignoring heterozygotes) for all SNPs.

For each SNP:
- A statistical test was performed
- P-values were computed
- Results were transformed into **-log10(P-value)** scores

---

# Manhattan Plot

A Manhattan plot was generated:

- X-axis: genomic position of SNPs
- Y-axis: -log10(P-value)

This visualization highlights genomic regions strongly associated with the phenotype.

---

# Key Output

- Identification of the most significant SNP
- Statistical significance evaluation
- Comparison between different genetic models

---

# Results Interpretation

The results demonstrate how different modeling assumptions (especially treatment of heterozygous markers) affect statistical significance and interpretation of genotype–phenotype associations.

---

# Technologies Used

- R / Python (depending on implementation)
- Linear regression (lm)
- ANOVA (aov)
- Statistical genetics methods
- Data visualization (Manhattan plot)

---

# Key Concepts

- SNP association analysis
- Genotype–phenotype modeling
- Linear regression
- ANOVA
- Genome-wide association study (GWAS)
- -log10(p-value) transformation

---

# Conclusion

This project demonstrates how statistical models can be used to uncover genetic associations underlying phenotypic variation in recombinant inbred mouse populations.

The Manhattan plot provides a genome-wide view of association signals, highlighting potentially meaningful genetic loci.

## Machine Learning Perspective

Although this project is rooted in statistical genetics, the core methodology relies on supervised learning techniques, primarily:

- Linear regression models for association learning
- Feature–outcome mapping (SNP → phenotype)
- Statistical hypothesis testing as model evaluation
- ANOVA as a multi-class learning framework

From an ML perspective, this project can be viewed as a predictive modeling task where genetic markers serve as features and phenotypic traits as targets.
