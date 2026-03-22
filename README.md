# hsCRP and Mortality Across PREVENT Risk Strata

**Additive Prognostic Value of High-Sensitivity C-Reactive Protein Across AHA PREVENT-Stratified Cardiovascular Risk Groups: A National Cohort Study**

Andrew Bouras, B.S. and Vikash Jaiswal, M.D.

## Overview

This repository contains the analysis code for our study examining the independent prognostic value of high-sensitivity C-reactive protein (hsCRP) across cardiovascular risk strata defined by the AHA PREVENT equation, using NHANES 2003-2016 with mortality follow-up through December 31, 2019.

## Data Source

All data come from the National Health and Nutrition Examination Survey (NHANES), publicly available from the CDC:
- Survey cycles: 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2015-2016
- Linked mortality files through 2019
- Note: hsCRP was not collected in 2011-2012 or 2013-2014

The analysis script automatically downloads all required NHANES files from the CDC website on first run and caches them locally in `data/nhanes_raw/`.

## Methods

1. **PREVENT risk scoring** -- 10-year CVD risk calculated using the AHA PREVENT equation (Khan et al., Circulation 2024), with coefficients validated against the `preventr` R package (v0.11.0)
2. **eGFR** -- CKD-EPI 2021 race-free creatinine equation (Inker et al., NEJM 2021)
3. **hsCRP categories** -- AHA/CDC classification: <1, 1-3, >3 mg/L (Pearson et al., Circulation 2003); also dichotomized at 2 mg/L per 2025 ESC/EAS focused update
4. **Survival analysis** -- Cox proportional hazards (all-cause mortality, cause-specific CV mortality), adjusted for age, sex, total cholesterol, HDL-C, SBP, BMI, eGFR, diabetes, smoking, and BP medications
5. **Discrimination** -- Harrell's C-statistic with 1,000-resample bootstrap

## Repository Structure

```
analysis/
  nhanes_crp_analysis.py    Main analysis script (single file, ~1,450 lines)
data/
  nhanes_raw/               Auto-downloaded NHANES XPT files (gitignored)
results/
  table1_demographics.csv   Baseline characteristics by PREVENT tier x CRP
  table2_cox_allcause.csv   Cox PH results, all-cause mortality
  table3_cox_cvmortality.csv   Cause-specific Cox PH, CV mortality
  table4_cox_crp2_*.csv     Dichotomized CRP (>=2 mg/L) models
  c_statistic_bootstrap.csv C-statistic with bootstrap 95% CI
  sample_sizes.csv          Sample sizes and event counts per stratum
  analytic_cohort.csv       Full analytic cohort (N=14,124)
figures/
  km_curves_allcause_by_crp_riskgroup.png
  forest_allcause_mortality.png
  forest_cv_mortality.png
  crp_distribution_by_tier.png
  mortality_rates_by_crp_tier.png
  prevent_risk_distribution.png
```

## How to Run

```bash
pip install -r requirements.txt
python analysis/nhanes_crp_analysis.py
```

The script will:
1. Download all NHANES data files from the CDC (~97 files, cached after first run)
2. Merge demographics, labs, questionnaires, and mortality data across 5 survey cycles
3. Apply inclusion/exclusion criteria (age 30-79, no baseline ASCVD, hsCRP <= 10 mg/L)
4. Calculate CKD-EPI 2021 eGFR and PREVENT 10-year CVD risk for each participant
5. Run all Cox PH models, generate tables and figures
6. Output everything to `results/` and `figures/`

## Key Findings

- Analytic cohort: N = 14,124 (mean follow-up 10.2 years; 1,364 all-cause deaths, 335 CV deaths)
- hsCRP >3 mg/L was independently associated with higher all-cause mortality in both low-risk (HR 1.39, 95% CI 1.06-1.83) and borderline-intermediate-risk (HR 1.43, 95% CI 1.18-1.74) PREVENT strata
- Dichotomized CRP >= 2 mg/L showed concordant results across American (AHA) and European (ESC/EAS) thresholds

## License

This analysis uses publicly available NHANES data. The code is provided for reproducibility of our published findings.
