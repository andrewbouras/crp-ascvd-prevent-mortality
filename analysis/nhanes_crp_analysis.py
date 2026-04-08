#!/usr/bin/env python3
"""
NHANES CRP and ASCVD Mortality Analysis
========================================
Additive Prognostic Value of hsCRP Across PREVENT-Stratified CV Risk Groups

Data:    NHANES 2003-2016 (primary), mortality follow-up through Dec 31, 2019
Methods: Cox PH (all-cause and CV mortality), cause-specific hazard approach
         Stratified by PREVENT 10-year CVD risk tiers × AHA hsCRP categories

Author:  Andrew Bouras, B.S. | Dr. Vikash Jaiswal (PI)
Date:    2026-02-26

References:
  PREVENT equation: Khan SS et al. Circulation 2024;149(6):430-449
  PREVENT R package: Mayer MG. preventr v0.11.0 (CRAN)
  CKD-EPI 2021: Inker LA et al. N Engl J Med 2021;385:1737-1749
  hsCRP thresholds: Pearson TA et al. Circulation 2003;107:499-511
  Template paper: Al-jarshawi MJ et al. Eur J Prev Cardiol 2026 (zwag037)
"""

import os, sys, warnings, re, io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import requests
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────
# SECTION 1: CONFIGURATION
# ─────────────────────────────────────────────────────────────

SCRIPT_DIR   = Path(__file__).parent
PROJECT_DIR  = SCRIPT_DIR.parent.parent
DATA_DIR     = PROJECT_DIR / 'data' / 'nhanes_raw'
RESULTS_DIR  = SCRIPT_DIR.parent / 'results'
FIGURES_DIR  = PROJECT_DIR / 'figures'

for d in [DATA_DIR, RESULTS_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

NHANES_BASE   = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
MORT_BASE     = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"

# Primary analysis cycles (hsCRP available and adequate follow-up through 2019).
# NOTE: hsCRP was NOT collected in 2011-2012 (G) or 2013-2014 (H).
PRIMARY_CYCLES = {
    'C': {'range': '2003-2004', 'y1': 2003, 'y2': 2004},
    'D': {'range': '2005-2006', 'y1': 2005, 'y2': 2006},
    'E': {'range': '2007-2008', 'y1': 2007, 'y2': 2008},
    'F': {'range': '2009-2010', 'y1': 2009, 'y2': 2010},
    'I': {'range': '2015-2016', 'y1': 2015, 'y2': 2016},
}
N_CYCLES = len(PRIMARY_CYCLES)

# Cycles where CRP is reported as LBXCRP in mg/dL (must ×10 to get mg/L)
CRP_MGDL_CYCLES = {'C', 'D', 'E', 'F'}

# Cycle-specific lab file name overrides (older NHANES naming)
# Key: (cycle_sfx, canonical_name) → actual_file_base
FILE_OVERRIDES = {
    ('C', 'TCHOL'): 'L13',       # Total chol + HDL in one file
    ('C', 'HDL'):   'L13',       # Same file
    ('C', 'CRP'):   'L11',       # CRP file
    ('C', 'BIOPRO'): 'L40',      # Biochemistry
    ('C', 'TRIGLY'): 'L13AM',    # LDL/TG
}

# ─────────────────────────────────────────────────────────────
# SECTION 2: DATA DOWNLOAD UTILITIES
# ─────────────────────────────────────────────────────────────

def download_xpt(sfx: str, file_base: str, yr_range: str, force: bool = False) -> object:
    """Download an NHANES XPT file using the CDC URL format, cache locally.
    Handles cycle-specific file name overrides via FILE_OVERRIDES dict.
    """
    # Apply cycle-specific overrides
    actual_base = FILE_OVERRIDES.get((sfx, file_base), file_base)
    fname_upper = f"{actual_base}_{sfx}.XPT"
    fname_lower = f"{actual_base}_{sfx}.xpt"
    local = DATA_DIR / fname_upper

    if local.exists() and local.stat().st_size > 500 and not force:
        try:
            return pd.read_sas(str(local), format='xport', encoding='iso-8859-1')
        except Exception:
            pass

    # CDC URL: /Nchs/Data/Nhanes/Public/{year}/DataFiles/{filename}
    year = yr_range.split('-')[0]
    for variant in [fname_upper, fname_lower]:
        url = f"{NHANES_BASE}/{year}/DataFiles/{variant}"
        try:
            r = requests.get(url, timeout=90, allow_redirects=True)
            # Validate it's actually an XPT file (not an HTML redirect)
            if r.status_code == 200 and r.content[:6] == b'HEADER':
                local.write_bytes(r.content)
                df = pd.read_sas(str(local), format='xport', encoding='iso-8859-1')
                print(f"  ✓ {fname_upper}")
                return df
        except Exception:
            continue
    print(f"  ✗ MISSING: {fname_upper}")
    return None


def download_mortality(y1: int, y2: int, force: bool = False) -> object:
    """Download and parse NHANES Public-Use Linked Mortality File (fixed-width)."""
    fname = f"NHANES_{y1}_{y2}_MORT_2019_PUBLIC.dat"
    local = DATA_DIR / fname
    if not local.exists() or force:
        url = f"{MORT_BASE}/{fname}"
        try:
            r = requests.get(url, timeout=90)
            if r.status_code == 200:
                local.write_bytes(r.content)
                print(f"  ✓ {fname}")
            else:
                print(f"  ✗ MISSING: {fname}")
                return None
        except Exception as e:
            print(f"  ✗ ERROR downloading {fname}: {e}")
            return None
    # Parse fixed-width format (per CDC documentation)
    colspecs = [(0,6),(14,15),(15,16),(16,19),(19,20),(20,21),(21,22),(22,26),
                (26,34),(34,42),(42,45),(45,48)]
    colnames = ['SEQN','ELIGSTAT','MORTSTAT','UCOD_LEADING','DIABETES_MORT',
                'HYPERTEN_MORT','DODQTR','DODYEAR','WGT_NEW','SA_WGT_NEW',
                'PERMTH_INT','PERMTH_EXM']
    try:
        df = pd.read_fwf(str(local), colspecs=colspecs, names=colnames,
                         dtype={'UCOD_LEADING': str})
        df['SEQN'] = pd.to_numeric(df['SEQN'], errors='coerce')
        df['MORTSTAT']     = pd.to_numeric(df['MORTSTAT'], errors='coerce')
        df['PERMTH_INT']   = pd.to_numeric(df['PERMTH_INT'], errors='coerce')
        df['UCOD_LEADING'] = df['UCOD_LEADING'].astype(str).str.strip()
        df['ELIGSTAT']     = pd.to_numeric(df['ELIGSTAT'], errors='coerce')
        return df
    except Exception as e:
        print(f"  ✗ Parse error {fname}: {e}")
        return None


# ─────────────────────────────────────────────────────────────
# SECTION 3: CLINICAL CALCULATIONS
# ─────────────────────────────────────────────────────────────

def ckd_epi_2021(scr: float, age: float, sex: int) -> float:
    """
    CKD-EPI 2021 race-free creatinine eGFR (Inker LA et al. NEJM 2021).
    sex: 1=male, 2=female
    Returns eGFR in mL/min/1.73m²
    """
    if pd.isna(scr) or pd.isna(age) or scr <= 0 or age <= 0:
        return np.nan
    if sex == 2:   # Female
        kappa, alpha = 0.7, -0.241
        sex_factor = 1.012
    else:           # Male
        kappa, alpha = 0.9, -0.302
        sex_factor = 1.0
    ratio = scr / kappa
    egfr = (142 * (min(ratio, 1.0) ** alpha) *
            (max(ratio, 1.0) ** (-1.200)) *
            (0.9938 ** age) * sex_factor)
    return egfr


def prevent_10yr_cvd(age, sex, tc, hdl, sbp, bp_rx, dm, smoker, bmi, egfr) -> float:
    """
    AHA PREVENT 10-year total CVD risk (base model, without UACR/HbA1c/SDI).
    Per Khan SS et al. Circulation 2024;149(6):430-449.

    Coefficients from the preventr R package (v0.11.0, CRAN), which implements
    the published PREVENT equations using sex-specific logistic regression with
    linearly centered and scaled predictors and piecewise linear splines.
    Validated against the AHA PREVENT online calculator.

    Inputs (all in standard clinical units):
      age   : years (30-79)
      sex   : 1=male, 2=female
      tc    : total cholesterol (mg/dL)
      hdl   : HDL-C (mg/dL)
      sbp   : systolic BP (mmHg)
      bp_rx : antihypertensive treatment (0/1)
      dm    : diabetes (treated or diagnosed) (0/1)
      smoker: current smoker (0/1)
      bmi   : kg/m²
      egfr  : CKD-EPI 2021 eGFR (mL/min/1.73m²)

    Returns 10-year CVD risk (0-1 scale)
    """
    vals = [age, tc, hdl, sbp, bmi, egfr]
    if any(pd.isna(v) for v in vals):
        return np.nan
    if any(v <= 0 for v in vals):
        return np.nan

    # Variable transformations (centering/scaling per preventr specification)
    CHOL_CONV = 0.02586  # mg/dL to mmol/L

    age_s       = (age - 55) / 10                              # centered at 55, per 10yr
    non_hdl_c   = (tc - hdl) * CHOL_CONV - 3.5                # non-HDL-C mmol/L, centered at 3.5
    hdl_s       = (hdl * CHOL_CONV - 1.3) / 0.3               # HDL mmol/L, centered at 1.3, per 0.3
    sbp_lt110   = (min(sbp, 110) - 110) / 20                  # SBP <110 mmHg, per 20
    sbp_gte110  = (max(sbp, 110) - 130) / 20                  # SBP >=110, centered at 130, per 20
    bmi_lt30    = (min(bmi, 30) - 25) / 5                     # BMI <30, per 5
    bmi_gte30   = (max(bmi, 30) - 30) / 5                     # BMI >=30, per 5
    egfr_lt60   = (min(egfr, 60) - 60) / (-15)                # eGFR <60, per -15
    egfr_gte60  = (max(egfr, 60) - 90) / (-15)                # eGFR >=60, centered at 90, per -15

    bp_rx  = float(bp_rx)
    dm     = float(dm)
    smoker = float(smoker)
    statin = 0.0  # statin status not available in NHANES questionnaire

    # Interaction terms
    bp_tx_sbp_gte110  = bp_rx * sbp_gte110
    statin_non_hdl_c  = statin * non_hdl_c
    age_non_hdl_c     = age_s * non_hdl_c
    age_hdl_c         = age_s * hdl_s
    age_sbp_gte110    = age_s * sbp_gte110
    age_dm            = age_s * dm
    age_smoking       = age_s * smoker
    age_bmi_gte30     = age_s * bmi_gte30
    age_egfr_lt60     = age_s * egfr_lt60

    if sex == 2:   # Women
        lp = (-3.3077280 +               # constant
               0.7939329 * age_s +
               0.0305239 * non_hdl_c +
              -0.1606857 * hdl_s +
              -0.2394003 * sbp_lt110 +
               0.3600781 * sbp_gte110 +
               0.8667604 * dm +
               0.5360739 * smoker +
               0.0000000 * bmi_lt30 +
               0.0000000 * bmi_gte30 +
               0.6045917 * egfr_lt60 +
               0.0433769 * egfr_gte60 +
               0.3151672 * bp_rx +
              -0.1477655 * statin +
              -0.0663612 * bp_tx_sbp_gte110 +
               0.1197879 * statin_non_hdl_c +
              -0.0819715 * age_non_hdl_c +
               0.0306769 * age_hdl_c +
              -0.0946348 * age_sbp_gte110 +
              -0.2705700 * age_dm +
              -0.0787150 * age_smoking +
               0.0000000 * age_bmi_gte30 +
              -0.1637806 * age_egfr_lt60)
    else:          # Men
        lp = (-3.0311680 +               # constant
               0.7688528 * age_s +
               0.0736174 * non_hdl_c +
              -0.0954431 * hdl_s +
              -0.4347345 * sbp_lt110 +
               0.3362658 * sbp_gte110 +
               0.7692857 * dm +
               0.4386871 * smoker +
               0.0000000 * bmi_lt30 +
               0.0000000 * bmi_gte30 +
               0.5378979 * egfr_lt60 +
               0.0164827 * egfr_gte60 +
               0.2888790 * bp_rx +
              -0.1337349 * statin +
              -0.0475924 * bp_tx_sbp_gte110 +
               0.1502730 * statin_non_hdl_c +
              -0.0517874 * age_non_hdl_c +
               0.0191169 * age_hdl_c +
              -0.1049477 * age_sbp_gte110 +
              -0.2251948 * age_dm +
              -0.0895067 * age_smoking +
               0.0000000 * age_bmi_gte30 +
              -0.1543702 * age_egfr_lt60)

    # Logistic link function: probability = exp(lp) / (1 + exp(lp))
    risk_10yr = 1.0 / (1.0 + np.exp(-lp))
    return float(np.clip(risk_10yr, 0.0001, 0.9999))


def validate_prevent() -> None:
    """
    Regression test: assert prevent_10yr_cvd matches preventr v0.11.0 (CRAN)
    on three reference patients spanning the risk continuum.

    Reference values were generated from preventr::estimate_risk() with
    model="base", time="10yr", and chol_unit="mg/dL". preventr internally
    rounds output to 3 decimal places, which is the agreement tolerance.

    This block is run unconditionally at the start of main() so the analysis
    will refuse to proceed if the PREVENT implementation has drifted.
    """
    cases = [
        # (label, kwargs, preventr_total_cvd_10yr)
        ("A: 40F, optimal labs, no risk factors",
         dict(age=40, sex=2, tc=180, hdl=60, sbp=110, bp_rx=0,
              dm=0, smoker=0, bmi=24, egfr=95),
         0.005),
        ("B: 60M, treated HTN, mild CKD (borderline-intermediate)",
         dict(age=60, sex=1, tc=220, hdl=42, sbp=140, bp_rx=1,
              dm=0, smoker=0, bmi=29, egfr=72),
         0.109),
        ("C: 70M, diabetes + smoker + CKD (high)",
         dict(age=70, sex=1, tc=240, hdl=38, sbp=158, bp_rx=1,
              dm=1, smoker=1, bmi=31, egfr=52),
         0.399),
    ]
    print("\n=== PREVENT validation against preventr v0.11.0 ===")
    all_pass = True
    for label, kwargs, expected in cases:
        actual = prevent_10yr_cvd(**kwargs)
        rounded = round(actual, 3)
        status = "PASS" if rounded == expected else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  [{status}] {label}")
        print(f"         Python = {actual:.6f} → {rounded:.3f}   preventr = {expected:.3f}")
    if not all_pass:
        raise RuntimeError(
            "PREVENT implementation does not match preventr v0.11.0. "
            "Refusing to run downstream analysis until coefficients are corrected."
        )
    print("  All 3 reference patients match preventr v0.11.0 to 3 decimal places.\n")


def prevent_risk_tier(risk: float) -> str:
    """Categorize PREVENT 10-yr CVD risk per AHA 2018/2019 guideline tiers."""
    if pd.isna(risk):
        return np.nan
    if risk < 0.05:
        return 'Low (<5%)'
    elif risk < 0.075:
        return 'Borderline (5-7.5%)'
    elif risk < 0.20:
        return 'Intermediate (7.5-20%)'
    else:
        return 'High (≥20%)'


def crp_category(crp: float) -> str:
    """AHA/CDC hsCRP categories (Pearson TA et al. Circulation 2003)."""
    if pd.isna(crp):
        return np.nan
    if crp < 1.0:
        return '<1 mg/L'
    elif crp <= 3.0:
        return '1-3 mg/L'
    else:
        return '>3 mg/L'


# ─────────────────────────────────────────────────────────────
# SECTION 4: LOAD ONE NHANES CYCLE
# ─────────────────────────────────────────────────────────────

def load_cycle(sfx: str, yr: dict) -> object:
    """Download, merge, and return a processed dataframe for one NHANES cycle."""
    yr_range = yr['range']
    print(f"\n--- Cycle {yr_range} (suffix _{sfx}) ---")

    # ── Demographics ──
    demo = download_xpt(sfx, 'DEMO', yr_range)
    if demo is None:
        return None
    demo_cols = ['SEQN','RIDAGEYR','RIAGENDR','WTMEC2YR','SDMVPSU','SDMVSTRA']
    # Add race variable (RIDRETH1 in early cycles, RIDRETH3 from 2011-2012+)
    for race_var in ['RIDRETH3','RIDRETH1']:
        if race_var in demo.columns:
            demo_cols.append(race_var)
            break
    demo = demo[[c for c in demo_cols if c in demo.columns]].copy()

    # ── Body Measures (BMI) ──
    bmx = download_xpt(sfx, 'BMX', yr_range)
    if bmx is not None and 'BMXBMI' in bmx.columns:
        demo = demo.merge(bmx[['SEQN','BMXBMI']].dropna(subset=['SEQN']),
                          on='SEQN', how='left')

    # ── Blood Pressure ──
    bpx = download_xpt(sfx, 'BPX', yr_range)
    if bpx is not None:
        bp_cols = [c for c in ['BPXSY2','BPXSY3','BPXSY4'] if c in bpx.columns]
        if bp_cols:
            bpx['SBP_MEAN'] = bpx[bp_cols].apply(
                lambda row: row.dropna().mean() if row.notna().sum() >= 1 else np.nan,
                axis=1)
            demo = demo.merge(bpx[['SEQN','SBP_MEAN']].dropna(subset=['SEQN']),
                              on='SEQN', how='left')

    # ── Total Cholesterol ──
    tchol = download_xpt(sfx, 'TCHOL', yr_range)
    if tchol is not None and 'LBXTC' in tchol.columns:
        demo = demo.merge(tchol[['SEQN','LBXTC']].dropna(subset=['SEQN']),
                          on='SEQN', how='left')

    # ── HDL Cholesterol ──
    # In 2003-2004 (C), HDL and TC are in the same file (L13_C, already loaded as TCHOL).
    # In later cycles, HDL is in HDL_{sfx}.xpt → LBDHDD.
    if sfx == 'C' and tchol is not None and 'LBXHDD' in tchol.columns:
        demo = demo.merge(tchol[['SEQN','LBXHDD']].rename(
            columns={'LBXHDD': 'HDL_C'}).dropna(subset=['SEQN']),
            on='SEQN', how='left')
    else:
        hdl_file = download_xpt(sfx, 'HDL', yr_range)
        if hdl_file is not None:
            for hdl_var in ['LBDHDD','LBXHDD','LBDHDL']:
                if hdl_var in hdl_file.columns:
                    demo = demo.merge(hdl_file[['SEQN', hdl_var]].rename(
                        columns={hdl_var: 'HDL_C'}).dropna(subset=['SEQN']),
                        on='SEQN', how='left')
                    break

    # ── hsCRP ──
    # 2003-2010: LBXCRP in CRP/L11 file (mg/dL) → multiply by 10 for mg/L
    # 2015-2016+: LBXHSCRP in HSCRP file (mg/L, no conversion needed)
    crp_loaded = False
    # Try primary CRP file
    crp_df = download_xpt(sfx, 'CRP', yr_range)
    if crp_df is not None:
        if 'LBXHSCRP' in crp_df.columns:
            demo = demo.merge(crp_df[['SEQN','LBXHSCRP']].dropna(subset=['SEQN']),
                              on='SEQN', how='left')
            crp_loaded = True
        elif 'LBXCRP' in crp_df.columns and sfx in CRP_MGDL_CYCLES:
            # Convert mg/dL → mg/L
            tmp = crp_df[['SEQN','LBXCRP']].copy()
            tmp['LBXHSCRP'] = tmp['LBXCRP'] * 10.0
            demo = demo.merge(tmp[['SEQN','LBXHSCRP']].dropna(subset=['SEQN']),
                              on='SEQN', how='left')
            crp_loaded = True
    if not crp_loaded:
        # Try HSCRP file (2015-2016+)
        crp_df2 = download_xpt(sfx, 'HSCRP', yr_range)
        if crp_df2 is not None and 'LBXHSCRP' in crp_df2.columns:
            demo = demo.merge(crp_df2[['SEQN','LBXHSCRP']].dropna(subset=['SEQN']),
                              on='SEQN', how='left')

    # ── Biochemistry (serum creatinine) ──
    bio = download_xpt(sfx, 'BIOPRO', yr_range)
    if bio is not None and 'LBXSCR' in bio.columns:
        demo = demo.merge(bio[['SEQN','LBXSCR']].dropna(subset=['SEQN']),
                          on='SEQN', how='left')

    # ── LDL Cholesterol (optional, for covariate adjustment) ──
    trigly = download_xpt(sfx, 'TRIGLY', yr_range)
    if trigly is not None:
        for ldl_var in ['LBDLDL','LBXLDL','LBDLDLSI']:
            if ldl_var in trigly.columns:
                demo = demo.merge(trigly[['SEQN', ldl_var]].rename(
                    columns={ldl_var: 'LDL_C'}).dropna(subset=['SEQN']),
                    on='SEQN', how='left')
                break

    # ── Diabetes ──
    diq = download_xpt(sfx, 'DIQ', yr_range)
    if diq is not None:
        diq_sub = diq[['SEQN']].copy()
        # Diabetes: told have diabetes OR taking insulin OR taking diabetic pills
        for col in ['DIQ010','DIQ050','DIQ070']:
            if col in diq.columns:
                diq_sub[col] = diq[col]
        demo = demo.merge(diq_sub.dropna(subset=['SEQN']), on='SEQN', how='left')

    # ── Smoking ──
    smq = download_xpt(sfx, 'SMQ', yr_range)
    if smq is not None:
        smq_sub = smq[['SEQN']].copy()
        for col in ['SMQ020','SMQ040']:
            if col in smq.columns:
                smq_sub[col] = smq[col]
        demo = demo.merge(smq_sub.dropna(subset=['SEQN']), on='SEQN', how='left')

    # ── Blood Pressure Medication ──
    bpq = download_xpt(sfx, 'BPQ', yr_range)
    if bpq is not None:
        bpq_sub = bpq[['SEQN']].copy()
        for col in ['BPQ040A','BPQ050A']:
            if col in bpq.columns:
                bpq_sub[col] = bpq[col]
        demo = demo.merge(bpq_sub.dropna(subset=['SEQN']), on='SEQN', how='left')

    # ── Medical Conditions (ASCVD history) ──
    mcq = download_xpt(sfx, 'MCQ', yr_range)
    if mcq is not None:
        mcq_cols = ['SEQN']
        for col in ['MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F']:
            if col in mcq.columns:
                mcq_cols.append(col)
        demo = demo.merge(mcq[mcq_cols].dropna(subset=['SEQN']), on='SEQN', how='left')

    # ── Mortality ──
    mort = download_mortality(yr['y1'], yr['y2'])
    if mort is not None:
        demo = demo.merge(mort[['SEQN','ELIGSTAT','MORTSTAT','UCOD_LEADING',
                                 'PERMTH_INT','PERMTH_EXM']].dropna(subset=['SEQN']),
                          on='SEQN', how='left')

    demo['CYCLE'] = yr_range
    return demo


# ─────────────────────────────────────────────────────────────
# SECTION 5: BUILD ANALYTIC COHORT
# ─────────────────────────────────────────────────────────────

def build_cohort(raw: pd.DataFrame) -> pd.DataFrame:
    """Apply inclusion/exclusion criteria and derive analytic variables."""
    df = raw.copy()
    n_start = len(df)
    print(f"\nCohort assembly: N = {n_start:,}")

    # ── Restrict to mortality-eligible adults (ELIGSTAT = 1) ──
    if 'ELIGSTAT' in df.columns:
        df = df[df['ELIGSTAT'] == 1]
    print(f"  After eligibility restriction: {len(df):,}")

    # ── Age 30-79 ──
    df = df[df['RIDAGEYR'].between(30, 79, inclusive='both')]
    print(f"  After age 30-79: {len(df):,}")

    # ── Exclude baseline ASCVD (CHF, CHD, angina, MI, stroke) ──
    ascvd_cols = ['MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F']
    for col in ascvd_cols:
        if col in df.columns:
            df = df[~(df[col] == 1)]
    print(f"  After excluding baseline ASCVD: {len(df):,}")

    # ── Require hsCRP ──
    df = df[df['LBXHSCRP'].notna()]
    print(f"  After requiring hsCRP: {len(df):,}")

    # ── Require PREVENT inputs ──
    required = ['BMXBMI','SBP_MEAN','LBXTC','HDL_C','LBXSCR']
    for col in required:
        if col in df.columns:
            df = df[df[col].notna()]
    print(f"  After requiring PREVENT inputs: {len(df):,}")

    # ── Require survey weight ──
    df = df[df['WTMEC2YR'].notna() & (df['WTMEC2YR'] > 0)]
    print(f"  After requiring survey weight: {len(df):,}")

    # ── Derive diabetes (treated or diagnosed) ──
    dm_cond = pd.Series(False, index=df.index)
    if 'DIQ010' in df.columns:
        dm_cond = dm_cond | (df['DIQ010'] == 1)
    if 'DIQ050' in df.columns:
        dm_cond = dm_cond | (df['DIQ050'] == 1)
    if 'DIQ070' in df.columns:
        dm_cond = dm_cond | (df['DIQ070'] == 1)
    df['DM'] = dm_cond.astype(int)

    # ── Derive current smoker ──
    smoker_cond = pd.Series(False, index=df.index)
    if 'SMQ020' in df.columns and 'SMQ040' in df.columns:
        smoker_cond = (df['SMQ020'] == 1) & (df['SMQ040'].isin([1, 2]))
    elif 'SMQ020' in df.columns:
        smoker_cond = (df['SMQ020'] == 1)
    df['SMOKER'] = smoker_cond.astype(int)

    # ── Derive BP medication ──
    bprx_cond = pd.Series(False, index=df.index)
    if 'BPQ050A' in df.columns:
        bprx_cond = (df['BPQ050A'] == 1)
    elif 'BPQ040A' in df.columns:
        bprx_cond = (df['BPQ040A'] == 1)
    df['BP_RX'] = bprx_cond.astype(int)

    # ── Calculate CKD-EPI 2021 eGFR ──
    df['EGFR'] = df.apply(
        lambda r: ckd_epi_2021(r['LBXSCR'], r['RIDAGEYR'], r['RIAGENDR']), axis=1)
    df = df[df['EGFR'].notna() & (df['EGFR'] > 0)]
    print(f"  After eGFR calculation: {len(df):,}")

    # ── Calculate PREVENT 10-year CVD risk ──
    df['PREVENT_10YR'] = df.apply(
        lambda r: prevent_10yr_cvd(
            age=r['RIDAGEYR'], sex=r['RIAGENDR'],
            tc=r['LBXTC'], hdl=r['HDL_C'],
            sbp=r['SBP_MEAN'], bp_rx=r['BP_RX'],
            dm=r['DM'], smoker=r['SMOKER'],
            bmi=r['BMXBMI'], egfr=r['EGFR']), axis=1)
    df = df[df['PREVENT_10YR'].notna()]
    print(f"  After PREVENT calculation: {len(df):,}")

    # ── Categorize PREVENT risk tier ──
    df['PREVENT_TIER'] = df['PREVENT_10YR'].apply(prevent_risk_tier)

    # ── Merge borderline + intermediate → borderline-intermediate ──
    df['PREVENT_TIER3'] = df['PREVENT_TIER'].replace({
        'Borderline (5-7.5%)': 'Borderline-Intermediate',
        'Intermediate (7.5-20%)': 'Borderline-Intermediate',
        'Low (<5%)': 'Low',
        'High (≥20%)': 'High',
    })
    tier_order = ['Low', 'Borderline-Intermediate', 'High']
    df['PREVENT_TIER3'] = pd.Categorical(df['PREVENT_TIER3'],
                                          categories=tier_order, ordered=True)

    # ── Categorize hsCRP ──
    # Primary analysis: exclude hsCRP >10 mg/L (acute inflammation)
    df_primary = df[df['LBXHSCRP'] <= 10.0].copy()
    df_primary['CRP_CAT'] = df_primary['LBXHSCRP'].apply(crp_category)
    crp_order = ['<1 mg/L', '1-3 mg/L', '>3 mg/L']
    df_primary['CRP_CAT'] = pd.Categorical(df_primary['CRP_CAT'],
                                            categories=crp_order, ordered=True)

    # hsCRP dichotomized at 2 mg/L (AHA/CDC threshold)
    df_primary['CRP_HIGH2'] = (df_primary['LBXHSCRP'] >= 2.0).astype(int)

    print(f"  After excluding CRP >10 mg/L: {len(df_primary):,}")

    # ── Survival variables ──
    # Follow-up time in months from examination (PERMTH_INT: interview date)
    # Use PERMTH_EXM (from exam date) when available; force numeric
    if 'PERMTH_EXM' in df_primary.columns:
        df_primary['FOLLOW_MONTHS'] = pd.to_numeric(df_primary['PERMTH_EXM'], errors='coerce')
    elif 'PERMTH_INT' in df_primary.columns:
        df_primary['FOLLOW_MONTHS'] = pd.to_numeric(df_primary['PERMTH_INT'], errors='coerce')
    else:
        df_primary['FOLLOW_MONTHS'] = np.nan
    df_primary['FOLLOW_YEARS'] = df_primary['FOLLOW_MONTHS'] / 12.0

    # Force mortality columns to numeric
    df_primary['MORTSTAT'] = pd.to_numeric(df_primary['MORTSTAT'], errors='coerce').fillna(0)
    df_primary['UCOD_LEADING'] = df_primary['UCOD_LEADING'].astype(str).str.strip().str.zfill(3)

    # All-cause mortality
    df_primary['ALL_CAUSE_DEATH'] = (df_primary['MORTSTAT'] == 1).astype(int)

    # CV mortality (heart disease [001] + cerebrovascular disease [005])
    df_primary['CV_DEATH'] = (
        (df_primary['MORTSTAT'] == 1) &
        (df_primary['UCOD_LEADING'].isin(['001', '005']))
    ).astype(int)

    # Non-CV death (competing event for Fine-Gray)
    df_primary['NONCV_DEATH'] = (
        (df_primary['MORTSTAT'] == 1) &
        (~df_primary['UCOD_LEADING'].isin(['001', '005']))
    ).astype(int)

    # Combined event indicator: 1=CV death, 2=non-CV death, 0=alive (for Fine-Gray)
    df_primary['CV_EVENT'] = 0
    df_primary.loc[df_primary['CV_DEATH'] == 1, 'CV_EVENT'] = 1
    df_primary.loc[df_primary['NONCV_DEATH'] == 1, 'CV_EVENT'] = 2

    # ── Combined survey weight (divided by number of cycles) ──
    df_primary['WGTMEC_COMBINED'] = df_primary['WTMEC2YR'] / N_CYCLES

    # ── Age group ──
    df_primary['AGE_GROUP'] = pd.cut(df_primary['RIDAGEYR'],
                                      bins=[29,44,59,79],
                                      labels=['30-44', '45-59', '60-79'])

    # ── Race/ethnicity ──
    race_col = None
    for rc in ['RIDRETH3','RIDRETH1']:
        if rc in df_primary.columns:
            race_col = rc
            break
    if race_col:
        race_map = {1: 'Mexican American', 2: 'Other Hispanic', 3: 'Non-Hispanic White',
                    4: 'Non-Hispanic Black', 5: 'Other/Multi', 6: 'Non-Hispanic Asian',
                    7: 'Other/Multi'}
        df_primary['RACE'] = df_primary[race_col].map(race_map)

    # Filter valid follow-up
    df_primary = df_primary[df_primary['FOLLOW_YEARS'] > 0]
    print(f"  Final analytic cohort: {len(df_primary):,}")

    return df_primary


# ─────────────────────────────────────────────────────────────
# SECTION 6: SURVIVAL ANALYSIS
# ─────────────────────────────────────────────────────────────

def cox_hr(df: pd.DataFrame, outcome_col: str, time_col: str = 'FOLLOW_YEARS',
           crp_ref: str = '<1 mg/L', covariates: list = None,
           weight_col = None) -> pd.DataFrame:
    """
    Cox proportional hazards regression for one outcome.
    Returns HR (95% CI) per CRP category vs reference.
    """
    from lifelines import CoxPHFitter

    results = []
    tiers = ['Low', 'Borderline-Intermediate', 'High']

    for tier in tiers:
        sub = df[df['PREVENT_TIER3'] == tier].copy()
        if len(sub) < 30:
            continue

        # Indicator encoding: CRP categories (ref = <1 mg/L)
        sub['CRP_MOD'] = (sub['CRP_CAT'] == '1-3 mg/L').astype(int)
        sub['CRP_HIGH'] = (sub['CRP_CAT'] == '>3 mg/L').astype(int)

        cols = [time_col, outcome_col, 'CRP_MOD', 'CRP_HIGH']
        if covariates:
            for c in covariates:
                if c in sub.columns:
                    cols.append(c)

        sub_model = sub[cols].dropna()

        if sub_model[outcome_col].sum() < 10:
            continue

        try:
            cph = CoxPHFitter(penalizer=0.01)
            fit_kwargs = {'duration_col': time_col, 'event_col': outcome_col}
            if weight_col and weight_col in sub.columns:
                sub_model[weight_col] = sub[weight_col]
                fit_kwargs['weights_col'] = weight_col
                fit_kwargs['robust'] = True
            cph.fit(sub_model, **fit_kwargs)

            for crp_label, crp_var in [('1-3 mg/L', 'CRP_MOD'), ('>3 mg/L', 'CRP_HIGH')]:
                if crp_var in cph.params_.index:
                    hr = np.exp(cph.params_[crp_var])
                    ci_lo = np.exp(cph.confidence_intervals_.loc[crp_var, '95% lower-bound'])
                    ci_hi = np.exp(cph.confidence_intervals_.loc[crp_var, '95% upper-bound'])
                    pval = cph.summary.loc[crp_var, 'p']
                    n_sub = len(sub_model)
                    n_events = int(sub_model[outcome_col].sum())
                    results.append({
                        'PREVENT_Tier': tier,
                        'CRP_Category': crp_label,
                        'HR': round(hr, 2),
                        'CI_Lower': round(ci_lo, 2),
                        'CI_Upper': round(ci_hi, 2),
                        'P_Value': '<0.0001' if pval < 0.0001 else f'{pval:.4f}',
                        'N': n_sub,
                        'Events': n_events,
                        'HR_CI': f"{hr:.2f} ({ci_lo:.2f}–{ci_hi:.2f})"
                    })
        except Exception as e:
            print(f"  Cox error in {tier}: {e}")

    return pd.DataFrame(results)


def cause_specific_cox(df: pd.DataFrame, outcome_col: str = 'CV_DEATH',
                       time_col: str = 'FOLLOW_YEARS',
                       covariates: list = None,
                       weight_col=None) -> pd.DataFrame:
    """
    Cause-specific Cox PH for CV mortality (standard competing-risks approach).
    Non-CV deaths are censored at time of death.
    Uses full covariate set matching the all-cause model.
    """
    from lifelines import CoxPHFitter

    results = []
    tiers = ['Low', 'Borderline-Intermediate', 'High']

    for tier in tiers:
        sub = df[df['PREVENT_TIER3'] == tier].copy()
        if len(sub) < 30:
            continue

        sub['CRP_MOD'] = (sub['CRP_CAT'] == '1-3 mg/L').astype(int)
        sub['CRP_HIGH'] = (sub['CRP_CAT'] == '>3 mg/L').astype(int)

        cols = [time_col, outcome_col, 'CRP_MOD', 'CRP_HIGH']
        if covariates:
            for c in covariates:
                if c in sub.columns:
                    cols.append(c)

        sub_model = sub[cols].dropna()

        if sub_model[outcome_col].sum() < 10:
            continue

        try:
            cph = CoxPHFitter(penalizer=0.01)
            fit_kwargs = {'duration_col': time_col, 'event_col': outcome_col}
            if weight_col and weight_col in sub.columns:
                sub_model[weight_col] = sub[weight_col]
                fit_kwargs['weights_col'] = weight_col
                fit_kwargs['robust'] = True
            cph.fit(sub_model, **fit_kwargs)

            for crp_label, crp_var in [('1-3 mg/L', 'CRP_MOD'), ('>3 mg/L', 'CRP_HIGH')]:
                if crp_var in cph.params_.index:
                    hr = np.exp(cph.params_[crp_var])
                    ci_lo = np.exp(cph.confidence_intervals_.loc[crp_var, '95% lower-bound'])
                    ci_hi = np.exp(cph.confidence_intervals_.loc[crp_var, '95% upper-bound'])
                    pval = cph.summary.loc[crp_var, 'p']
                    n_sub = len(sub_model)
                    n_events = int(sub_model[outcome_col].sum())
                    results.append({
                        'PREVENT_Tier': tier,
                        'CRP_Category': crp_label,
                        'HR': round(hr, 2),
                        'CI_Lower': round(ci_lo, 2),
                        'CI_Upper': round(ci_hi, 2),
                        'P_Value': '<0.0001' if pval < 0.0001 else f'{pval:.4f}',
                        'N': n_sub,
                        'Events': n_events,
                        'HR_CI': f"{hr:.2f} ({ci_lo:.2f}–{ci_hi:.2f})"
                    })
        except Exception as e:
            print(f"  Cause-specific Cox error in {tier}: {e}")

    return pd.DataFrame(results)


def cox_hr_dichotomized(df: pd.DataFrame, outcome_col: str,
                        time_col: str = 'FOLLOW_YEARS',
                        covariates: list = None) -> pd.DataFrame:
    """
    Cox PH with dichotomized CRP >=2 mg/L (AHA risk-enhancer threshold).
    Returns HR per PREVENT tier for CRP_HIGH2 vs <2 mg/L.
    """
    from lifelines import CoxPHFitter

    results = []
    tiers = ['Low', 'Borderline-Intermediate', 'High']

    for tier in tiers:
        sub = df[df['PREVENT_TIER3'] == tier].copy()
        if len(sub) < 30:
            continue

        cols = [time_col, outcome_col, 'CRP_HIGH2']
        if covariates:
            for c in covariates:
                if c in sub.columns:
                    cols.append(c)

        sub_model = sub[cols].dropna()
        if sub_model[outcome_col].sum() < 10:
            continue

        try:
            cph = CoxPHFitter(penalizer=0.01)
            cph.fit(sub_model, duration_col=time_col, event_col=outcome_col)

            if 'CRP_HIGH2' in cph.params_.index:
                hr = np.exp(cph.params_['CRP_HIGH2'])
                ci_lo = np.exp(cph.confidence_intervals_.loc['CRP_HIGH2', '95% lower-bound'])
                ci_hi = np.exp(cph.confidence_intervals_.loc['CRP_HIGH2', '95% upper-bound'])
                pval = cph.summary.loc['CRP_HIGH2', 'p']
                results.append({
                    'PREVENT_Tier': tier,
                    'Comparison': 'CRP >=2 vs <2 mg/L',
                    'HR': round(hr, 2),
                    'CI_Lower': round(ci_lo, 2),
                    'CI_Upper': round(ci_hi, 2),
                    'P_Value': '<0.0001' if pval < 0.0001 else f'{pval:.4f}',
                    'N': len(sub_model),
                    'Events': int(sub_model[outcome_col].sum()),
                    'HR_CI': f"{hr:.2f} ({ci_lo:.2f}–{ci_hi:.2f})"
                })
        except Exception as e:
            print(f"  CRP_HIGH2 Cox error in {tier}: {e}")

    return pd.DataFrame(results)


def cox_interaction_test(df: pd.DataFrame, outcome_col: str,
                         time_col: str = 'FOLLOW_YEARS',
                         covariates: list = None) -> dict:
    """
    Formal test of CRP × PREVENT-tier interaction.

    Fits unweighted Cox PH models on the full cohort:
      - Reduced: main effects for CRP category, PREVENT tier, and covariates
      - Full:    reduced + four CRP × tier product terms
    Reports a likelihood-ratio test (4 df) for the joint null that all four
    interaction coefficients equal zero, plus a confirmatory Wald χ² and the
    individual interaction coefficient estimates.

    Survey-weighted Cox models give numerically unstable joint Wald statistics
    in lifelines (variance matrix scaling artifacts), so the interaction test
    is reported on the unweighted model. The reference HR estimates and primary
    effect estimates remain survey-weighted in the main tables.
    """
    from lifelines import CoxPHFitter
    from scipy import stats
    import warnings

    sub = df.copy()
    sub['CRP_MOD'] = (sub['CRP_CAT'] == '1-3 mg/L').astype(int)
    sub['CRP_HIGH'] = (sub['CRP_CAT'] == '>3 mg/L').astype(int)
    sub['TIER_BI'] = (sub['PREVENT_TIER3'] == 'Borderline-Intermediate').astype(int)
    sub['TIER_HIGH'] = (sub['PREVENT_TIER3'] == 'High').astype(int)
    sub['CRP_MOD_x_BI'] = sub['CRP_MOD'] * sub['TIER_BI']
    sub['CRP_MOD_x_HIGH'] = sub['CRP_MOD'] * sub['TIER_HIGH']
    sub['CRP_HIGH_x_BI'] = sub['CRP_HIGH'] * sub['TIER_BI']
    sub['CRP_HIGH_x_HIGH'] = sub['CRP_HIGH'] * sub['TIER_HIGH']

    interaction_terms = ['CRP_MOD_x_BI', 'CRP_MOD_x_HIGH',
                         'CRP_HIGH_x_BI', 'CRP_HIGH_x_HIGH']
    main_terms = ['CRP_MOD', 'CRP_HIGH', 'TIER_BI', 'TIER_HIGH']

    cov_cols = [c for c in (covariates or []) if c in sub.columns]
    cols_full = [time_col, outcome_col] + main_terms + interaction_terms + cov_cols
    cols_red = [time_col, outcome_col] + main_terms + cov_cols

    m_full = sub[cols_full].dropna()
    m_red = sub[cols_red].dropna()
    n = len(m_full)
    n_events = int(m_full[outcome_col].sum())

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cph_full = CoxPHFitter(penalizer=0.01)
        cph_full.fit(m_full, duration_col=time_col, event_col=outcome_col,
                     show_progress=False)
        cph_red = CoxPHFitter(penalizer=0.01)
        cph_red.fit(m_red, duration_col=time_col, event_col=outcome_col,
                    show_progress=False)

    LL_full = cph_full.log_likelihood_
    LL_red = cph_red.log_likelihood_
    df_test = len(interaction_terms)
    lrt = float(2 * (LL_full - LL_red))
    lrt_p = float(1 - stats.chi2.cdf(lrt, df=df_test)) if lrt > 0 else float('nan')

    # Confirmatory Wald χ² on the same unweighted fit
    try:
        V = cph_full.variance_matrix_.loc[interaction_terms, interaction_terms].values
        beta = cph_full.params_.loc[interaction_terms].values
        wald_chi2 = float(beta @ np.linalg.inv(V) @ beta)
        wald_p = float(1 - stats.chi2.cdf(wald_chi2, df=df_test))
    except (np.linalg.LinAlgError, KeyError):
        wald_chi2 = float('nan')
        wald_p = float('nan')

    coef_rows = []
    for term in interaction_terms + main_terms:
        if term in cph_full.params_.index:
            hr = float(np.exp(cph_full.params_[term]))
            ci_lo = float(np.exp(cph_full.confidence_intervals_.loc[term, '95% lower-bound']))
            ci_hi = float(np.exp(cph_full.confidence_intervals_.loc[term, '95% upper-bound']))
            pval = float(cph_full.summary.loc[term, 'p'])
            coef_rows.append({
                'Term': term,
                'HR': round(hr, 3),
                'CI_Lower': round(ci_lo, 3),
                'CI_Upper': round(ci_hi, 3),
                'P_Value': pval,
            })

    return {
        'outcome': outcome_col,
        'n': n,
        'events': n_events,
        'lrt': lrt,
        'lrt_df': df_test,
        'lrt_p': lrt_p,
        'wald_chi2': wald_chi2,
        'wald_df': df_test,
        'wald_p': wald_p,
        'coefficients': pd.DataFrame(coef_rows),
    }


def compute_c_statistic(df: pd.DataFrame, covars: list, n_bootstrap: int = 1000,
                        weight_col: str = None) -> list:
    """
    Compute Harrell's C-statistic for all-cause mortality models with and without CRP.
    Bootstrap 95% CI with optimism correction.
    """
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index

    np.random.seed(42)
    results = []
    tiers = ['Low', 'Borderline-Intermediate', 'High']

    use_weights = bool(weight_col)

    def _fit_kwargs(frame_cols):
        kw = {'duration_col': 'FOLLOW_YEARS', 'event_col': 'ALL_CAUSE_DEATH'}
        if use_weights and weight_col in frame_cols:
            kw['weights_col'] = weight_col
            kw['robust'] = True
        return kw

    for tier in tiers:
        sub = df[df['PREVENT_TIER3'] == tier].copy()
        sub['CRP_MOD'] = (sub['CRP_CAT'] == '1-3 mg/L').astype(int)
        sub['CRP_HIGH'] = (sub['CRP_CAT'] == '>3 mg/L').astype(int)

        cols_base = ['FOLLOW_YEARS', 'ALL_CAUSE_DEATH'] + [c for c in covars if c in sub.columns]
        cols_crp = cols_base + ['CRP_MOD', 'CRP_HIGH']
        if use_weights and weight_col in sub.columns:
            cols_base = cols_base + [weight_col]
            cols_crp = cols_crp + [weight_col]

        sub_base = sub[cols_base].dropna()
        sub_crp = sub[cols_crp].dropna()

        if sub_base['ALL_CAUSE_DEATH'].sum() < 20:
            continue

        try:
            # Fit base model (covariates only) — weighted if requested
            cph_base = CoxPHFitter(penalizer=0.01)
            cph_base.fit(sub_base, **_fit_kwargs(sub_base.columns))
            base_features = [c for c in sub_base.columns
                             if c not in ('FOLLOW_YEARS', 'ALL_CAUSE_DEATH', weight_col or '')]
            c_base = concordance_index(sub_base['FOLLOW_YEARS'],
                                       -cph_base.predict_partial_hazard(sub_base[base_features]),
                                       sub_base['ALL_CAUSE_DEATH'])

            # Fit CRP model (covariates + CRP)
            cph_crp = CoxPHFitter(penalizer=0.01)
            cph_crp.fit(sub_crp, **_fit_kwargs(sub_crp.columns))
            crp_features = [c for c in sub_crp.columns
                            if c not in ('FOLLOW_YEARS', 'ALL_CAUSE_DEATH', weight_col or '')]
            c_crp = concordance_index(sub_crp['FOLLOW_YEARS'],
                                      -cph_crp.predict_partial_hazard(sub_crp[crp_features]),
                                      sub_crp['ALL_CAUSE_DEATH'])

            # Bootstrap for 95% CI of delta-C — fits also weighted
            delta_boot = []
            n = len(sub_crp)
            for _ in range(n_bootstrap):
                idx = np.random.choice(n, n, replace=True)
                boot = sub_crp.iloc[idx]
                if boot['ALL_CAUSE_DEATH'].sum() < 5:
                    continue
                try:
                    cb = CoxPHFitter(penalizer=0.01)
                    cb.fit(boot[cols_base], **_fit_kwargs(cols_base))
                    c_b = concordance_index(boot['FOLLOW_YEARS'],
                                            -cb.predict_partial_hazard(boot[base_features]),
                                            boot['ALL_CAUSE_DEATH'])
                    cc = CoxPHFitter(penalizer=0.01)
                    cc.fit(boot[cols_crp], **_fit_kwargs(cols_crp))
                    c_c = concordance_index(boot['FOLLOW_YEARS'],
                                            -cc.predict_partial_hazard(boot[crp_features]),
                                            boot['ALL_CAUSE_DEATH'])
                    delta_boot.append(c_c - c_b)
                except Exception:
                    continue

            delta_c = c_crp - c_base
            if delta_boot:
                delta_ci_lo = np.percentile(delta_boot, 2.5)
                delta_ci_hi = np.percentile(delta_boot, 97.5)
            else:
                delta_ci_lo = delta_ci_hi = np.nan

            results.append({
                'PREVENT_Tier': tier,
                'C_base': round(c_base, 3),
                'C_with_CRP': round(c_crp, 3),
                'Delta_C': round(delta_c, 4),
                'Delta_CI_Lower': round(delta_ci_lo, 4) if not np.isnan(delta_ci_lo) else '',
                'Delta_CI_Upper': round(delta_ci_hi, 4) if not np.isnan(delta_ci_hi) else '',
                'N': n,
                'Events': int(sub_crp['ALL_CAUSE_DEATH'].sum())
            })
            print(f"  {tier}: C-base={c_base:.3f}, C+CRP={c_crp:.3f}, Delta={delta_c:.4f} "
                  f"({delta_ci_lo:.4f}–{delta_ci_hi:.4f})")
        except Exception as e:
            print(f"  C-statistic error in {tier}: {e}")

    return results


# ─────────────────────────────────────────────────────────────
# SECTION 7: TABLE GENERATION
# ─────────────────────────────────────────────────────────────

def generate_table1(df: pd.DataFrame) -> pd.DataFrame:
    """Survey-weighted descriptive statistics stratified by PREVENT tier × CRP category."""
    tiers = ['Low', 'Borderline-Intermediate', 'High']
    crp_cats = ['<1 mg/L', '1-3 mg/L', '>3 mg/L']

    rows = []
    def wt_mean_sd(series, weights):
        w = weights[series.notna()]
        v = series[series.notna()]
        if len(v) == 0:
            return 'N/A'
        wm = np.average(v, weights=w)
        ws = np.sqrt(np.average((v - wm)**2, weights=w))
        return f"{wm:.1f} ± {ws:.1f}"

    def wt_pct(mask, weights):
        return f"{100 * (weights[mask].sum() / weights.sum()):.1f}%"

    for tier in tiers:
        sub_t = df[df['PREVENT_TIER3'] == tier]
        for crp in crp_cats:
            sub = sub_t[sub_t['CRP_CAT'] == crp]
            w = sub['WGTMEC_COMBINED']
            n = len(sub)
            n_deaths = int(sub['ALL_CAUSE_DEATH'].sum())
            n_cv = int(sub['CV_DEATH'].sum())
            follow = wt_mean_sd(sub['FOLLOW_YEARS'], w) if n > 0 else 'N/A'
            age = wt_mean_sd(sub['RIDAGEYR'], w) if n > 0 else 'N/A'
            pct_male = wt_pct(sub['RIAGENDR'] == 1, w) if n > 0 else 'N/A'
            tc = wt_mean_sd(sub['LBXTC'], w) if n > 0 else 'N/A'
            hdl = wt_mean_sd(sub['HDL_C'], w) if n > 0 else 'N/A'
            sbp = wt_mean_sd(sub['SBP_MEAN'], w) if n > 0 else 'N/A'
            bmi = wt_mean_sd(sub['BMXBMI'], w) if n > 0 else 'N/A'
            egfr = wt_mean_sd(sub['EGFR'], w) if n > 0 else 'N/A'
            crp_med = (f"{sub['LBXHSCRP'].median():.2f}" if n > 0 else 'N/A')
            pct_dm = wt_pct(sub['DM'] == 1, w) if n > 0 else 'N/A'
            pct_smk = wt_pct(sub['SMOKER'] == 1, w) if n > 0 else 'N/A'
            pct_bprx = wt_pct(sub['BP_RX'] == 1, w) if n > 0 else 'N/A'
            rows.append({
                'PREVENT Tier': tier, 'hsCRP Category': crp, 'N': n,
                'All-Cause Deaths': n_deaths, 'CV Deaths': n_cv,
                'Follow-up (yrs)': follow, 'Age (yrs)': age,
                'Male (%)': pct_male, 'TC (mg/dL)': tc,
                'HDL-C (mg/dL)': hdl, 'SBP (mmHg)': sbp,
                'BMI (kg/m²)': bmi, 'eGFR (mL/min)': egfr,
                'hsCRP median (mg/L)': crp_med,
                'Diabetes (%)': pct_dm, 'Smoker (%)': pct_smk,
                'BP Meds (%)': pct_bprx,
            })

    return pd.DataFrame(rows)


def generate_sample_sizes(df: pd.DataFrame) -> pd.DataFrame:
    """N, events, mortality rate by stratum."""
    tiers = ['Low', 'Borderline-Intermediate', 'High']
    crp_cats = ['<1 mg/L', '1-3 mg/L', '>3 mg/L']
    rows = []
    for tier in tiers:
        for crp in crp_cats:
            sub = df[(df['PREVENT_TIER3'] == tier) & (df['CRP_CAT'] == crp)]
            n = len(sub)
            n_ac = int(sub['ALL_CAUSE_DEATH'].sum())
            n_cv = int(sub['CV_DEATH'].sum())
            rate_ac = f"{100*n_ac/n:.1f}%" if n > 0 else 'N/A'
            rate_cv = f"{100*n_cv/n:.1f}%" if n > 0 else 'N/A'
            follow = sub['FOLLOW_YEARS'].mean()
            rows.append({
                'PREVENT Tier': tier, 'hsCRP Category': crp,
                'N': n, 'All-Cause Deaths': n_ac, 'AC Mortality Rate': rate_ac,
                'CV Deaths': n_cv, 'CV Mortality Rate': rate_cv,
                'Mean Follow-up (yrs)': f"{follow:.1f}" if n > 0 else 'N/A'
            })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────
# SECTION 8: FIGURES
# ─────────────────────────────────────────────────────────────

TIER_COLORS = {
    'Low': '#2196F3',
    'Borderline-Intermediate': '#FF9800',
    'High': '#F44336',
}
CRP_LINESTYLES = {
    '<1 mg/L': '-',
    '1-3 mg/L': '--',
    '>3 mg/L': ':',
}
CRP_COLORS = {
    '<1 mg/L': '#4CAF50',
    '1-3 mg/L': '#FF9800',
    '>3 mg/L': '#F44336',
}


def plot_km_curves(df: pd.DataFrame) -> None:
    """Kaplan-Meier all-cause survival curves stratified by CRP within each PREVENT tier."""
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    tiers = ['Low', 'Borderline-Intermediate', 'High']
    crp_cats = ['<1 mg/L', '1-3 mg/L', '>3 mg/L']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('All-Cause Survival by hsCRP Category Within PREVENT Risk Strata\n'
                 'NHANES 2003–2016 (Primary Prevention Cohort)',
                 fontsize=12, y=1.02)

    for ax, tier in zip(axes, tiers):
        sub_t = df[df['PREVENT_TIER3'] == tier]
        kmf_list = []
        for crp in crp_cats:
            sub = sub_t[sub_t['CRP_CAT'] == crp]
            if len(sub) < 10:
                continue
            kmf = KaplanMeierFitter()
            kmf.fit(sub['FOLLOW_YEARS'], event_observed=sub['ALL_CAUSE_DEATH'],
                    label=crp)
            kmf.plot_survival_function(ax=ax,
                                       color=CRP_COLORS.get(crp, 'gray'),
                                       linestyle=CRP_LINESTYLES.get(crp, '-'),
                                       linewidth=1.5, ci_show=True, ci_alpha=0.1)
            kmf_list.append((sub, crp))

        # Log-rank test (low vs high)
        low = sub_t[sub_t['CRP_CAT'] == '<1 mg/L']
        high = sub_t[sub_t['CRP_CAT'] == '>3 mg/L']
        if len(low) > 5 and len(high) > 5:
            lr = logrank_test(low['FOLLOW_YEARS'], high['FOLLOW_YEARS'],
                              low['ALL_CAUSE_DEATH'], high['ALL_CAUSE_DEATH'])
            pstr = (f"p = {lr.p_value:.3f}" if lr.p_value >= 0.001
                    else "p < 0.001")
            ax.text(0.05, 0.05, f"Log-rank {pstr}", transform=ax.transAxes,
                    fontsize=8, style='italic')

        ax.set_title(f'{tier}', fontsize=10, fontweight='bold')
        ax.set_xlabel('Follow-up (years)')
        ax.set_ylabel('Survival probability')
        ax.set_xlim(0, 20)
        ax.set_ylim(0.5, 1.01)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = FIGURES_DIR / 'km_curves_allcause_by_crp_riskgroup.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_forest(cox_results: pd.DataFrame, title: str, outname: str) -> None:
    """Forest plot: HR (95% CI) by CRP category × PREVENT tier."""
    if cox_results.empty:
        print(f"  Skipping forest plot (no results): {outname}")
        return

    fig, ax = plt.subplots(figsize=(9, max(4, 0.8 * len(cox_results) + 1)))

    tiers = ['Low', 'Borderline-Intermediate', 'High']
    crp_cats = ['1-3 mg/L', '>3 mg/L']

    y_labels = []
    y_pos = []
    y = 0

    for tier in reversed(tiers):
        for crp in crp_cats:
            row = cox_results[
                (cox_results['PREVENT_Tier'] == tier) &
                (cox_results['CRP_Category'] == crp)
            ]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            hr_col = 'HR' if 'HR' in cox_results.columns else 'sHR'
            hr = row[hr_col]
            ci_lo = row['CI_Lower']
            ci_hi = row['CI_Upper']
            color = TIER_COLORS.get(tier, 'gray')

            ax.errorbar(hr, y, xerr=[[hr - ci_lo], [ci_hi - hr]],
                        fmt='s', color=color, ecolor=color,
                        capsize=3, markersize=7, linewidth=1.5)

            y_labels.append(f"{tier}\nhsCRP {crp}")
            y_pos.append(y)
            y += 1

        y += 0.5  # gap between tiers

    ax.axvline(1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=10)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.grid(True, axis='x', alpha=0.3)
    ax.set_xlim(left=0.3)

    # Legend
    patches = [mpatches.Patch(color=c, label=t) for t, c in TIER_COLORS.items()]
    ax.legend(handles=patches, title='PREVENT Tier', fontsize=8,
              loc='lower right')

    plt.tight_layout()
    outpath = FIGURES_DIR / outname
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_crp_distribution(df: pd.DataFrame) -> None:
    """hsCRP distribution by PREVENT tier."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    tiers = ['Low', 'Borderline-Intermediate', 'High']

    for ax, tier in zip(axes, tiers):
        sub = df[df['PREVENT_TIER3'] == tier]
        crp_vals = sub['LBXHSCRP'].clip(upper=10)
        ax.hist(crp_vals, bins=40, color=TIER_COLORS[tier], alpha=0.7, edgecolor='white')
        ax.axvline(1, color='gray', linestyle='--', linewidth=1, label='1 mg/L')
        ax.axvline(3, color='black', linestyle='--', linewidth=1, label='3 mg/L')
        ax.set_title(tier, fontsize=9, fontweight='bold')
        ax.set_xlabel('hsCRP (mg/L)')
        ax.set_ylabel('Count')
        if ax == axes[0]:
            ax.legend(fontsize=8)

    fig.suptitle('hsCRP Distribution by PREVENT Risk Tier (capped at 10 mg/L)',
                 fontsize=10)
    plt.tight_layout()
    outpath = FIGURES_DIR / 'crp_distribution_by_tier.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_mortality_rates(df: pd.DataFrame) -> None:
    """Bar chart: mortality rates by PREVENT tier × CRP category."""
    tiers = ['Low', 'Borderline-Intermediate', 'High']
    crp_cats = ['<1 mg/L', '1-3 mg/L', '>3 mg/L']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, (outcome, label) in zip(axes, [('ALL_CAUSE_DEATH', 'All-Cause Mortality'),
                                            ('CV_DEATH', 'CV Mortality')]):
        x = np.arange(len(crp_cats))
        width = 0.25
        for i, tier in enumerate(tiers):
            rates = []
            for crp in crp_cats:
                sub = df[(df['PREVENT_TIER3'] == tier) & (df['CRP_CAT'] == crp)]
                if len(sub) > 0:
                    rate = 100 * sub[outcome].sum() / len(sub)
                else:
                    rate = 0
                rates.append(rate)
            bars = ax.bar(x + i * width, rates, width,
                          label=tier, color=TIER_COLORS[tier], alpha=0.8)

        ax.set_xticks(x + width)
        ax.set_xticklabels(crp_cats, fontsize=9)
        ax.set_ylabel('Mortality Rate (%)')
        ax.set_title(f'{label} by hsCRP Category\nand PREVENT Risk Tier', fontsize=10)
        ax.legend(title='PREVENT Tier', fontsize=8)
        ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    outpath = FIGURES_DIR / 'mortality_rates_by_crp_tier.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_prevent_distribution(df: pd.DataFrame) -> None:
    """PREVENT 10-year CVD risk distribution by tier."""
    fig, ax = plt.subplots(figsize=(9, 5))
    tiers = ['Low', 'Borderline-Intermediate', 'High']
    for tier in tiers:
        sub = df[df['PREVENT_TIER3'] == tier]['PREVENT_10YR'] * 100
        ax.hist(sub, bins=40, alpha=0.6, label=tier, color=TIER_COLORS[tier],
                edgecolor='white')
    ax.axvline(5, color='gray', linestyle='--', linewidth=1)
    ax.axvline(20, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel('PREVENT 10-year CVD Risk (%)')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of PREVENT 10-year CVD Risk by Tier')
    ax.legend()
    plt.tight_layout()
    outpath = FIGURES_DIR / 'prevent_risk_distribution.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


# ─────────────────────────────────────────────────────────────
# SECTION 9: MAIN PIPELINE
# ─────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("NHANES CRP and ASCVD Mortality Analysis")
    print("Project: crp-ascvd-mortality")
    print("=" * 60)

    # ── Step 0: PREVENT validation against preventr v0.11.0 ──
    validate_prevent()

    # ── Step 1: Download and stack all cycles ──
    print("\n[Step 1] Downloading NHANES data (2003–2016)...")
    frames = []
    for sfx, yr in PRIMARY_CYCLES.items():
        frame = load_cycle(sfx, yr)
        if frame is not None:
            frames.append(frame)

    if not frames:
        print("ERROR: No data loaded. Check internet connection and NHANES URLs.")
        sys.exit(1)

    raw = pd.concat(frames, ignore_index=True)
    print(f"\nTotal raw records across all cycles: {len(raw):,}")

    # ── Step 2: Build analytic cohort ──
    print("\n[Step 2] Applying inclusion/exclusion criteria...")
    cohort = build_cohort(raw)

    if len(cohort) < 1000:
        print(f"WARNING: Very small cohort (N={len(cohort)}). Check data loading.")

    # Save cohort to CSV for audit
    cohort_path = RESULTS_DIR / 'analytic_cohort.csv'
    cohort.to_csv(cohort_path, index=False)
    print(f"  Analytic cohort saved: {cohort_path}")

    # ── Step 3: Sample sizes and descriptive stats ──
    print("\n[Step 3] Generating descriptive statistics...")

    # Sample sizes
    ss = generate_sample_sizes(cohort)
    ss.to_csv(RESULTS_DIR / 'sample_sizes.csv', index=False)
    print("  Saved: results/sample_sizes.csv")
    print(ss[['PREVENT Tier','hsCRP Category','N','All-Cause Deaths','CV Deaths']].to_string(index=False))

    # Table 1
    t1 = generate_table1(cohort)
    t1.to_csv(RESULTS_DIR / 'table1_demographics.csv', index=False)
    print("  Saved: results/table1_demographics.csv")

    # Overall cohort summary
    print(f"\n  Overall cohort summary:")
    print(f"    N total: {len(cohort):,}")
    print(f"    All-cause deaths: {cohort['ALL_CAUSE_DEATH'].sum():,} "
          f"({100*cohort['ALL_CAUSE_DEATH'].mean():.1f}%)")
    print(f"    CV deaths: {cohort['CV_DEATH'].sum():,} "
          f"({100*cohort['CV_DEATH'].mean():.1f}%)")
    print(f"    Mean follow-up: {cohort['FOLLOW_YEARS'].mean():.1f} years")
    print(f"    PREVENT tier distribution:")
    print(cohort['PREVENT_TIER3'].value_counts().to_string())
    print(f"    CRP category distribution:")
    print(cohort['CRP_CAT'].value_counts().to_string())

    # ── Step 4: Cox regression (all-cause mortality) ──
    print("\n[Step 4] Cox PH regression - all-cause mortality...")

    covars = ['RIDAGEYR', 'RIAGENDR', 'LBXTC', 'HDL_C', 'SBP_MEAN',
              'BMXBMI', 'EGFR', 'DM', 'SMOKER', 'BP_RX']

    # Unadjusted
    cox_unadj = cox_hr(cohort, 'ALL_CAUSE_DEATH', covariates=None)
    cox_unadj['Model'] = 'Unadjusted'

    # Adjusted
    cox_adj = cox_hr(cohort, 'ALL_CAUSE_DEATH', covariates=covars)
    cox_adj['Model'] = 'Adjusted'

    # Weighted adjusted
    cox_adj_wt = cox_hr(cohort, 'ALL_CAUSE_DEATH', covariates=covars,
                        weight_col='WGTMEC_COMBINED')
    cox_adj_wt['Model'] = 'Adjusted (Weighted)'

    cox_all = pd.concat([cox_unadj, cox_adj, cox_adj_wt], ignore_index=True)
    cox_all.to_csv(RESULTS_DIR / 'table2_cox_allcause.csv', index=False)
    print("  Saved: results/table2_cox_allcause.csv")
    print(cox_adj[['PREVENT_Tier','CRP_Category','HR_CI','P_Value','Events']].to_string(index=False))

    # ── Step 5: Cause-specific Cox PH (CV mortality) ──
    print("\n[Step 5] Cause-specific Cox PH - CV mortality...")
    cox_cv_adj = cause_specific_cox(cohort, 'CV_DEATH', covariates=covars)
    cox_cv_adj['Model'] = 'Cause-specific (Adjusted)'
    cox_cv_adj.to_csv(RESULTS_DIR / 'table3_cox_cvmortality.csv', index=False)
    print("  Saved: results/table3_cox_cvmortality.csv")
    if not cox_cv_adj.empty:
        print(cox_cv_adj[['PREVENT_Tier','CRP_Category','HR_CI','P_Value','Events']].to_string(index=False))

    # ── Step 5b: Formal CRP × PREVENT-tier interaction test ──
    print("\n[Step 5b] Formal CRP × PREVENT-tier interaction test (unweighted LRT)...")
    interaction_rows = []
    for outcome_label, outcome_col in [('All-cause', 'ALL_CAUSE_DEATH'),
                                       ('CV', 'CV_DEATH')]:
        try:
            res = cox_interaction_test(cohort, outcome_col, covariates=covars)
            print(f"  {outcome_label} mortality (N={res['n']:,}, events={res['events']}):")
            print(f"    LRT  χ²({res['lrt_df']}) = {res['lrt']:.3f}, "
                  f"p = {res['lrt_p']:.4f}")
            print(f"    Wald χ²({res['wald_df']}) = {res['wald_chi2']:.3f}, "
                  f"p = {res['wald_p']:.4f}")
            res['coefficients'].to_csv(
                RESULTS_DIR / f"interaction_coefs_{outcome_col.lower()}.csv",
                index=False)
            interaction_rows.append({
                'Outcome': outcome_label,
                'N': res['n'],
                'Events': res['events'],
                'LRT_chi2': round(res['lrt'], 4),
                'LRT_df': res['lrt_df'],
                'LRT_P': res['lrt_p'],
                'Wald_chi2': round(res['wald_chi2'], 4),
                'Wald_P': res['wald_p'],
            })
        except Exception as e:
            print(f"  Interaction test failed for {outcome_label}: {e}")
    if interaction_rows:
        pd.DataFrame(interaction_rows).to_csv(
            RESULTS_DIR / 'cox_interaction_test.csv', index=False)
        print("  Saved: results/cox_interaction_test.csv")

    # ── Step 6: Dichotomized CRP (>=2 mg/L) Cox models ──
    print("\n[Step 6] Dichotomized CRP (>=2 mg/L) Cox PH - all-cause mortality...")
    cox_crp2 = cox_hr_dichotomized(cohort, 'ALL_CAUSE_DEATH', covariates=covars)
    if not cox_crp2.empty:
        cox_crp2.to_csv(RESULTS_DIR / 'table4_cox_crp2_allcause.csv', index=False)
        print("  Saved: results/table4_cox_crp2_allcause.csv")
        print(cox_crp2.to_string(index=False))
    cox_crp2_cv = cox_hr_dichotomized(cohort, 'CV_DEATH', covariates=covars)
    if not cox_crp2_cv.empty:
        cox_crp2_cv.to_csv(RESULTS_DIR / 'table4_cox_crp2_cvmortality.csv', index=False)
        print("  Saved: results/table4_cox_crp2_cvmortality.csv")

    # ── Step 7: C-statistic (discrimination), survey-weighted ──
    print("\n[Step 7] C-statistic with bootstrap (survey-weighted Cox fits)...")
    cstat_results = compute_c_statistic(cohort, covars,
                                        n_bootstrap=300,
                                        weight_col='WGTMEC_COMBINED')
    if cstat_results:
        cstat_df = pd.DataFrame(cstat_results)
        cstat_df.to_csv(RESULTS_DIR / 'c_statistic_bootstrap.csv', index=False)
        print("  Saved: results/c_statistic_bootstrap.csv")
        print(cstat_df.to_string(index=False))

    # ── Step 8: Sensitivity analysis (include CRP >10 mg/L) ──
    print("\n[Step 8] Sensitivity analysis (without CRP >10 exclusion)...")
    raw_frames_full = []
    for sfx, yr in PRIMARY_CYCLES.items():
        frame = load_cycle(sfx, yr)  # Will use cached files
        if frame is not None:
            raw_frames_full.append(frame)
    if raw_frames_full:
        raw_full = pd.concat(raw_frames_full, ignore_index=True)
        cohort_full = build_cohort_all_crp(raw_full)
        if cohort_full is not None and len(cohort_full) > 0:
            cox_sens = cox_hr(cohort_full, 'ALL_CAUSE_DEATH', covariates=covars)
            cox_sens['Model'] = 'Sensitivity (all CRP)'
            cox_sens.to_csv(RESULTS_DIR / 'sensitivity_cox_all_crp.csv', index=False)
            print(f"  Sensitivity cohort N={len(cohort_full):,}")
            print("  Saved: results/sensitivity_cox_all_crp.csv")

    # ── Step 9: Figures ──
    print("\n[Step 9] Generating figures...")
    plot_km_curves(cohort)
    plot_forest(cox_adj, 'Adjusted HR (95% CI) for All-Cause Mortality\nhsCRP vs <1 mg/L Reference',
                'forest_allcause_mortality.png')
    if not cox_cv_adj.empty:
        plot_forest(cox_cv_adj,
                    'Cause-specific HR (95% CI) for CV Mortality\nhsCRP vs <1 mg/L Reference',
                    'forest_cv_mortality.png')
    plot_crp_distribution(cohort)
    plot_mortality_rates(cohort)
    plot_prevent_distribution(cohort)

    # ── Step 10: Summary statistics for manuscript ──
    print("\n[Step 10] Key findings summary...")
    summary_lines = []
    summary_lines.append("# Key Analysis Results\n")
    summary_lines.append(f"**Analytic cohort N:** {len(cohort):,}")
    summary_lines.append(f"**All-cause deaths:** {cohort['ALL_CAUSE_DEATH'].sum():,} "
                         f"({100*cohort['ALL_CAUSE_DEATH'].mean():.1f}%)")
    summary_lines.append(f"**CV deaths:** {cohort['CV_DEATH'].sum():,} "
                         f"({100*cohort['CV_DEATH'].mean():.1f}%)")
    summary_lines.append(f"**Mean follow-up:** {cohort['FOLLOW_YEARS'].mean():.1f} years")
    summary_lines.append(f"\n**PREVENT tier distribution:**")
    for tier, cnt in cohort['PREVENT_TIER3'].value_counts().items():
        summary_lines.append(f"  - {tier}: {cnt:,} ({100*cnt/len(cohort):.1f}%)")
    summary_lines.append(f"\n**hsCRP distribution:**")
    for cat, cnt in cohort['CRP_CAT'].value_counts().items():
        summary_lines.append(f"  - {cat}: {cnt:,} ({100*cnt/len(cohort):.1f}%)")

    if not cox_adj.empty:
        summary_lines.append(f"\n**Cox PH - Adjusted HRs (all-cause mortality):**")
        for _, row in cox_adj.iterrows():
            summary_lines.append(
                f"  - {row['PREVENT_Tier']} / {row['CRP_Category']}: "
                f"HR {row['HR_CI']} (p={row['P_Value']}, N events={row['Events']})")

    summary_path = RESULTS_DIR / 'key_findings_summary.md'
    with open(summary_path, 'w') as f:
        f.write('\n'.join(summary_lines))
    print(f"  Saved: {summary_path}")

    print("\n" + "=" * 60)
    print("Analysis complete.")
    print(f"Results: {RESULTS_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    print("=" * 60)


def build_cohort_all_crp(raw: pd.DataFrame) -> object:
    """Build cohort WITHOUT the CRP >10 mg/L exclusion (sensitivity analysis)."""
    df = raw.copy()
    if 'ELIGSTAT' in df.columns:
        df = df[df['ELIGSTAT'] == 1]
    df = df[df['RIDAGEYR'].between(30, 79, inclusive='both')]
    ascvd_cols = ['MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F']
    for col in ascvd_cols:
        if col in df.columns:
            df = df[~(df[col] == 1)]
    df = df[df['LBXHSCRP'].notna()]
    required = ['BMXBMI','SBP_MEAN','LBXTC','HDL_C','LBXSCR']
    for col in required:
        if col in df.columns:
            df = df[df[col].notna()]
    df = df[df['WTMEC2YR'].notna() & (df['WTMEC2YR'] > 0)]

    dm_cond = pd.Series(False, index=df.index)
    if 'DIQ010' in df.columns: dm_cond = dm_cond | (df['DIQ010'] == 1)
    if 'DIQ050' in df.columns: dm_cond = dm_cond | (df['DIQ050'] == 1)
    if 'DIQ070' in df.columns: dm_cond = dm_cond | (df['DIQ070'] == 1)
    df['DM'] = dm_cond.astype(int)

    smoker_cond = pd.Series(False, index=df.index)
    if 'SMQ020' in df.columns and 'SMQ040' in df.columns:
        smoker_cond = (df['SMQ020'] == 1) & (df['SMQ040'].isin([1, 2]))
    df['SMOKER'] = smoker_cond.astype(int)

    bprx_cond = pd.Series(False, index=df.index)
    if 'BPQ050A' in df.columns: bprx_cond = (df['BPQ050A'] == 1)
    elif 'BPQ040A' in df.columns: bprx_cond = (df['BPQ040A'] == 1)
    df['BP_RX'] = bprx_cond.astype(int)

    df['EGFR'] = df.apply(lambda r: ckd_epi_2021(r['LBXSCR'], r['RIDAGEYR'], r['RIAGENDR']), axis=1)
    df = df[df['EGFR'].notna() & (df['EGFR'] > 0)]

    df['PREVENT_10YR'] = df.apply(lambda r: prevent_10yr_cvd(
        r['RIDAGEYR'], r['RIAGENDR'], r['LBXTC'], r['HDL_C'],
        r['SBP_MEAN'], r['BP_RX'], r['DM'], r['SMOKER'],
        r['BMXBMI'], r['EGFR']), axis=1)
    df = df[df['PREVENT_10YR'].notna()]
    df['PREVENT_TIER'] = df['PREVENT_10YR'].apply(prevent_risk_tier)
    df['PREVENT_TIER3'] = df['PREVENT_TIER'].replace({
        'Borderline (5-7.5%)': 'Borderline-Intermediate',
        'Intermediate (7.5-20%)': 'Borderline-Intermediate',
        'Low (<5%)': 'Low', 'High (≥20%)': 'High'})
    tier_order = ['Low', 'Borderline-Intermediate', 'High']
    df['PREVENT_TIER3'] = pd.Categorical(df['PREVENT_TIER3'], categories=tier_order, ordered=True)

    df['CRP_CAT'] = df['LBXHSCRP'].apply(crp_category)
    crp_order = ['<1 mg/L', '1-3 mg/L', '>3 mg/L']
    df['CRP_CAT'] = pd.Categorical(df['CRP_CAT'], categories=crp_order, ordered=True)
    df['WGTMEC_COMBINED'] = df['WTMEC2YR'] / N_CYCLES

    if 'PERMTH_EXM' in df.columns:
        df['FOLLOW_MONTHS'] = pd.to_numeric(df['PERMTH_EXM'], errors='coerce')
    elif 'PERMTH_INT' in df.columns:
        df['FOLLOW_MONTHS'] = pd.to_numeric(df['PERMTH_INT'], errors='coerce')
    else:
        df['FOLLOW_MONTHS'] = np.nan
    df['FOLLOW_YEARS'] = df['FOLLOW_MONTHS'] / 12.0
    df['MORTSTAT'] = pd.to_numeric(df['MORTSTAT'], errors='coerce').fillna(0)
    df['UCOD_LEADING'] = df['UCOD_LEADING'].astype(str).str.strip().str.zfill(3)
    df['ALL_CAUSE_DEATH'] = (df['MORTSTAT'] == 1).astype(int)
    df['CV_DEATH'] = ((df['MORTSTAT'] == 1) & (df['UCOD_LEADING'].isin(['001','005']))).astype(int)
    df['NONCV_DEATH'] = ((df['MORTSTAT'] == 1) & (~df['UCOD_LEADING'].isin(['001','005']))).astype(int)
    df['CV_EVENT'] = 0
    df.loc[df['CV_DEATH'] == 1, 'CV_EVENT'] = 1
    df.loc[df['NONCV_DEATH'] == 1, 'CV_EVENT'] = 2
    df = df[df['FOLLOW_YEARS'] > 0]
    return df


if __name__ == '__main__':
    main()
