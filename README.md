# Survival Analysis of Heart Failure Patients

## Project Overview

This project analyzes survival outcomes in heart failure patients using the [Heart Failure Clinical Records Dataset](https://www.kaggle.com/datasets/rabieelkharoua/predict-survival-of-patients-with-heart-failure) from Kaggle. The dataset contains 299 patients from Pakistan with a follow-up period of up to 9 months.

**Team Members:** Cameron Chesbrough, Jing Lyu, Meitong Zhou, Zhengkun Ou, Xuange Liang

## Project Objectives

1. Identify significant predictors of mortality using adjusted hazard ratios
2. Compare survival across subgroups using Kaplan–Meier curves and log-rank tests
3. Build Cox proportional hazards models and assess model performance

## Dataset Description

- **Sample Size:** 299 heart failure patients
- **Outcome:** Time to death (right-censored for patients alive at end of follow-up)
- **Event Rate:** 32% died during follow-up, 68% alive/censored

### Key Variables

| Variable | Description |
|----------|-------------|
| `time` | Follow-up time (days) |
| `DEATH_EVENT` | Death indicator (1=died, 0=censored) |
| `age` | Age (years) |
| `ejection_fraction` | % of blood leaving heart at each contraction |
| `serum_creatinine` | Serum creatinine level (mg/dL) |
| `serum_sodium` | Serum sodium level (mEq/L) |
| `high_blood_pressure` | Hypertension (1=yes, 0=no) |
| `anaemia` | Anaemia (1=yes, 0=no) |
| `diabetes` | Diabetes (1=yes, 0=no) |
| `smoking` | Smoking status (1=yes, 0=no) |
| `sex` | Sex (1=male, 0=female) |

## Analysis Pipeline

### 1. Descriptive Analysis (`Descriptive.Rmd`) - Author: Jing Lyu

- Overall cohort summary statistics
- Comparison of characteristics by survival status (Wilcoxon/Chi-square tests)
- Distribution visualization (boxplots, histograms)
- Univariate Kaplan-Meier curves with log-rank tests

**Key Findings:**
- Patients who died were significantly older (median 65 vs 60 years, p<0.001)
- Lower ejection fraction in deceased (median 30% vs 38%, p<0.001)
- Higher serum creatinine in deceased (median 1.30 vs 1.00 mg/dL, p<0.001)
- Lower serum sodium in deceased (median 136 vs 137 mEq/L, p<0.001)
- No significant differences for sex, anaemia, diabetes, smoking

### 2. Cox Regression & Log-rank Tests (`analysis_1.Rmd`) - Author: Zhengkun Ou

- Overall Kaplan-Meier survival curve
- KM curves stratified by clinical covariates (high blood pressure, sex, smoking, diabetes)
- Univariate Cox models for binary predictors
- Multivariable Cox proportional hazards model
- Peto-Peto weighted log-rank test

**Key Findings:**
- High blood pressure shows significant survival difference (log-rank p=0.036)
- Multivariable Cox model significant predictors:
  - Age: HR per year increase
  - Ejection fraction: HR per % increase (protective)
  - Serum creatinine: HR per unit increase

### 3. PH Assumption & Interaction Effects (`analysis_2.Rmd`) - Author: Xuange Liang

- **Proportional Hazards Testing:**
  - Schoenfeld residuals test (`cox.zph()`)
  - Log-log survival plots
  - Stratified Cox model for PH violations

- **Addressing PH Violation with Interaction Model:**
  - Age × Ejection Fraction interaction (to address EF's PH violation)

**Key Findings:**
- Ejection fraction marginally violates PH assumption (p=0.039)
- **Age × EF interaction is statistically significant (p=0.015)**
  - The significant interaction explains the non-proportionality of EF
  - Protective effect of higher EF is more pronounced in older patients
  - In patients >70 years, low EF (≤30%) is associated with dramatically worse survival

### 4. Final Model Selection

Based on diagnostic analyses, we recommend the **Model with Age × EF Interaction**:

```r
cox_final <- coxph(
  Surv(time, DEATH_EVENT) ~ age * ejection_fraction + sex + 
    serum_creatinine + high_blood_pressure + diabetes + smoking,
  data = heart
)
```

**Rationale:**
1. The Age × EF interaction is statistically significant (p=0.015)
2. The interaction explains the marginal PH violation of EF (p=0.039)
3. Provides clinically meaningful insights: the protective effect of higher EF varies by age
4. Addresses non-proportionality of EF by allowing its effect to vary with age

## Model Comparison

| Model | AIC | C-index |
|-------|-----|---------|
| Base Model (No Interaction) | - | ~0.73 |
| Model with Age × EF Interaction | Lower | ~0.74 |
| Stratified by EF Group | - | ~0.72 |

## File Structure

```
Survival-Project-Heart-Failure/
├── README.md                              # Project documentation
├── heart_failure_clinical_records_dataset.csv  # Raw data
├── Descriptive.Rmd                        # EDA and descriptive analysis
├── Descriptive.html                       # Rendered EDA report
├── analysis_1.Rmd                         # Cox models and log-rank tests
├── analysis_1.html                        # Rendered analysis report
├── analysis_2.Rmd                         # PH assumption and interactions
├── analysis_2.html                        # Rendered analysis report
└── files/                                 # Additional project files
    └── Survival Analysis Project Proposal (Group 3).md
```

## Requirements

```r
# Required R packages
install.packages(c(
  "survival",
  "survminer", 
  "dplyr",
  "ggplot2",
  "gtsummary",
  "tableone",
  "scales",
  "tidyr"
))
```

## Key Clinical Implications

1. **Age and Ejection Fraction interact**: The protective effect of higher EF is more pronounced in older patients
2. **Strongest predictors of mortality**: Age, ejection fraction, and serum creatinine
3. **Serum creatinine** reflects renal function and independently predicts cardiovascular mortality
4. **Traditional risk factors** (smoking, diabetes) showed weaker associations in this cohort

## References

- Chicco, D., Jurman, G. Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone. *BMC Med Inform Decis Mak* 20, 16 (2020).

## Timeline

| Phase | Task | Deadline |
|-------|------|----------|
| Week 1 | EDA and Begin Modelling | Nov 22 |
| Week 2 | Complete Modelling and Interpretation | Nov 29 |
| Week 3 | Complete Interpretation and Slides | Dec 5 (Presentation) |
| Week 4 | Finalize Report | Dec 12 (Report Due) |

## License

This project is for educational purposes as part of a Survival Analysis course project