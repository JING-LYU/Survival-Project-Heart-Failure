## Final Project Proposal 

Cameron Chesbrough, Jing Lyu, Meitong Zhou, Ou Zhengkun, Xuange Liang

**Topic: Predicting Survival Outcomes in Heart Failure Patients Using Clinical Risk Factors**

### **Objective**

### Aim: Modelling and predicting time-to-death among patients with heart failure.

1. Identify significant predictors using adjusted hazard ratios.  
2. Compare survival across subgroups using Kaplan–Meier curves and log-rank tests.  
3. Build Cox and Random Survival Forest models, assess them using C-index, KM separation, and calibration, and select the better-performing model.

### **Data Description**

This [dataset](https://www.kaggle.com/datasets/rabieelkharoua/predict-survival-of-patients-with-heart-failure) describes heart failure patients in Pakistan over a time period of 9 months. It includes demographic and clinical information about the patients, their time until death/censoring, and their survival status. 299 patients are included in the dataset with 11 continuous and categorical variables, as well as time to event and survival status. No information is missing.

**Statistical Analysis Plan**

1. **Exploratory Data Analysis**  
   * Descriptive statistics and visualization by survival status.  
   * Kaplan–Meier survival curves by selected covariates (e.g., sex, anemia, diabetes).  
   * Log-rank tests to compare survival distributions.  
2. **Modeling Approaches**  
   * **Univariate Cox proportional hazards models** to identify significant predictors.  
   * **Multivariate Cox model** including significant covariates; check proportional-hazards assumption (Schoenfeld residuals).  
3. **Model Interpretation**  
   * Hazard ratios with confidence intervals.  
   * Partial dependence or variable-importance plots for RSF.  
   * Discussion of clinical implications and limitations.

**Timeline:**

| Week 1 (Nov 16–22) | EDA and Begin Modelling | Nov 22 |
| :---: | :---: | :---: |
| **Week 2 (Nov 23–29)** | **Complete Modelling and Begin Interpretation** | **Nov 29** |
| **Week 3 (Nov 30–Dec 5\)** | **Complete Interpretation and Prepare Slides** | **Dec 5 – Presentation Due** |
| **Week 4 (Dec 6–12)** | **Finalize Report** | **Dec 12 – Report Due** |

**Role of Each Member:**  
Cameron Chesbrough: Drafting and editing the final report    
Jing Lyu: Descriptive analysis on the data and design the analysis plan  
Meitong Zhou: Filtering the cited paper and finishing the method part  
Ou Zhengkun: building models and interpreting the model.   
Xuange Liang: Modeling and Interpreting the result

Coding: Zhengkun & Jing  
Manuscript: Meitong & Cameron & Alex  
