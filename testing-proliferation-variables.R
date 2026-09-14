library(survival)
library(broom)

# list of variables to test corresponding to columns in df called patient_dat
uva_vars <- c("AR_HM", "ER_E_HM", "SPOP_mutated_Liu", "TP53_loss_Chipidza", "CCP", "met_stage")

# UVA of each variable against OS
uva_results <- map_df(uva_vars, function(v) {
  f <- as.formula(paste("Surv(days_from_biopsy_to_last_contact, death) ~", v))
  fit <- coxph(f, data = patient_dat)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(variable = v, HR = estimate, lower_CI = conf.low, 
              upper_CI = conf.high, p_value = p.value)
})

uva_results <- uva_results %>% arrange(p_value)
print(uva_results)

# MVA including significant variables from the UVA (met_stage, CCP, TP53_loss_Chipidza)
mva_formula <- as.formula("Surv(days_from_biopsy_to_last_contact, death) ~ met_stage + CCP + TP53_loss_Chipidza")
fit_mva <- coxph(mva_formula, data = patient_dat)
summary(fit_mva)

# Check proportional hazards assumption
cox.zph(fit_mva)
