library(dplyr)
library(gtsummary)
library(gt)


# ── 1. Load data ─────────────────────────────────────────────

epic     <- readRDS("data/epic.rds")     %>% mutate(system = "EPIC")
pre_epic <- readRDS("data/pre_epic.rds") %>% mutate(system = "Pre-EPIC")
ehr_raw  <- bind_rows(epic, pre_epic)

cat("Raw dimensions:\n")
print(dim(ehr_raw))
print(table(ehr_raw$system))


# ── 2. Variable name checks + export log ─────────────────────

dir.create("output", showWarnings = FALSE)

sink("output/00_variable_checks.txt")
cat("--- stillbirth / miscarriage ---\n")
print(grep("miscarr|stillbirth|abort", names(ehr_raw), ignore.case = TRUE, value = TRUE))
cat("\n--- diabetes ---\n")
print(grep("diab|gdm|gestdiab",        names(ehr_raw), ignore.case = TRUE, value = TRUE))
cat("\n--- HTN ---\n")
print(grep("htn|hypert|chronic|preexist|preeclamp",
           names(ehr_raw), ignore.case = TRUE, value = TRUE))
cat("\n--- geography (county / FIPS / zip / geocode / lat / lon) ---\n")
print(grep("county|fips|zip|geocode|address|lat|lon",
           names(ehr_raw), ignore.case = TRUE, value = TRUE))
sink()

cat("--- stillbirth / miscarriage ---\n")
print(grep("miscarr|stillbirth|abort", names(ehr_raw), ignore.case = TRUE, value = TRUE))
cat("--- diabetes ---\n")
print(grep("diab|gdm|gestdiab",        names(ehr_raw), ignore.case = TRUE, value = TRUE))
cat("--- HTN ---\n")
print(grep("htn|hypert|chronic|preexist|preeclamp",
           names(ehr_raw), ignore.case = TRUE, value = TRUE))
cat("--- geography ---\n")
print(grep("county|fips|zip|geocode|address|lat|lon",
           names(ehr_raw), ignore.case = TRUE, value = TRUE))


# ── 3. Cohort exclusions ─────────────────────────────────────

## Step 1: Exclude miscarriage / stillbirth
ehr_step1 <- ehr_raw %>%
  filter(
    if ("miscarriage"  %in% names(.)) miscarriage  != 1 else TRUE,
    if ("stillbirth"   %in% names(.)) stillbirth   != 1 else TRUE,
    if ("stillbirth_s" %in% names(.)) stillbirth_s != 1 else TRUE
  )

## Step 2: Restrict to delivery years 2016-2023

ehr_step2 <- ehr_step1 %>%
  mutate(delivery_year = as.integer(format(date_delivery, "%Y"))) %>%
  filter(delivery_year >= 2016 & delivery_year <= 2023)

cat("\nDelivery year distribution after 2016-2023 restriction:\n")
print(table(ehr_step2$delivery_year))

## Step 3: 5-county Atlanta metropolitan core

atlanta_counties <- c("fulton", "dekalb", "gwinnett", "cobb", "clayton")

cat("\npatient_address_county_desc value distribution (top 20, pre-filter):\n")
print(sort(table(ehr_step2$patient_address_county_desc, useNA = "ifany"),
           decreasing = TRUE)[1:20])

ehr_step3 <- ehr_step2 %>%
  filter(
    tolower(trimws(patient_address_county_desc)) %in% atlanta_counties |
    is.na(patient_address_county_desc)
  )

cat(sprintf("\nn after 5-county restriction: %d (excluded: %d)\n",
            nrow(ehr_step3), nrow(ehr_step2) - nrow(ehr_step3)))
cat("County distribution after restriction:\n")
print(table(ehr_step3$patient_address_county_desc, useNA = "ifany"))

## Step 4: Maternal age - summary only, no restriction applied

cat("\nMaternal age summary (no restriction applied):\n")
print(summary(ehr_step3$age))
cat("Age distribution (5-year bins):\n")
print(table(cut(ehr_step3$age, breaks = c(0,14,19,24,29,34,39,44,49,Inf),
                right = TRUE, include.lowest = TRUE), useNA = "ifany"))

## Step 5: Define main analytic cohort

ehr_clean <- ehr_step3

## Step 6: Sensitivity cohort - exclude pre-existing DM

cat("\ndiab distribution before sensitivity cohort split:\n")
print(table(ehr_step3$diab, useNA = "ifany"))

ehr_clean_no_preDM <- ehr_step3 %>% filter(diab == 0)
saveRDS(ehr_clean_no_preDM, "data/ehr_clean_no_preDM.rds")

## Cohort flow table
flow <- tibble(
  stage = c(
    "01_raw",
    "02_exclude_miscarriage_stillbirth",
    "03_restrict_delivery_year_2016_2023",
    "04_restrict_5county_atlanta",
    "05_main_cohort_no_age_restriction",
    "06_sensitivity_exclude_preexisting_DM"
  ),
  n = c(
    nrow(ehr_raw),
    nrow(ehr_step1),
    nrow(ehr_step2),
    nrow(ehr_step3),
    nrow(ehr_clean),
    nrow(ehr_clean_no_preDM)
  ),
  n_excluded = c(
    NA,
    nrow(ehr_raw)    - nrow(ehr_step1),
    nrow(ehr_step1)  - nrow(ehr_step2),
    nrow(ehr_step2)  - nrow(ehr_step3),
    nrow(ehr_step3)  - nrow(ehr_clean),
    nrow(ehr_clean)  - nrow(ehr_clean_no_preDM)
  ),
  note = c(
    "Combined EPIC + Pre-EPIC",
    "Excluded miscarriage and stillbirth records",
    "Restricted to 2016-2023 delivery years",
    "5-county filter ACTIVE: Pre-EPIC restricted to Fulton/DeKalb/Gwinnett/Cobb/Clayton (case-insensitive); EPIC records retained regardless (county field unpopulated in EPIC system)",
    "No maternal age restriction applied; see age summary above",
    "Sensitivity cohort: restricted to confirmed diab==0 records; diab==1 and diab==NA excluded"
  )
)
print(flow)
write.csv(flow, "output/01_cohort_flow.csv", row.names = FALSE)
saveRDS(ehr_clean, "data/ehr_clean.rds")


# ── 5. Create analysis variables ─────────────────────────────

ehr_clean <- ehr_clean %>%
  mutate(
    delivery_year = as.integer(format(date_delivery, "%Y")),
    system        = factor(system, levels = c("Pre-EPIC", "EPIC")),
    race          = factor(race),
    obesity       = factor(obesity,  levels = c(0, 1), labels = c("No", "Yes")),
    diab          = factor(diab,     levels = c(0, 1), labels = c("No", "Yes")),
    gestdiab      = factor(gestdiab, levels = c(0, 1), labels = c("No", "Yes")),
    htn           = factor(htn,      levels = c(0, 1), labels = c("No", "Yes")),
    preecl        = factor(
      as.integer(preecl_s == 1 | preecl_n == 1),
      levels = c(0, 1), labels = c("No", "Yes")
    )
  )

cat("\ngestdiab distribution (including NA):\n")
print(table(ehr_clean$gestdiab, useNA = "ifany"))

cat("\nrace distribution (top 10):\n")
print(sort(table(ehr_clean$race, useNA = "ifany"), decreasing = TRUE)[1:10])

ehr_clean <- ehr_clean %>%
  mutate(race = na_if(trimws(as.character(race)), ""),
         race = factor(race))

cat("\nrace distribution after empty-string cleanup (top 10):\n")
print(sort(table(ehr_clean$race, useNA = "ifany"), decreasing = TRUE)[1:10])


# ── 6. Missing data summary ───────────────────────────────────

missing_summary <- ehr_clean %>%
  summarise(
    n_total               = n(),
    age_missing           = sum(is.na(age)),
    age_missing_pct       = round(mean(is.na(age)) * 100, 2),
    race_missing          = sum(is.na(race)),
    race_missing_pct      = round(mean(is.na(race)) * 100, 2),
    obesity_missing       = sum(is.na(obesity)),
    obesity_missing_pct   = round(mean(is.na(obesity)) * 100, 2),
    diab_missing          = sum(is.na(diab)),
    diab_missing_pct      = round(mean(is.na(diab)) * 100, 2),
    gest_age_missing      = sum(is.na(gest_age)),
    gest_age_missing_pct  = round(mean(is.na(gest_age)) * 100, 2),
    gestdiab_missing      = sum(is.na(gestdiab)),
    gestdiab_missing_pct  = round(mean(is.na(gestdiab)) * 100, 2),
    htn_missing           = sum(is.na(htn)),
    htn_missing_pct       = round(mean(is.na(htn)) * 100, 2),
    preecl_s_missing      = sum(is.na(preecl_s)),
    preecl_s_missing_pct  = round(mean(is.na(preecl_s)) * 100, 2),
    preecl_n_missing      = sum(is.na(preecl_n)),
    preecl_n_missing_pct  = round(mean(is.na(preecl_n)) * 100, 2),
    date_delivery_missing     = sum(is.na(date_delivery)),
    date_delivery_missing_pct = round(mean(is.na(date_delivery)) * 100, 2)
  )

print(t(missing_summary))   # transpose for easier reading in console
write.csv(missing_summary, "output/00_missing_summary.csv", row.names = FALSE)
cat("\nMissing data summary saved to output/00_missing_summary.csv\n")


# ── 7. Outcome rates by system and year ──────────────────────

rates_by_system <- ehr_clean %>%
  group_by(system) %>%
  summarise(
    n           = n(),
    gdm_pct     = round(mean(gestdiab == "Yes", na.rm = TRUE) * 100, 1),
    htn_pct     = round(mean(htn      == "Yes", na.rm = TRUE) * 100, 1),
    preecl_pct  = round(mean(preecl   == "Yes", na.rm = TRUE) * 100, 1),
    obesity_pct = round(mean(obesity  == "Yes", na.rm = TRUE) * 100, 1),
    diab_pct    = round(mean(diab     == "Yes", na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
print(rates_by_system)
write.csv(rates_by_system, "output/02_rates_by_system.csv", row.names = FALSE)

rates_by_year_system <- ehr_clean %>%
  group_by(system, delivery_year) %>%
  summarise(
    n        = n(),
    gdm_pct  = round(mean(gestdiab == "Yes", na.rm = TRUE) * 100, 1),
    htn_pct  = round(mean(htn      == "Yes", na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(system, delivery_year)
print(rates_by_year_system)
write.csv(rates_by_year_system, "output/03_rates_by_year_system.csv", row.names = FALSE)


# ── 8. Table 1: Overall descriptives ─────────────────────────

t1 <- ehr_clean %>%
  select(system, age, race, obesity, diab, gest_age, gestdiab, htn, preecl) %>%
  tbl_summary(
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits  = all_continuous() ~ 1,
    missing = "ifany"
  ) %>%
  bold_labels()

t1
t1 %>% as_gt() %>% gtsave("output/04_table1_overall.html")


# ── 9. Table 2: Stratified by GDM ────────────────────────────

t2 <- ehr_clean %>%
  select(system, age, race, obesity, diab, gest_age, htn, preecl, gestdiab) %>%
  tbl_summary(
    by = gestdiab,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits  = all_continuous() ~ 1,
    missing = "ifany"
  ) %>%
  add_overall() %>%
  bold_labels()

t2
t2 %>% as_gt() %>% gtsave("output/05_table2_by_gdm.html")


# ── 10. Table 3: Stratified by HTN ───────────────────────────

t3 <- ehr_clean %>%
  select(system, age, race, obesity, diab, gest_age, gestdiab, preecl, htn) %>%
  tbl_summary(
    by = htn,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits  = all_continuous() ~ 1,
    missing = "ifany"
  ) %>%
  add_overall() %>%
  bold_labels()

t3
t3 %>% as_gt() %>% gtsave("output/06_table3_by_htn.html")


# ── 11. P-values ──────────────────────────────────────────────

## Continuous variables: Welch two-sample t-test
cont_vars <- c("age", "gest_age")

pvals_cont <- lapply(cont_vars, function(v) {
  x  <- ehr_clean[[v]]
  g  <- ehr_clean$gestdiab
  x1 <- x[g == "Yes" & !is.na(x)]
  x0 <- x[g == "No"  & !is.na(x)]
  if (length(x1) < 2 || length(x0) < 2) {
    return(data.frame(variable = v, test = "t.test", p = NA))
  }
  data.frame(variable = v, test = "t.test", p = signif(t.test(x1, x0)$p.value, 4))
}) %>% bind_rows()

cat_vars <- c("system", "race", "obesity", "diab", "htn", "preecl")

pvals_cat <- lapply(cat_vars, function(v) {
  tbl <- table(ehr_clean[[v]], ehr_clean$gestdiab, useNA = "no")
  chi <- suppressWarnings(chisq.test(tbl))
  use_fisher <- any(chi$expected < 5)

  if (use_fisher) {
    if (all(dim(tbl) == c(2, 2))) {
      p    <- fisher.test(tbl)$p.value
      test <- "fisher.test_exact"
    } else {
      set.seed(2024)
      p    <- fisher.test(tbl, simulate.p.value = TRUE, B = 10000)$p.value
      test <- "fisher.test_sim"
    }
  } else {
    p    <- chi$p.value
    test <- "chisq.test"
  }
  data.frame(variable = v, test = test, p = signif(p, 4))
}) %>% bind_rows()

pvals_all <- bind_rows(pvals_cont, pvals_cat)
print(pvals_all)
write.csv(pvals_all, "output/07_pvalues_by_gdm.csv", row.names = FALSE)

cat("\nDescriptive outputs saved to output/ folder.\n")


# ── 12. Logistic regression helper function ───────────────────

or_table <- function(model, label = "") {
  s     <- summary(model)$coef
  terms <- rownames(s)

  int_idx       <- which(terms == "(Intercept)")
  int_logodds   <- s[int_idx, "Estimate"]
  int_se        <- s[int_idx, "Std. Error"]
  baseline_prob <- plogis(int_logodds)

  baseline_lo   <- plogis(int_logodds - 1.96 * int_se)
  baseline_hi   <- plogis(int_logodds + 1.96 * int_se)

  cat(sprintf(
    "\n[%s] Intercept:\n  log-odds = %.4f\n  Baseline probability = %.4f (%.2f%%)\n  95%% CI: (%.2f%%, %.2f%%)\n  (ref: Black/AA, mean age, year 2016, no obesity, no pre-DM, Pre-EPIC)\n",
    label, int_logodds, baseline_prob, baseline_prob * 100,
    baseline_lo * 100, baseline_hi * 100
  ))

  intercept_row <- data.frame(
    term             = "(Intercept)",
    estimate_logodds = round(int_logodds,   4),
    baseline_prob    = round(baseline_prob, 4),
    baseline_prob_lo = round(baseline_lo,   4),
    baseline_prob_hi = round(baseline_hi,   4),
    OR               = NA_real_,
    OR_lo            = NA_real_,
    OR_hi            = NA_real_,
    p                = signif(s[int_idx, "Pr(>|z|)"], 4)
  )

  coef_rows <- data.frame(
    term             = terms[-int_idx],
    estimate_logodds = NA_real_,
    baseline_prob    = NA_real_,
    baseline_prob_lo = NA_real_,
    baseline_prob_hi = NA_real_,
    OR    = round(exp(s[-int_idx, "Estimate"]),                                    3),
    OR_lo = round(exp(s[-int_idx, "Estimate"] - 1.96 * s[-int_idx, "Std. Error"]), 3),
    OR_hi = round(exp(s[-int_idx, "Estimate"] + 1.96 * s[-int_idx, "Std. Error"]), 3),
    p     = signif(s[-int_idx, "Pr(>|z|)"], 4),
    row.names = NULL
  )

  bind_rows(intercept_row, coef_rows)
}


# ── 13. Main cohort: baseline logistic regression ─────────────

age_mean <- mean(ehr_clean$age, na.rm = TRUE)
cat(sprintf("\nCentering: mean age = %.2f years; delivery_year centered at 2016\n", age_mean))

ehr_m <- ehr_clean %>%
  mutate(
    race       = relevel(factor(race), ref = "Black or African American"),
    age_c      = age - age_mean,
    year_c     = delivery_year - 2016,
    gdm        = as.integer(gestdiab == "Yes"),
    htn_y      = as.integer(htn      == "Yes"),
    preecl_any = as.integer(preecl_s == 1 | preecl_n == 1),
    hdp_any    = as.integer((htn     == "Yes") | (preecl_s == 1) | (preecl_n == 1))
  )

cat("\npreecl_any distribution:\n"); print(table(ehr_m$preecl_any, useNA = "ifany"))
cat("\nhdp_any distribution:\n");    print(table(ehr_m$hdp_any,    useNA = "ifany"))

cat("\nMissingness in preecl_s / preecl_n (NA propagates to preecl_any and hdp_any):\n")
print(colSums(is.na(ehr_m[, c("preecl_s", "preecl_n")])))

m_gdm    <- glm(gdm        ~ age_c + race + obesity + diab + system + year_c,
                data = ehr_m, family = binomial())
m_htn    <- glm(htn_y      ~ age_c + race + obesity + diab + system + year_c,
                data = ehr_m, family = binomial())
m_preecl <- glm(preecl_any ~ age_c + race + obesity + diab + system + year_c,
                data = ehr_m, family = binomial())
m_hdp    <- glm(hdp_any    ~ age_c + race + obesity + diab + system + year_c,
                data = ehr_m, family = binomial())

write.csv(or_table(m_gdm,    "GDM"),    "output/08_baseline_OR_gdm.csv",    row.names = FALSE)
write.csv(or_table(m_htn,    "HTN"),    "output/09_baseline_OR_htn.csv",    row.names = FALSE)
write.csv(or_table(m_preecl, "Preecl"), "output/12_baseline_OR_preecl.csv", row.names = FALSE)
write.csv(or_table(m_hdp,    "HDP"),    "output/13_baseline_OR_hdp.csv",    row.names = FALSE)

cat("\nMain cohort models complete.\n")


# ── 14. Sensitivity cohort: exclude pre-existing DM ──────────

ehr_sens <- readRDS("data/ehr_clean_no_preDM.rds") %>%
  mutate(
    delivery_year = as.integer(format(date_delivery, "%Y")),
    system        = factor(system,   levels = c("Pre-EPIC", "EPIC")),
    race          = relevel(factor(race), ref = "Black or African American"),
    obesity       = factor(obesity,  levels = c(0, 1), labels = c("No", "Yes")),
    gestdiab      = factor(gestdiab, levels = c(0, 1), labels = c("No", "Yes")),
    htn           = factor(htn,      levels = c(0, 1), labels = c("No", "Yes")),
    gdm           = as.integer(gestdiab == "Yes"),
    htn_y         = as.integer(htn      == "Yes"),

    age_c         = age - age_mean,
    year_c        = delivery_year - 2016
  )

m_gdm_sens <- glm(gdm   ~ age_c + race + obesity + system + year_c,
                  data = ehr_sens, family = binomial())
m_htn_sens <- glm(htn_y ~ age_c + race + obesity + system + year_c,
                  data = ehr_sens, family = binomial())

write.csv(or_table(m_gdm_sens, "GDM_sensitivity"),
          "output/10_sensitivity_OR_gdm_no_preDM.csv", row.names = FALSE)
write.csv(or_table(m_htn_sens, "HTN_sensitivity"),
          "output/11_sensitivity_OR_htn_no_preDM.csv", row.names = FALSE)

cat("\nSensitivity cohort models complete.\n")


# ── 15. Model sample sizes ────

model_n <- tibble(
  model   = c("m_gdm", "m_htn", "m_preecl", "m_hdp",
               "m_gdm_sens", "m_htn_sens"),
  outcome = c("GDM", "Gestational HTN", "Preeclampsia", "HDP composite",
               "GDM (sensitivity)", "HTN (sensitivity)"),
  cohort  = c(rep("Main", 4), rep("Sensitivity (no pre-DM)", 2)),
  n_used  = c(
    nobs(m_gdm),
    nobs(m_htn),
    nobs(m_preecl),
    nobs(m_hdp),
    nobs(m_gdm_sens),
    nobs(m_htn_sens)
  )
)

print(model_n)
write.csv(model_n, "output/15_model_n.csv", row.names = FALSE)

cat("\nAll outputs saved to output/ folder.\n")
