library(dplyr)
library(broom)

data_dir   <- "R:/hhchang/Epidemiology/NHAPPS/Analysis/Feiran/Emory EHR/20260217/data"
output_dir <- "R:/hhchang/Epidemiology/NHAPPS/Analysis/Feiran/Emory EHR/20260309/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load data

dat <- readRDS(file.path(data_dir, "epic_final.rds"))
cat(sprintf("Loaded: %d records, %d columns\n", nrow(dat), ncol(dat)))

age_mean  <- mean(dat$age, na.rm = TRUE)
year_mean <- 2016

dat <- dat %>%
  mutate(
    # Outcomes (binary)
    gdm_bin    = as.integer(gestdiab == 1),
    gh_bin     = as.integer(htn == 1),
    preecl_bin = as.integer(preecl_s == 1 | preecl_n == 1),
    hdp_bin    = as.integer(htn == 1 | preecl_s == 1 | preecl_n == 1),

    # Race: merge small groups into Other/Unknown to avoid separation
    race_grp = case_when(
      race == "Black or African American" ~ "Black",
      race == "White"                     ~ "White",
      race == "Asian"                     ~ "Asian",
      TRUE                                ~ "Other/Unknown"
    ),
    race_grp = relevel(factor(race_grp), ref = "Black"),

    # Other covariates
    obesity    = factor(obesity, levels = c(0, 1)),
    diab       = factor(diab,    levels = c(0, 1)),
    age_c      = age - age_mean,
    year_c     = delivery_year - year_mean,
    income_10k = median_household_income / 10000
  )

cat(sprintf("age centered at: %.2f years\n", age_mean))
cat(sprintf("delivery_year centered at: %d\n", year_mean))
cat("race_grp distribution:\n")
print(table(dat$race_grp))


# 2. Define exposures, outcomes, covariates 

exposure_vars <- c(
  "t1_T2_corr_mean",      "t2_T2_corr_mean",      "t3_T2_corr_mean",
  "t1_HI_corr_mean",      "t2_HI_corr_mean",      "t3_HI_corr_mean",
  "t1_WBT_corr_mean",     "t2_WBT_corr_mean",     "t3_WBT_corr_mean",
  "t1_HUMIDEX_corr_mean", "t2_HUMIDEX_corr_mean", "t3_HUMIDEX_corr_mean",
  "t1_T2_min_mean",       "t2_T2_min_mean",       "t3_T2_min_mean"
)

outcomes <- c("gdm_bin", "gh_bin", "preecl_bin", "hdp_bin")

covars <- c("age_c", "race_grp", "obesity", "diab", "year_c", "income_10k")


# 3. Pre-modeling QC 

cat("\n── Pre-modeling missingness check ──\n")
check_vars <- c(outcomes, covars, "t1_T2_corr_mean")
miss_check <- sapply(dat[check_vars], function(x) sum(is.na(x)))
print(miss_check)

cat(sprintf("\nComplete cases (all covariates + income): %d\n",
            sum(complete.cases(dat[c(covars, outcomes)]))))


# 4. Helper function

extract_or <- function(fit, exposure, outcome, n) {
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::filter(term == exposure) %>%
    dplyr::mutate(
      outcome  = outcome,
      exposure = exposure,
      n        = n,
      OR       = round(estimate,  3),
      CI_lo    = round(conf.low,  3),
      CI_hi    = round(conf.high, 3),
      p        = signif(p.value,  3)
    ) %>%
    dplyr::select(outcome, exposure, n, OR, CI_lo, CI_hi, p)
}


# 5. Run models 

cat("\n── Running logistic regression models ──\n")

results_crude    <- list()
results_adjusted <- list()
k <- 1

for (y in outcomes) {
  for (x in exposure_vars) {

    dsub_crude <- dat %>%
      dplyr::select(all_of(c(y, x))) %>%
      na.omit()

    dsub_adj <- dat %>%
      dplyr::select(all_of(c(y, x, covars))) %>%
      na.omit()

    if (nrow(dsub_crude) == 0 | nrow(dsub_adj) == 0) next

    # Crude model
    fit_crude <- glm(as.formula(paste(y, "~", x)),
                     data = dsub_crude, family = binomial())
    results_crude[[k]] <- extract_or(fit_crude, x, y, nrow(dsub_crude))

    # Adjusted model
    fml_adj <- as.formula(
      paste(y, "~", paste(c(x, covars[covars %in% names(dsub_adj)]),
                          collapse = " + "))
    )
    fit_adj <- glm(fml_adj, data = dsub_adj, family = binomial())
    results_adjusted[[k]] <- extract_or(fit_adj, x, y, nrow(dsub_adj))

    k <- k + 1
  }
}

results_crude_df    <- dplyr::bind_rows(results_crude)
results_adjusted_df <- dplyr::bind_rows(results_adjusted)

cat(sprintf("Total models run: %d crude, %d adjusted\n",
            nrow(results_crude_df), nrow(results_adjusted_df)))


# 6. Save results 

write.csv(results_crude_df,
          file.path(output_dir, "22_models_crude.csv"),
          row.names = FALSE)
write.csv(results_adjusted_df,
          file.path(output_dir, "23_models_adjusted.csv"),
          row.names = FALSE)

cat("\nResults saved.\n")
cat("\n── Adjusted model results: GDM ──\n")
print(results_adjusted_df %>% dplyr::filter(outcome == "gdm_bin"))
cat("\n── Adjusted model results: GH ──\n")
print(results_adjusted_df %>% dplyr::filter(outcome == "gh_bin"))
