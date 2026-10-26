
# =============================================================================
# CORAL STUDY: DELAYED TREATMENT-SEEKING FOR VISCERAL LEISHMANIASIS
# Prospective quantitative data analysis
#
# Primary outcome:
# Delayed care-seeking = >=30 days from symptom onset to first care
#
# Workflow:
# 1. Import data
# 2. Assess data quality
# 3. Clean data and define outcome
# 4. Recode independent variables
# 5. Rename REDCap checkbox variables
# 6. Descriptive analysis
# 7. Crude logistic regression
# 8. Adjusted logistic regression
# 9. Export results
# =============================================================================


# 01. SETUP -------------------------------------------------------------------

rm(list = ls())

setwd("~/CORAL KIIs/CORAL Analysis/CORAL_Delayed_Rx_Seeking_Meron/CORAL_Quantitative/CORAL Quantitative")

# Load packages
packages <- c(
  "tidyverse",
  "gtsummary",
  "broom",
  "broom.helpers"
)

missing_packages <- setdiff(packages, rownames(installed.packages()))

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

# set setseed
set.seed(2026)
# set output dir
output_dir <- "CORAL_results"
dir.create(output_dir, showWarnings = FALSE)

save_csv <- function(data, filename) {
  readr::write_csv(
    data,
    file.path(output_dir, filename),
    na = ""
  )
}


# 02. IMPORT DATA -------------------------------------------------------------
#load the data downloaded from the RedCap database 
raw <- read_csv(
  "CORAL_Delayed_Rx.csv",
  show_col_types = FALSE
)

colnames(raw)

#remove irrelevant variables from the data. this include time, date, etc variables
raw <- raw[ , c(-2, -3, -4, -6, -7, -8, -9, -10, -14, -16, -31, -37, -39,
             -48, -51, -54, -62, -63, -65, -75, -84, -85, -92, -96, -97, -107, -108, -110, 
             -116, -120, -121)]

# Inspect dataset
glimpse(raw)
dim(raw)
names(raw)

# Confirm required variables
required <- c(
  "record_id",
  "age",
  "gender",
  "illness_duration_pre_care",
  "referred_from_another",
  "how_long_to_referral"
)

stopifnot(all(required %in% names(raw)))

# Save original data
saveRDS(raw, file.path(output_dir, "01_original_data.rds"))


# 03. DATA QUALITY ASSESSMENT -------------------------------------------------

# A. Uniqueness: duplicate IDs
duplicates <- raw %>%
  filter(!is.na(record_id)) %>%
  group_by(record_id) %>%
  filter(n() > 1) %>%
  ungroup()

save_csv(duplicates, "02_duplicate_records.csv")


# B. Variable-level completeness
completeness <- tibble(
  variable = names(raw),
  n_total = nrow(raw),
  n_missing = sapply(raw, function(x) {
    sum(is.na(x) | (is.character(x) & trimws(x) == ""))
  })
) %>%
  mutate(
    n_observed = n_total - n_missing,
    pct_complete = round(100 * n_observed / n_total, 1)
  ) %>%
  arrange(pct_complete)

print(completeness)

save_csv(completeness, "03_variable_completeness.csv")


# C. Record-level completeness
record_completeness <- raw %>%
  mutate(
    n_missing = rowSums(is.na(raw)),
    pct_complete = round(
      100 * (ncol(raw) - n_missing) / ncol(raw), 1
    )
  ) %>%
  select(record_id, n_missing, pct_complete)

save_csv(record_completeness, "04_record_completeness.csv")


# D. Overall completeness
overall_completeness <- round(
  100 * sum(!is.na(raw)) / (nrow(raw) * ncol(raw)),
  1
) #94.5% completness

cat("Overall completeness:", overall_completeness, "%\n")

# NOTE:
# REDCap branching logic can produce legitimate missing values.
# Review structural missingness before interpreting completeness.


# E. Validity and skip-logic checks
qc_flags <- raw %>%
  transmute(
    record_id,
    
    missing_outcome = is.na(illness_duration_pre_care),
    
    invalid_duration =
      !is.na(illness_duration_pre_care) &
      illness_duration_pre_care <= 0,
    
    very_long_duration =
      !is.na(illness_duration_pre_care) &
      illness_duration_pre_care > 365,
    
    invalid_age =
      !is.na(age) & (age < 0 | age > 120),
    
    negative_distance =
      !is.na(distance_health_facility) &
      distance_health_facility < 0,
    
    negative_cost =
      !is.na(treatment_cost_before_here) &
      treatment_cost_before_here < 0,
    
    referred_missing_time =
      !is.na(referred_from_another) &
      referred_from_another == 1 &
      is.na(how_long_to_referral)
  )

save_csv(qc_flags, "05_quality_flags.csv")


# 04. CLEAN DATA AND DEFINE OUTCOMES -----------------------------------------

# Preserve all participants.
# Invalid values are set to NA only in derived variables.

dat <- raw %>%
  mutate(
    
    # Primary duration
    duration_days = if_else(
      !is.na(illness_duration_pre_care) &
        illness_duration_pre_care > 0,
      as.numeric(illness_duration_pre_care),
      NA_real_
    ),
    
    # Primary outcome
    delayed_care = case_when(
      is.na(duration_days) ~ NA_integer_,
      duration_days >= 30 ~ 1L,
      TRUE ~ 0L
    ),
    
    delayed_care_factor = factor(
      delayed_care,
      levels = c(0, 1),
      labels = c("Not delayed", "Delayed")
    ),
    
    # Valid age
    age_years = if_else(
      !is.na(age) & age >= 0 & age <= 120,
      as.numeric(age),
      NA_real_
    ),
    
    # Secondary outcome
    referral_duration_days = if_else(
      !is.na(referred_from_another) &
        referred_from_another == 1 &
        !is.na(how_long_to_referral) &
        how_long_to_referral >= 0,
      as.numeric(how_long_to_referral),
      NA_real_
    ),
    
    delayed_referral = case_when(
      is.na(referred_from_another) ~ NA_integer_,
      referred_from_another != 1 ~ NA_integer_,
      is.na(referral_duration_days) ~ NA_integer_,
      referral_duration_days > 7 ~ 1L,
      TRUE ~ 0L
    )
  )

# IMPORTANT:
# Confirm that illness_duration_pre_care and
# how_long_to_referral are recorded in DAYS.
# Also verify whether zero days is valid before excluding it.

summary(raw$distance_health_facility)

# 05. Data preprocessing for analysis & RECODE INDEPENDENT VARIABLES WHERE REQUIRED AND RELEVANT-------------------------------------------

dat <- dat %>%
  mutate(
    
    # A. Age categories (<18, 18-29, 30+) with NAs replaced by the mode
    age_cat = case_when(
      age_years < 18 ~ "<18",
      age_years >= 18 & age_years <= 29 ~ "18-29",
      age_years >= 30 & age_years <= 120 ~ "30+",
      TRUE ~ NA_character_
    ),
    
    # Replace NA age categories with the modal category
    age_cat = {
      valid_age <- age_cat[!is.na(age_cat)]
      if (length(valid_age) > 0) {
        uniqv <- unique(valid_age)
        mod_val <- uniqv[which.max(tabulate(match(valid_age, uniqv)))]
        if_else(is.na(age_cat), mod_val, age_cat)
      } else {
        age_cat
      }
    },
    
    age_cat = factor(
      age_cat,
      levels = c("<18", "18-29", "30+")
    ), 
    
    # B. Distance tertiles
    # B. Distance categorized by 3 hours (180 minutes)
    distance_valid = if_else(
      distance_health_facility >= 0,
      as.numeric(distance_health_facility),
      NA_real_
    ),
    
    distance_cat = case_when(
      distance_valid < 120 ~ "< 2 hours (< 120 mins)",
      distance_valid >= 120 ~ ">= 2 hours (>= 120 mins)",
      TRUE ~ NA_character_
    ),
    
    distance_cat = factor(
      distance_cat,
      levels = c("< 2 hours (< 120 mins)", ">= 2 hours (>= 120 mins)")
    ),
    
    # C. Family size categories
    family_size_cat = case_when(
      family_size >= 1 & family_size <= 3 ~ "Small (1-3)",
      family_size >= 4 & family_size <= 6 ~ "Medium (4-6)",
      family_size >= 7 ~ "Large (7+)",
      TRUE ~ NA_character_
    ),
    
    family_size_cat = factor(
      family_size_cat,
      levels = c(
        "Small (1-3)",
        "Medium (4-6)",
        "Large (7+)"
      )
    ),
    
    # D. Household income is coded categorically
    household_income = as.factor(household_income)
  )

# Check recoding
table(dat$age_cat, useNA = "ifany")
table(dat$distance_cat, useNA = "ifany")
table(dat$family_size_cat, useNA = "ifany")


# 06. Check if all vars have the right class  ---------------------------------------

as.factor(dat$health_center) #because it was number


# 07. RECODE, MERGE "OTHER"/"DON'T KNOW" TO MODE, AND FACTOR VARIABLES ------------

# Define a helper function to replace "Other" and "Don't know"/"Not sure" 
# categories with the statistical mode of the remaining valid responses.
# 06_B. RECODE, MERGE "OTHER"/"DON'T KNOW" TO MODE, AND FACTOR VARIABLES -----

# Define a helper function to replace "Other" and "Don't know"/"Not sure" 
# categories with the statistical mode of the remaining valid responses.
merge_to_mode <- function(x, codes_to_replace) {
  valid_x <- x[!is.na(x) & !(x %in% codes_to_replace)]
  if (length(valid_x) == 0) return(x)
  uniqv <- unique(valid_x)
  mod_val <- uniqv[which.max(tabulate(match(valid_x, uniqv)))]
  x[x %in% codes_to_replace] <- mod_val
  return(x)
}

dat <- dat %>%
  mutate(
    # 1. Marital status (98 = Other)[cite: 1]
    marital_status = merge_to_mode(marital_status, c(98)),
    marital_status = factor(
      marital_status,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Single", "Married", "Divorced", "Widowed", "Not applicable (child)")
    ),
    
    # 2. Occupational status (98 = Other)[cite: 1]
    occupational_status = merge_to_mode(occupational_status, c(98)),
    occupational_status = factor(
      occupational_status,
      levels = c(1, 2, 3, 4, 5, 6, 7),
      labels = c("Farmer/pastoralist", "Trader", "Daily laborer", "Government employee", "Unemployed", "Student", "Not applicable (child)")
    ),
    
    # 3. Household income (5 = Don't know)[cite: 1]
    household_income = merge_to_mode(household_income, c(5)),
    household_income = factor(
      household_income,
      levels = c(1, 2, 3, 4),
      labels = c("< 10000 ETB", "10000-30000 ETB", "30000-50000 ETB", ">50000 ETB")
    ),
    
    # 4. Reason first care choice (98 = Other)[cite: 1]
    reason_first_care_choice = merge_to_mode(reason_first_care_choice, c(98)),
    reason_first_care_choice = factor(
      reason_first_care_choice,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Proximity", "Cost", "Trust in provider", "Referred", "Advice from others")
    ),
    
    # 5. Who advised delay reason (98 = Other)[cite: 1]
    who_advised_delay_reason = merge_to_mode(who_advised_delay_reason, c(98)),
    who_advised_delay_reason = factor(
      who_advised_delay_reason,
      levels = c(1, 2, 3),
      labels = c("Family member", "Community leader", "Traditional healer")
    ),
    
    # 6. Source VL info (98 = Other)[cite: 1]
    source_vl_info = merge_to_mode(source_vl_info, c(98)),
    source_vl_info = factor(
      source_vl_info,
      levels = c(1, 2, 3, 4, 5, 6),
      labels = c("Health worker", "Family/friends", "Radio/TV", "Community meeting", "School", "I previously had VL")
    ),
    
    # 7. VL preventable (3 = Don't know)[cite: 1]
    vl_preventable = merge_to_mode(vl_preventable, c(3)),
    vl_preventable = factor(
      vl_preventable,
      levels = c(1, 2),
      labels = c("Yes", "No")
    ),
    
    # 8. VL serious (3 = Don't know)[cite: 1]
    vl_serious = merge_to_mode(vl_serious, c(3)),
    vl_serious = factor(
      vl_serious,
      levels = c(1, 2),
      labels = c("Yes", "No")
    ),
    
    # 9. VL curable (3 = Don't know)[cite: 1]
    vl_curable = merge_to_mode(vl_curable, c(3)),
    vl_curable = factor(
      vl_curable,
      levels = c(1, 2),
      labels = c("Yes", "No")
    ),
    
    # 10. Transport mode (98 = Other)[cite: 1]
    transport_mode = merge_to_mode(transport_mode, c(98)),
    transport_mode = factor(
      transport_mode,
      levels = c(1, 2, 3, 4, 5),
      labels = c("On foot", "Public transport", "Motorcycle", "Ambulance", "Animal (donkey/horse)")
    ),
    
    # # 11. When plan go (5 = Not sure)[cite: 1]
    # when_plan_go = merge_to_mode(when_plan_go, c(5)),
    # when_plan_go = factor(
    #   when_plan_go,
    #   levels = c(1, 2, 3, 4),
    #   labels = c("Today", "Tomorrow", "In 2-4 days", "Later than 7 days (1 Week)")
    # ),
    # 
    # 12. Perceive need urgent care (3 = Not sure)[cite: 1]
    perceive_need_urgent_care = merge_to_mode(perceive_need_urgent_care, c(3)),
    perceive_need_urgent_care = factor(
      perceive_need_urgent_care,
      levels = c(1, 2),
      labels = c("Yes - it's serious", "No - I feel fine now")
    ),
    
    # 13. Referred before similar (4 = Don't remember)[cite: 1]
    referred_before_similar = merge_to_mode(referred_before_similar, c(4)),
    referred_before_similar = factor(
      referred_before_similar,
      levels = c(1, 2, 3),
      labels = c("Yes - and I went", "Yes - but I did not go", "No, first time")
    ),
    
    # Standard variables without "Other" or "Don't know" codes
    health_center = factor(health_center, levels = c(1, 2, 3, 4), labels = c("Alduba HC", "Turmi HC", "Kangatin HC", "Assile HC")),
    gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    residence = factor(residence, levels = c(1, 2, 3), labels = c("Rural", "Urban", "Semi-urban")),
    educational_status = factor(educational_status, levels = c(1, 2, 3, 4, 5), labels = c("No formal education", "Informal/religious education", "Primary (Grade 1-8)", "Secondary (Grade 9-12)", "Tertiary (college/university)")),
    #past_previous_sympts = factor(past_previous_sympts, levels = c(0, 1), labels = c("No", "Yes")),
    #referred_from_another = factor(referred_from_another, levels = c(0, 1), labels = c("No", "Yes")),
    advised_delay_formal = factor(advised_delay_formal, levels = c(0, 1), labels = c("No", "Yes")),
    heard_of_vl = factor(heard_of_vl, levels = c(0, 1), labels = c("No", "Yes")),
    transport_available = factor(transport_available, levels = c(0, 1), labels = c("No", "Yes")),
    transport_affordable = factor(transport_affordable, levels = c(0, 1), labels = c("No", "Yes")),
    know_where_to_go = factor(know_where_to_go, levels = c(0, 1), labels = c("No", "Yes")),
    facility_delay = factor(facility_delay, levels = c(0, 1), labels = c("No", "Yes")),
    borrow_or_sell_for_treatment = factor(borrow_or_sell_for_treatment, levels = c(0, 1), labels = c("No", "Yes")),
    illness_affect_work = factor(illness_affect_work, levels = c(1, 2, 3), labels = c("Yes", "No", "Not applicable (child)")),
    family_stop_work_for_care = factor(family_stop_work_for_care, levels = c(0, 1), labels = c("No", "Yes")),
    #plan_go_referral_facility = factor(plan_go_referral_facility, levels = c(1, 2, 3), labels = c("Yes", "Not sure", "No")),
    told_travel_reimbursement = factor(told_travel_reimbursement, levels = c(0, 1), labels = c("No", "Yes"))
  )

checkbox_names <- c(
  belief_sandfly_bite = "belief_cause_vl___1",
  belief_poor_hygiene = "belief_cause_vl___2",
  belief_drinking_dirty_water = "belief_cause_vl___3",
  belief_witchcraft = "belief_cause_vl___4",
  belief_other = "belief_cause_vl___98",
  
  vl_risk_children = "vl_risk_group___1",
  vl_risk_poor_people = "vl_risk_group___2",
  vl_risk_forest_rivers = "vl_risk_group___3",
  vl_risk_everyone = "vl_risk_group___4",
  vl_risk_dont_know = "vl_risk_group___5",
  
  list_symptoms_fever = "list_symptoms___1",
  list_symptoms_weight_loss = "list_symptoms___2",
  list_symptoms_fatigue = "list_symptoms___3",
  list_symptoms_appetite_loss = "list_symptoms___4",
  list_symptoms_swelling = "list_symptoms___5",
  list_symptoms_others = "list_symptoms___98"
)


# Report any missing checkbox fields
missing_checkbox <- setdiff(
  unname(checkbox_names),
  names(dat)
)

print(missing_checkbox)

# Copy existing fields to their descriptive names
for (new_name in names(checkbox_names)) {
  
  old_name <- checkbox_names[[new_name]]
  
  if (old_name %in% names(dat)) {
    dat[[new_name]] <- dat[[old_name]]
  }
}


# 08. DESCRIPTIVE ANALYSIS ----------------------------------------------------
levels(dat$transport_affordable)

# A. Primary outcome
primary_summary <- dat %>%
  summarise(
    n_total = n(),
    n_evaluable = sum(!is.na(delayed_care)),
    n_delayed = sum(delayed_care == 1, na.rm = TRUE),
    
    pct_delayed = round(
      100 * mean(delayed_care == 1, na.rm = TRUE), 1
    ),
    
    median_days = median(duration_days, na.rm = TRUE),
    q1 = quantile(duration_days, 0.25, na.rm = TRUE),
    q3 = quantile(duration_days, 0.75, na.rm = TRUE)
  )

print(primary_summary)
save_csv(primary_summary, "06_primary_outcome.csv")


# B. Secondary outcome
referral_summary <- dat %>%
  filter(referred_from_another == 1) %>%
  summarise(
    n_referred = n(),
    n_evaluable = sum(!is.na(delayed_referral)),
    n_delayed = sum(delayed_referral == 1, na.rm = TRUE),
    pct_delayed = if_else(
      n_evaluable > 0,
      round(100 * n_delayed / n_evaluable, 1),
      NA_real_
    )
  )

print(referral_summary)
save_csv(referral_summary, "07_referral_outcome.csv")

colnames(dat)
# C. Table 1: participant characteristics
descriptive_vars <- c(
  "age_cat",
  "health_center",
  "gender",
  "marital_status",
  "kebele_village",
  "residence",
  "educational_status",
  "household_income",
  "advised_delay_formal",
  "family_size_cat",
  "places_sought_care",
  "distance_cat",
  "transport_mode",
  "transport_available",
  "transport_affordable",
  "heard_of_vl"
)

descriptive_vars <- intersect(
  descriptive_vars, names(dat)
)

table1 <- dat %>%
  filter(!is.na(delayed_care)) %>%
  select(
    delayed_care_factor,
    all_of(descriptive_vars)
  ) %>%
  tbl_summary(
    by = delayed_care_factor,
    missing = "ifany"
  ) %>%
  add_overall() %>%
  bold_labels()

table1

save_csv(
  as_tibble(table1, col_labels = TRUE),
  "08_descriptive_table1.csv"
)


# D. Histogram of treatment-seeking duration
p <- ggplot(dat, aes(x = duration_days)) +
  geom_histogram(
    binwidth = 10,
    boundary = 0,
    colour = "white",
    fill = "steelblue"
  ) +
  geom_vline(
    xintercept = 30,
    linetype = "dashed",
    colour = "red"
  ) +
  labs(
    title = "Time from symptom onset to first care",
    x = "Duration (days)",
    y = "Number of participants"
  ) +
  theme_minimal()

print(p)

ggsave(
  file.path(output_dir, "09_duration_histogram.png"),
  p,
  width = 8,
  height = 5,
  dpi = 300
)


# RECOLLAPSING/Merging SPARSE CATEGORIES TO REDUCE SPARSITY ---------------------------


dat <- dat %>%
  mutate(
    # 1. Marital status: Keep Single, Married, Not applicable; merge Divorced/Widowed into Married
    marital_status = case_when(
      marital_status %in% c("Divorced", "Widowed") ~ "Married",
      TRUE ~ as.character(marital_status)
    ),
    marital_status = factor(marital_status, levels = c("Single", "Married", "Not applicable (child)")),
    
    # 2. Residence: Keep Rural, merge Urban and Semi-urban together
    residence = case_when(
      residence %in% c("Urban", "Semi-urban") ~ "Urban / Semi-urban",
      TRUE ~ as.character(residence)
    ),
    residence = factor(residence, levels = c("Rural", "Urban / Semi-urban")),
    
    # 3. Educational status: Merge into No formal education vs. Formal education
    educational_status = case_when(
      educational_status == "No formal education" ~ "No formal education",
      TRUE ~ "Formal education"
    ),
    educational_status = factor(educational_status, levels = c("No formal education", "Formal education")),
    
    # 4. Household income: <= 10000 ETB vs > 10000 ETB
    household_income = case_when(
      household_income == "< 10000 ETB" ~ "<= 10000 ETB",
      TRUE ~ "> 10000 ETB"
    ),
    household_income = factor(household_income, levels = c("<= 10000 ETB", "> 10000 ETB")),
    
    # 5. Family size: <= 6 vs > 6 (using the numeric family_size column)
    family_size_cat = case_when(
      family_size <= 6 ~ "<= 6",
      family_size > 6 ~ "> 6",
      TRUE ~ NA_character_
    ),
    family_size_cat = factor(family_size_cat, levels = c("<= 6", "> 6")),
    
    # 6. Transport mode: On foot vs Public transport / Other
    transport_mode = case_when(
      transport_mode == "On foot" ~ "On foot",
      TRUE ~ "Public transport / Other"
    ),
    transport_mode = factor(transport_mode, levels = c("On foot", "Public transport / Other"))
  )

# 9. Prepare data for BINARY LOGISTIC REGRESSION --------------------------------------------

candidate_vars <- c(
  "gender",
  "health_center",
  "heard_of_vl",
  "vl_preventable",
  "transport_mode",
  "transport_available",
  "transport_affordable",
  "borrow_or_sell_for_treatment",
  "family_stop_work_for_care",
  "belief_sandfly_bite",
  "belief_poor_hygiene",
  "belief_drinking_dirty_water",
  "vl_risk_everyone",
  "age_cat",
  "distance_cat"
)


# Check variables available
missing_candidates <- setdiff(candidate_vars, names(dat))

cat("Missing candidate variables:\n")
print(missing_candidates)

candidate_vars <- intersect(candidate_vars, names(dat))

# Define numeric predictors
numeric_predictors <- c("age_years")

# All other candidate variables are categorical
categorical_predictors <- setdiff(
  candidate_vars,
  numeric_predictors
)

# Regression dataset
analysis <- dat %>%
  filter(!is.na(delayed_care)) %>%
  mutate(
    across(
      all_of(categorical_predictors),
      as.factor
    )
  )

# Exclude predictors with no variation in the analysis sample
candidate_vars <- candidate_vars[
  sapply(candidate_vars, function(v) {
    n_distinct(analysis[[v]], na.rm = TRUE) > 1
  })
]

save_csv(analysis, "10_regression_ready_data.csv")


# 10. Run Binary LOGISTIC REGRESSION analysis -----------------------------------------------

# Each independent variable is fitted separately.
# Outcome: delayed_care (1 = delayed, 0 = not delayed).

# Fit one model per predictor, with error handling
crude_models <- lapply(candidate_vars, function(v) {
  
  model_data <- analysis %>%
    select(delayed_care, all_of(v)) %>%
    drop_na()
  
  tryCatch(
    {
      if (n_distinct(model_data$delayed_care) < 2 ||
          n_distinct(model_data[[v]]) < 2) {
        return(NULL)
      }
      
      glm(
        reformulate(v, response = "delayed_care"),
        data = model_data,
        family = binomial()
      )
    },
    error = function(e) NULL
  )
})

names(crude_models) <- candidate_vars

# Keep successfully fitted models
crude_models <- Filter(Negate(is.null), crude_models)

# Extract crude OR, 95% CI and p-value
crude_results <- purrr::imap_dfr(
  crude_models,
  function(model, variable) {
    
    broom::tidy(
      model,
      exponentiate = TRUE,
      conf.int = TRUE,
      conf.method = "wald"
    ) %>%
      filter(term != "(Intercept)") %>%
      mutate(variable = variable, .before = 1)
  }
)

print(crude_results)
crude_results

save_csv(
  crude_results,
  "11_crude_logistic_regression.csv"
)


# generate crude regression table

crude_tbl <- analysis %>%
  select(
    delayed_care,
    all_of(names(crude_models))
  ) %>%
  tbl_uvregression(
    method = glm,
    y = delayed_care,
    method.args = list(family = binomial),
    exponentiate = TRUE
  ) %>%
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels()

crude_tbl

save_csv(
  as_tibble(crude_tbl, col_labels = TRUE),
  "12_crude_regression_table.csv"
)

# 11. Multivariable LOGISTIC REGRESSION --------------------------------------------

# final used 
adjusted_vars <- c(
  "health_center",
  "transport_mode",
  "borrow_or_sell_for_treatment",
  "age_cat", # because the third category in the crude is significant and pvalue is 0.054
  "distance_cat"
)


# Verify variable availability
missing_adjusted <- setdiff(adjusted_vars, names(analysis))

cat("Missing adjusted variables:\n")
print(missing_adjusted)

# Stop rather than silently fit a different model
if (length(missing_adjusted) > 0) {
  stop("Resolve missing adjusted-model variables before fitting.")
}

# Convert categorical predictors to factors
analysis <- analysis %>%
  mutate(
    across(all_of(adjusted_vars), as.factor)
  )


# Complete cases for adjusted model
adjusted_data <- analysis %>%
  select(delayed_care, all_of(adjusted_vars)) %>%
  drop_na() %>%
  droplevels()

# Check sample size and events
cat("Adjusted model sample:", nrow(adjusted_data), "\n")
print(table(adjusted_data$delayed_care))

# Check whether each factor has at least two levels
factor_levels <- sapply(
  adjusted_data[adjusted_vars],
  function(x) n_distinct(x)
)

print(factor_levels)

if (any(factor_levels < 2)) {
  stop("One or more adjusted predictors have only one observed level.")
}

# Fit adjusted logistic regression
adjusted_formula <- reformulate(
  adjusted_vars,
  response = "delayed_care"
)

adjusted_model <- glm(
  adjusted_formula,
  data = adjusted_data,
  family = binomial(link = "logit")
)

summary(adjusted_model)

# Extract AOR, 95% CI and p-values
adjusted_results <- broom::tidy(
  adjusted_model,
  exponentiate = TRUE,
  conf.int = TRUE,
  conf.method = "wald"
) %>%
  filter(term != "(Intercept)")

print(adjusted_results)

save_csv(
  adjusted_results,
  "13_adjusted_logistic_regression.csv"
)


# generate adjusted regression table
adjusted_tbl <- tbl_regression(
  adjusted_model,
  exponentiate = TRUE
) %>%
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_header(
    estimate = "**Adjusted OR**"
  )

adjusted_tbl

save_csv(
  as_tibble(adjusted_tbl, col_labels = TRUE),
  "14_adjusted_regression_table.csv"
)

# 12. BASIC MODEL DIAGNOSTICS -------------------------------------------------

library(car)

# Calculate Variance Inflation Factors for the adjusted model
vif_values <- vif(adjusted_model)
print(vif_values)

# Save VIF results to the output directory
write.csv(as.data.frame(vif_values), file.path(output_dir, "15_model_vif.csv"))

# Check convergence

cat("Model converged:", adjusted_model$converged, "\n")

# Check extreme fitted probabilities
summary(fitted(adjusted_model))

# Number of model parameters
length(coef(adjusted_model)) - 1

# Check missingness in adjusted predictors
adjusted_missing <- analysis %>%
  select(all_of(adjusted_vars)) %>%
  summarise(
    across(
      everything(),
      ~sum(is.na(.))
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_missing"
  )

save_csv(
  adjusted_missing,
  "15_adjusted_model_missingness.csv"
)

#Plot ROC-AUC curve to assess model discrimination 

# Load required libraries for performance evaluation
library(pROC)
library(ResourceSelection)

# 1. Calculate AUC value
roc_obj <- roc(adjusted_model$y, fitted(adjusted_model))
auc_value <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

cat("AUC:", round(auc_value, 3), "\n")
cat("95% CI:", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), "\n")

# 2. Hosmer-Lemeshow Goodness-of-Fit Test (p > 0.05 indicates good fit)
hl_test <- hoslem.test(adjusted_model$y, fitted(adjusted_model), g = 10)
print(hl_test)

# 3. Generate a Publication-Ready ROC Curve Plot
png(file.path(output_dir, "16_roc_curve_plot.png"), width = 800, height = 800, res = 120)
ROC_curve <- plot(
  roc_obj, 
  col = "#2b5c8f", 
  lwd = 3, 
  main = "ROC Curve for Delayed Care-Seeking Model",
  print.auc = TRUE, 
  auc.polygon = TRUE, 
  auc.polygon.col = "aliceblue",
  grid = c(0.2, 0.2), 
  grid.col = "gray80"
)
dev.off()

plot(ROC_curve)

#======================================

#PLOTS

# -----------------------------------------------------------------------------
# ROC CURVE WITH YOUDEN'S INDEX
# -----------------------------------------------------------------------------
library(pROC)
library(ggplot2)

# 1. Calculate ROC object and confidence intervals
roc_obj <- roc(adjusted_model$y, fitted(adjusted_model))
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

# 2. Extract coordinates for ggplot2
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities,
  thresholds = roc_obj$thresholds
)

# 3. Find optimal threshold using Youden's Index (Sensitivity + Specificity - 1)
coords_best <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
opt_fpr <- 1 - coords_best$specificity
opt_tpr <- coords_best$sensitivity

# 4. Generate Publication-Ready ROC Plot
p_roc_pub <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  # Diagonal reference line (Chance line)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  # ROC Curve line
  geom_line(color = "#2b5c8f", linewidth = 1.2) +
  # Shaded area under the curve
  geom_ribbon(aes(ymin = 0, ymax = TPR), fill = "#2b5c8f", alpha = 0.1) +
  # Optimal cutoff point (Youden's Index)
  geom_point(aes(x = opt_fpr, y = opt_tpr), color = "#d95f02", size = 3.5) +
  geom_text(aes(x = opt_fpr, y = opt_tpr, label = paste0("  Optimal Cutoff\n  (Sens: ", round(opt_tpr, 2), ", Spec: ", round(coords_best$specificity, 2), ")")),
            hjust = 0, vjust = 0, size = 3.5, color = "#d95f02", fontface = "bold") +
  # Labels and subtitle with AUC and 95% CI
  labs(
    title = "ROC Curve for Multivariable Delayed Care-Seeking Model",
    subtitle = sprintf("Area Under the Curve (AUC) = %.3f (95%% CI: %.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  coord_fixed() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0.5)
  )

print(p_roc_pub)

# 5. Save high-resolution publication image
ggsave(
  file.path(output_dir, "16_roc_curve_publication.png"),
  p_roc_pub,
  width = 7,
  height = 7,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Forest Plot of Adjusted Odds Ratios
# -----------------------------------------------------------------------------
library(ggstats)
library(ggplot2)

p_forest <- ggcoef_model(
  adjusted_model,
  exponentiate = TRUE,
  variable_labels = c(
    health_center = "Health Center",
    transport_mode = "Transport Mode",
    borrow_or_sell_for_treatment = "Sold Assets",
    age_cat = "Age",
    distance_cat = "Distance"
  ),
  show_p_values = TRUE,
  signif_stars = FALSE
) +
  labs(
    title = "Multivariable Odds Ratios for Delayed Treatment-Seeking Among VL Patients in South Ethiopia",
    x = "Adjusted Odds Ratio (95% CI)"
  ) +
 theme_bw(base_size = 12) +
          theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

p_forest

ggsave(file.path(output_dir, "17_forest_plot_aor.png"), p_forest, width = 8, height = 5, dpi = 300)


# -----------------------------------------------------------------------------
# Predicted Probability Density Plot
# -----------------------------------------------------------------------------

adjusted_data$pred_prob <- fitted(adjusted_model)

p_prob <- ggplot(adjusted_data, aes(x = pred_prob, fill = factor(delayed_care, labels = c("Not Delayed", "Delayed")))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("magenta", "blue4"), name = "Treatment-seeking Status") +
  labs(
    title = "Distribution of Predicted Risk Probabilities: Delayed Treatment-Seeking Among VL-Patients in South Ethiopia",
    x = "Predicted Probability of Delay",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
)

p_prob
ggsave(file.path(output_dir, "19_predicted_probability_density.png"), p_prob, width = 7, height = 4.5, dpi = 300)


# -----------------------------------------------------------------------------
# Proportion of Treatment Delay by Distance and Transport Mode
# -----------------------------------------------------------------------------
p_prop <- dat %>%
  filter(!is.na(delayed_care_factor), !is.na(distance_cat), !is.na(transport_mode)) %>%
  group_by(distance_cat, transport_mode, delayed_care_factor) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(distance_cat, transport_mode) %>%
  mutate(pct = n / sum(n) * 100) %>%
  filter(delayed_care_factor == "Delayed") %>%
  ggplot(aes(x = distance_cat, y = pct, fill = transport_mode)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("#4682B4", "#D95F02"), name = "Transport Mode") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Percentage of vL Patients Experiencing Delay (>=30 days)",
    x = "Distance to Health Facility",
    y = "Delayed Care Seeking (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", 
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5))


ggsave(file.path(output_dir, "21_delay_proportion_stratified.png"), p_prop, width = 7.5, height = 5, dpi = 300)

p_prop


# -----------------------------------------------------------------------------
# Proportion of treatment seeking across heath centers
# -----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggalluvial)

p_hc_alluvial <- dat %>%
  filter(!is.na(health_center), !is.na(delayed_care_factor)) %>%
  count(health_center, delayed_care_factor) %>%
  ggplot(aes(y = n, axis1 = health_center, axis2 = delayed_care_factor)) +
  # Flowing ribbons connecting centers to status
  geom_alluvium(aes(fill = delayed_care_factor), width = 1/8, alpha = 0.85, curve_type = "sigmoid") +
  # Vertical blocks (strata) for categories
  geom_stratum(width = 1/5, fill = "grey2", color = "white", linewidth = 0.5) +
  # Labels inside the stratum blocks
  geom_text(
    stat = "stratum", 
    aes(label = after_stat(stratum)), 
    color = "white", 
    size = 3.5, 
    fontface = "bold"
  ) +
  scale_x_discrete(limits = c("Health Center", "Care Status"), expand = c(0.15, 0.15)) +
  # Custom palette: Soft Coral for Delayed, Deep Teal/Cyan for On-Time
  scale_fill_manual(
    values = c("red4", "blue4"),  #"#E76F51", "#2A9D8F")
    name = "Care Status"
  ) +
  labs(
    title = "Proportion of Delayed Treatment-Seeking Across Health Centers",
    subtitle = "Patient flow breakdown by facility catchment area and care seeking status",
    x = NULL,
    y = "Number of Patients"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    # Centering title and subtitle (hjust = 0.5)
    plot.title = element_text(face = "bold", size = 16, color = "#1D3557", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 11, color = "#6C757D", hjust = 0.5, margin = margin(b = 20)),
    
    # Axis & Panel styling
    axis.text.x = element_text(face = "bold", size = 12, color = "#1D3557"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    
    # Legend position
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_text(face = "bold", size = 10),
    
    # Outer margins for breathing room
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25)
  )

# Save high-resolution plot
ggsave(
  file.path(output_dir, "22_health_center_delay_alluvial.png"), 
  p_hc_alluvial, 
  width = 8.5, 
  height = 6, 
  dpi = 300
)

# Display plot
p_hc_alluvial

# -----------------------------------------------------------------------------
# CALIBRATION PLOT (Logistic Smoothing)
# -----------------------------------------------------------------------------

library(ResourceSelection)
library(ggplot2)

# 1. Prepare data
calib_data <- adjusted_data %>%
  mutate(
    pred_prob = fitted(adjusted_model),
    outcome = as.numeric(delayed_care)
  )

hl <- hoslem.test(calib_data$outcome, calib_data$pred_prob, g = 5) # g=5 or g=6 better for N=92 to avoid empty bins
hl
# 2. Group into quintiles (g=5) or deciles for summary points
calib_data <- calib_data %>%
  mutate(risk_group = ntile(pred_prob, 5))

calib_summary <- calib_data %>%
  group_by(risk_group) %>%
  summarise(
    mean_pred = mean(pred_prob),
    obs_rate = mean(outcome),
    n = n(),
    .groups = "drop"
  )

# 3. Generate clean calibration plot with logistic calibration fit
p_calib_fixed <- ggplot(calib_summary, aes(x = mean_pred, y = obs_rate)) +
  # Perfect calibration reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  # Logistic calibration curve from raw data
  geom_smooth(data = calib_data, aes(x = pred_prob, y = outcome), 
              method = "glm", method.args = list(family = "binomial"), # when the method was loess bad 'geom_smooth(method = "loess",'
              se = TRUE, color = "#2b5c8f", fill = "aliceblue", linewidth = 1.2) +
  # Grouped points (sized by group N)
  geom_point(aes(size = n), color = "#d95f02", alpha = 0.8) +
  scale_size_continuous(range = c(3, 7), guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Calibration Plot for Delayed Treatment-Seeking for VL Model",
    subtitle = paste0("Hosmer-Lemeshow p-value = ", round(hl$p.value, 3)),
    x = "Predicted Probability of Delay",
    y = "Observed Proportion of Delay"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  )

print(p_calib_fixed)

ggsave(file.path(output_dir, "17_calibration_plot_fixed.png"), p_calib_fixed, width = 7, height = 6, dpi = 300)
# summary of the calibration across quintiles. 
View(calib_summary)
print(calib_summary, n = Inf)


#----------------------------
# SAVE Important documents and EVERYTHING in the global env
# save data and models
saveRDS(
  dat,
  file.path(output_dir, "CORAL_clean_data.rds")
)

save_csv(
  dat,
  "CORAL_clean_data.csv"
)

saveRDS(
  adjusted_model,
  file.path(output_dir, "adjusted_model.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt")
)

cat("\nAnalysis completed. Results saved in:", output_dir, "\n")


#save.image(file = "Delayed Rx seeking_All quantitative Objects.RDATA")