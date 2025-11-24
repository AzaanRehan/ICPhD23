meth_1114 <- read.delim("new_meth_1114.txt")
colnames(meth_1114)

pollution_cols <- grep("^M(0[1-9]|1[0-2])_(NO2|O3|PM10|PM25)\\.20(08|09|10|11|12)$",
                       names(meth_1114),
                       value = TRUE)

length(pollution_cols)
head(pollution_cols)
tail(pollution_cols)

library(tidyverse)
library(lubridate)

# ---------- 0. defensive ID detection ----------
id_candidates <- c("SID", "part_id", "BARCODE", "id")
id_col <- id_candidates[id_candidates %in% names(meth_1114)][1]
if(is.na(id_col)) {
  stop("No candidate ID column found. Please set 'id_col' manually to the name of your participant ID column.")
} else message("Using ID column: ", id_col)

# ---------- 1. pollution columns already detected earlier ----------
pollution_cols <- grep("^M(0[1-9]|1[0-2])_(NO2|O3|PM10|PM25)\\.20(08|09|10|11|12)$",
                       names(meth_1114),
                       value = TRUE)

# sanity
if(length(pollution_cols) != 240) {
  warning("Expected 240 pollution cols; found ", length(pollution_cols), ". Proceeding anyway.")
}

# ---------- 2. pivot to long ----------
pollution_long <- meth_1114 %>%
  select(all_of(id_col), all_of(pollution_cols), everything()) %>% 
  pivot_longer(
    cols = all_of(pollution_cols),
    names_to = c("month", "pollutant", "year"),
    names_pattern = "M(0[1-9]|1[0-2])_([A-Za-z0-9]+)\\.(20[0-9]{2})",
    values_to = "value",
    values_drop_na = FALSE
  ) %>%
  mutate(
    month = as.integer(month),
    year  = as.integer(year),
    date  = as.Date(sprintf("%04d-%02d-01", year, month)),
    pollutant = as.character(pollutant)
  ) %>%
  # keep a tidy column order for clarity
  relocate(all_of(id_col), date, year, month, pollutant, value)

# ---------- 3. quick integrity checks ----------
n_rows <- nrow(pollution_long)
n_unique_ids <- n_distinct(pollution_long[[id_col]])
n_pollutant_types <- length(unique(pollution_long$pollutant))
months_present <- sort(unique(pollution_long$date))

message("Long-format created.")
message("Rows: ", n_rows)
message("Unique IDs: ", n_unique_ids)
message("Pollutants found: ", paste(sort(unique(pollution_long$pollutant)), collapse = ", "))
message("Date range (first 3): ", paste(head(months_present,3), collapse = ", "), " ... (last 3): ", paste(tail(months_present,3), collapse = ", "))
message("Rows per person expected (if complete): ", length(pollution_cols))
message("Typical rows per person (median): ", median(table(pollution_long[[id_col]])))

meth_1114 <- meth_1114 %>%
  mutate(screen_year = lubridate::year(first_clinic_date))

pollution_annual <- pollution_long %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(SID, pollutant, year) %>%
  summarise(
    annual_avg = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

meth_pollution <- meth_1114 %>%
  mutate(screen_year = year(first_clinic_date)) %>%
  left_join(
    pollution_annual,
    by = c("SID" = "SID", "screen_year" = "year")
  )

meth_pollution_wide <- meth_pollution %>%
  pivot_wider(
    names_from = pollutant,
    values_from = annual_avg,
    names_prefix = "",
    names_sep = "_year"
  )

meth_pollution_wide_filtered <- meth_pollution_wide %>%
  filter(screen_year >= 2008)

pollutant_IQRs <- meth_pollution_wide_filtered %>%
  summarise(
    IQR_NO2 = IQR(NO2, na.rm = TRUE),
    IQR_O3 = IQR(O3, na.rm = TRUE),
    IQR_PM10 = IQR(PM10, na.rm = TRUE),
    IQR_PM25 = IQR(PM25, na.rm = TRUE)
  )
print(pollutant_IQRs)


vars_to_factor <- c("gender", "smoking3", "alc_code", "hhincome_code",
                    "edu_code", "phy_code_chung", "anx_code")

meth_pollution_wide_filtered <- meth_pollution_wide_filtered %>%
  mutate(across(all_of(vars_to_factor), as.factor))

meth_pollution_wide_filtered <- meth_pollution_wide_filtered %>%
  mutate(
    total_fruit_veg_portions_per_week = as.numeric(total_fruit_veg_portions_per_week)
  )

# MAIN MODELS with 4 pollutants and 5 clocks
clocks_to_scale <- c("AgeAccelPheno", "AgeAccelGrim", "bAgeAccel", "epiTOC1", "DunedinPACE")

meth_pollution_scaled <- meth_pollution_wide_filtered %>%
  mutate(across(all_of(clocks_to_scale), ~ scale(.)[, 1], .names = "{.col}_z"))

library(dplyr)
library(broom)
library(purrr)
library(tidyr)

# 1. Define z-scored clocks and pollutants
clock_vars <- c("AgeAccelPheno_z", "AgeAccelGrim_z", "bAgeAccel_z", "epiTOC1_z", "DunedinPACE_z")
pollutants <- c("NO2", "O3", "PM10", "PM25")

# 2. Covariates for adjustment
covariates <- c("age_when_screened", "gender", "body_mass_index", "smoking3", "alc_code",
                "hhincome_code", "edu_code", "phy_code_chung", "anx_code", "total_fruit_veg_portions_per_week")

# 3. Function to fit a model for one clock-pollutant pair
fit_model <- function(clock, pollutant) {
  formula <- as.formula(
    paste(clock, "~", pollutant, "+", paste(covariates, collapse = " + "))
  )
  
  lm(formula, data = meth_pollution_scaled) %>%
    tidy(conf.int = TRUE) %>%
    filter(term == pollutant)  # no need to add clock/pollutant again here
}


# 4. Run models across all clock × pollutant combinations
model_results <- cross_df(list(clock = clock_vars, pollutant = pollutants)) %>%
  mutate(results = map2(clock, pollutant, fit_model)) %>%
  unnest(results)

# Create a named vector of IQRs
iqr_lookup <- c(
  NO2 = pollutant_IQRs$IQR_NO2,
  O3 = pollutant_IQRs$IQR_O3,
  PM10 = pollutant_IQRs$IQR_PM10,
  PM25 = pollutant_IQRs$IQR_PM25
)

# Scale results
model_results_scaled <- model_results %>%
  mutate(
    IQR = iqr_lookup[pollutant],
    estimate_IQR = estimate * IQR,
    conf.low_IQR = conf.low * IQR,
    conf.high_IQR = conf.high * IQR
  )

model_results_scaled <- model_results_scaled %>%
  mutate(clock_clean = recode(clock,
                              AgeAccelPheno_z = "PhenoAA",
                              AgeAccelGrim_z  = "GrimAA",
                              bAgeAccel_z     = "bernabeuAA",
                              epiTOC1_z       = "Epitoc",
                              DunedinPACE_z   = "DundeinPACE"
  ))

library(ggplot2)

# Flip the levels to have PhenoAA at the top
clock_levels <- c("PhenoAA", "GrimAA", "bernabeuAA", "Epitoc", "DundeinPACE")

model_results_scaled <- model_results_scaled %>%
  mutate(clock_clean = factor(clock_clean, levels = rev(clock_levels)))


# Plot per pollutant
plots <- model_results_scaled %>%
  split(.$pollutant) %>%
  lapply(function(df) {
    ggplot(df, aes(x = estimate_IQR, y = clock_clean)) +
      geom_point() +
      geom_errorbar(
        aes(xmin = conf.low_IQR, xmax = conf.high_IQR),
        height = 0.2,
        orientation = "y"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        x = "SD change in EAA per IQR increase in pollutant",
        y = NULL,
        title = paste0("Annual Average ", df$pollutant[1], " Exposure & Age Acceleration")
      ) +
      theme_minimal()
  })

plots$NO2
plots$PM25
plots$PM10
plots$O3

# Display effect estimates per IQR increase with confidence intervals and p-values
library(dplyr)

model_results_scaled %>%
  select(
    Clock = clock_clean,
    Pollutant = pollutant,
    Estimate_per_IQR = estimate_IQR,
    CI_Lower = conf.low_IQR,
    CI_Upper = conf.high_IQR,
    P_Value = p.value
  ) %>%
  arrange(Pollutant, factor(Clock, levels = c("PhenoAA", "GrimAA", "bernabeuAA", "Epitoc", "DundeinPACE")))


# NEW MITOTIC CLOCK MODELS and PM2.5 & PM10

library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(ggplot2)

# ---- 1. Define new mitotic clocks ----
new_mitotic_clocks <- c("stemTOC", "RepliTali", "HypoClock")

# ---- 2. Standardize them ----
meth_pollution_scaled <- meth_pollution_scaled %>%
  mutate(across(all_of(new_mitotic_clocks), ~ scale(.)[, 1], .names = "{.col}_z"))

# ---- 3. Define pollutants of interest ----
pollutants_new <- c("PM10", "PM25")

# ---- 4. Adjustment set (same as before) ----
covariates <- c(
  "age_when_screened", "gender", "body_mass_index", "smoking3", "alc_code",
  "hhincome_code", "edu_code", "phy_code_chung", "anx_code",
  "total_fruit_veg_portions_per_week"
)

# ---- 5. Model fitting function ----
fit_model <- function(clock, pollutant) {
  formula <- as.formula(
    paste(clock, "~", pollutant, "+", paste(covariates, collapse = " + "))
  )
  lm(formula, data = meth_pollution_scaled) %>%
    tidy(conf.int = TRUE) %>%
    filter(term == pollutant) %>%
    mutate(clock = clock, pollutant = pollutant)
}

# ---- 6. Run models (3 clocks × 2 pollutants = 6 models) ----
clock_vars_new <- paste0(new_mitotic_clocks, "_z")

model_results_new <- expand_grid(clock = clock_vars_new, pollutant = pollutants_new) %>%
  mutate(results = map2(clock, pollutant, fit_model)) %>%
  unnest(results, names_sep = "_")

# ---- 7. Apply IQR scaling for PM10 and PM25 ----
iqr_lookup <- c(PM10 = pollutant_IQRs$IQR_PM10, PM25 = pollutant_IQRs$IQR_PM25)

model_results_scaled_new <- model_results_new %>%
  mutate(
    IQR = iqr_lookup[pollutant],
    estimate_IQR = results_estimate * IQR,
    conf.low_IQR = results_conf.low * IQR,
    conf.high_IQR = results_conf.high * IQR
  )

# ---- 8. Clean up clock names for plotting ----
clock_levels_new <- c("stemTOC", "RepliTali", "HypoClock")

model_results_scaled_new <- model_results_scaled_new %>%
  mutate(
    clock_clean = recode(clock,
                         stemTOC_z = "stemTOC",
                         RepliTali_z = "RepliTali",
                         HypoClock_z = "HypoClock"),
    clock_clean = factor(clock_clean, levels = rev(clock_levels_new))
  )

# ---- 9. Create plots ----
plots_new <- model_results_scaled_new %>%
  split(.$pollutant) %>%
  lapply(function(df) {
    ggplot(df, aes(x = estimate_IQR, y = clock_clean)) +
      geom_point() +
      geom_errorbar(
        aes(xmin = conf.low_IQR, xmax = conf.high_IQR),
        height = 0.2,
        orientation = "y"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        x = "SD change in Mitotic EAA per IQR increase in pollutant",
        y = NULL,
        title = paste0("Annual Average ", df$pollutant[1], " Exposure & Age Acceleration")
      ) +
      scale_x_continuous(breaks = seq(0, 0.3, by = 0.1)) +   # adjust tick spacing here
      theme_minimal()
  })

# ---- 10. View plots ----
plots_new$PM10
plots_new$PM25

# ---- 11. View results ----
model_results_scaled_new %>%
  select(
    Clock = clock_clean,
    Pollutant = pollutant,
    Estimate_per_IQR = estimate_IQR,
    CI_Lower = conf.low_IQR,
    CI_Upper = conf.high_IQR,
    P_Value = results_p.value
  ) %>%
  arrange(Pollutant, factor(Clock, levels = c("stemTOC", "RepliTali", "HypoClock")))

# Apply FDR correction for multiple testing
model_results_scaled_new <- model_results_scaled_new %>%
  mutate(FDR_Adjusted_P = p.adjust(results_p.value, method = "fdr"))

# View results with both raw and adjusted p-values
model_results_scaled_new %>%
  select(
    Clock = clock_clean,
    Pollutant = pollutant,
    Estimate_per_IQR = estimate_IQR,
    CI_Lower = conf.low_IQR,
    CI_Upper = conf.high_IQR,
    P_Value = results_p.value,
    FDR_Adjusted_P
  )

# MAIN MODELS with cell-type adjustment
# --- Define covariates including cell proportions ---
covariates_cell <- c(
  "age_when_screened", "gender", "body_mass_index", "smoking3", "alc_code",
  "hhincome_code", "edu_code", "phy_code_chung", "anx_code",
  "total_fruit_veg_portions_per_week",
  "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "PlasmaBlast"
)

# --- Fit models (using previously scaled dataset) ---
fit_model_cell <- function(clock, pollutant) {
  lm(as.formula(paste(clock, "~", pollutant, "+", paste(covariates_cell, collapse = " + "))),
     data = meth_pollution_scaled) %>%
    tidy(conf.int = TRUE) %>%
    filter(term == pollutant)
}

# --- Run models across 5 clocks × 4 pollutants ---
model_results_cell <- expand_grid(
  clock = paste0(c("AgeAccelPheno", "AgeAccelGrim", "bAgeAccel", "epiTOC1", "DunedinPACE"), "_z"),
  pollutant = c("NO2", "O3", "PM10", "PM25")
) %>%
  mutate(results = map2(clock, pollutant, fit_model_cell)) %>%
  unnest(results, names_sep = "_")

# --- Apply IQR scaling safely ---
model_results_scaled_cell <- model_results_cell %>%
  mutate(
    IQR = as.numeric(unlist(pollutant_IQRs[paste0("IQR_", pollutant)])),
    estimate_IQR = results_estimate * IQR,
    conf.low_IQR = results_conf.low * IQR,
    conf.high_IQR = results_conf.high * IQR,
    clock_clean = recode(clock,
                         AgeAccelPheno_z = "PhenoAA",
                         AgeAccelGrim_z  = "GrimAA",
                         bAgeAccel_z     = "bernabeuAA",
                         epiTOC1_z       = "Epitoc",
                         DunedinPACE_z   = "DundeinPACE")
  )

# --- Plot results ---
clock_levels <- c("PhenoAA", "GrimAA", "bernabeuAA", "Epitoc", "DundeinPACE")

plots_cell <- model_results_scaled_cell %>%
  mutate(clock_clean = factor(clock_clean, levels = rev(clock_levels))) %>%
  split(.$pollutant) %>%
  map(~ ggplot(.x, aes(x = estimate_IQR, y = clock_clean)) +
        geom_point() +
        geom_errorbar(aes(xmin = conf.low_IQR, xmax = conf.high_IQR), height = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(
          x = "SD change in EAA per IQR increase in pollutant (cell-type adjusted)",
          y = NULL,
          title = paste0("Annual Average ", .x$pollutant[1], " Exposure & Age Acceleration")
        ) + 
        theme_minimal()+
        xlim(-0.3,0.3))

# --- View plots ---
plots_cell$NO2
plots_cell$O3
plots_cell$PM10
plots_cell$PM25

# ---- View cell-type adjusted results with FDR correction ----
model_results_scaled_cell <- model_results_scaled_cell %>%
  mutate(FDR_Adjusted_P = p.adjust(results_p.value, method = "fdr"))

# Print results for manual table creation
model_results_scaled_cell %>%
  select(
    Clock = clock_clean,
    Pollutant = pollutant,
    Estimate_per_IQR = estimate_IQR,
    CI_Lower = conf.low_IQR,
    CI_Upper = conf.high_IQR,
    P_Value = results_p.value,
    FDR_Adjusted_P
  ) %>%
  arrange(Pollutant, factor(Clock, levels = c("PhenoAA", "GrimAA", "bernabeuAA", "Epitoc", "DundeinPACE")))

# Cell-type adjusted mitotic clock models (PM10 & PM25)
new_mitotic_clocks <- c("stemTOC", "RepliTali", "HypoClock")
pollutants_new <- c("PM10", "PM25")
covariates_cell <- c("age_when_screened","gender","body_mass_index","smoking3","alc_code",
                     "hhincome_code","edu_code","phy_code_chung","anx_code","total_fruit_veg_portions_per_week",
                     "CD8T","CD4T","NK","Bcell","Mono","Gran","PlasmaBlast")

meth_pollution_scaled <- meth_pollution_scaled %>%
  mutate(across(all_of(new_mitotic_clocks), ~ scale(.)[,1], .names="{.col}_z"))

fit_model_mitotic_cell <- function(clock, pollutant)
  lm(as.formula(paste(clock, "~", pollutant, "+", paste(covariates_cell, collapse=" + "))),
     data=meth_pollution_scaled) %>%
  tidy(conf.int=TRUE) %>% filter(term==pollutant) %>%
  mutate(clock=clock, pollutant=pollutant)

model_results_mitotic_scaled_cell <- expand_grid(clock=paste0(new_mitotic_clocks,"_z"),
                                                 pollutant=pollutants_new) %>%
  mutate(results=map2(clock,pollutant,fit_model_mitotic_cell)) %>%
  unnest(results, names_sep="_") %>%
  mutate(IQR=c(PM10=pollutant_IQRs$IQR_PM10, PM25=pollutant_IQRs$IQR_PM25)[pollutant],
         estimate_IQR=results_estimate*IQR,
         conf.low_IQR=results_conf.low*IQR,
         conf.high_IQR=results_conf.high*IQR,
         clock_clean=recode(clock, stemTOC_z="stemTOC", RepliTali_z="RepliTali", HypoClock_z="HypoClock"),
         clock_clean=factor(clock_clean, levels=rev(c("stemTOC","RepliTali","HypoClock"))),
         FDR_Adjusted_P=p.adjust(results_p.value, method="fdr"))

plots_mitotic_cell <- model_results_mitotic_scaled_cell %>%
  split(.$pollutant) %>%
  map(~ ggplot(.x, aes(x=estimate_IQR, y=clock_clean)) +
        geom_point() +
        geom_errorbar(aes(xmin=conf.low_IQR, xmax=conf.high_IQR), height=0.2) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        scale_x_continuous(breaks=seq(0,0.3,by=0.1)) +
        labs(x="SD change in Mitotic EAA per IQR increase in pollutant (cell-type adjusted)",
             y=NULL, title=paste0("Annual Average ", .x$pollutant[1], " Exposure & Mitotic Clocks")) +
        theme_minimal())

plots_mitotic_cell$PM10
plots_mitotic_cell$PM25

model_results_mitotic_scaled_cell %>%
  select(Clock=clock_clean, Pollutant=pollutant, Estimate_per_IQR=estimate_IQR,
         CI_Lower=conf.low_IQR, CI_Upper=conf.high_IQR, P_Value=results_p.value,
         FDR_Adjusted_P)

