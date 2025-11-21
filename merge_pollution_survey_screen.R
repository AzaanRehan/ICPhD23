# merge datasets
setwd("~/OneDrive - Imperial College London/PhD/Project work/AIRWAVE/1129_meth")

# pollution
pollution_score <- read.csv("POLLUTION_SCORES_29-MAY.CSV")
pollution_part <- read.csv("POLLUTION_PARTICIPANTS_29-MAY.CSV")

pollution_score_wide <- reshape(
  pollution_score,
  timevar = "YEAR_VALID",
  idvar = "POSTCODE_X",
  direction = "wide"
)
summary(duplicated(pollution_score_wide$POSTCODE_X))
pollution_score_wide$BARCODE <- pollution_part$BARCODE[match(pollution_score_wide$POSTCODE_X, pollution_part$POSTCODE_X)]

#head(pollution_score_wide)

#rm(pollution_score)
#rm(pollution_part)

## clinical survey v2.0
survey <- read.csv("Airwave_clinic_surveyV02.csv")
colnames(survey)[2] <- "BARCODE"
merge_dat <- merge(pollution_score_wide, survey, by = "BARCODE")

## screen
screen <- read.delim("airwave-screening-v3.tsv")
colnames(screen)
merge_dat <- merge(merge_dat, screen, by = "BARCODE")

write.table(merge_dat,"updated_whole_pol.txt", sep = "\t", quote = F, row.names = F)

# Read the merged methylation data
merge_dat <- read.delim("survey_screen_pol.txt")

meth_1114 <- read.delim("meth_1114.txt")

colnames(sample_mapping)
colnames(meth_1114)

# Make barcodes same case to avoid mismatches
sample_mapping$barcode <- toupper(sample_mapping$barcode)
meth_1114$BARCODE      <- toupper(meth_1114$BARCODE)

# Merge SID into meth_1114
meth_1114 <- merge(
  meth_1114,
  sample_mapping[, c("barcode", "sample_id")],
  by.x = "BARCODE",
  by.y = "barcode",
  all.x = TRUE
)

names(meth_1114)[names(meth_1114) == "sample_id"] <- "SID"


# Merge clock results (left-joins, preserve order)
meth_1114 <- merge(meth_1114, Mit_clocks,     by = "SID", all.x = TRUE, sort = FALSE)
meth_1114 <- merge(meth_1114, DunedinPACE,    by = "SID", all.x = TRUE, sort = FALSE)
meth_1114 <- merge(meth_1114, bernabeu,       by = "SID", all.x = TRUE, sort = FALSE)

# Only strip a leading 'X' from calculator_res$SID if you've confirmed it's spurious; otherwise skip.
meth_1114 <- merge(meth_1114, calculator_res, by = "SID", all.x = TRUE, sort = FALSE)

# Restore original sample order and drop helper column
#meth_1129 <- meth_1129[order(meth_1129$.ord), ]
#meth_1129$.ord <- NULL

# Write out the final table as tab-delimited text
write.table(meth_1114, file = "new_meth_1114.txt", sep = "\t", quote = FALSE, row.names = FALSE)

update_merge <- read.delim("new_meth_1114.txt")

# smoking
update_merge$smoking3[update_merge$is_smoker == "NO"]<-"never smoker"
update_merge$smoking3[update_merge$f_220 == "Yes"]<- "former smoker"
update_merge$smoking3[update_merge$is_smoker == "YES"]<-"current smoker"
update_merge$smoking3<-as.factor(update_merge$smoking3)
table(update_merge$smoking3)

# alcohol
update_merge$drinksT <- update_merge$f_459 + update_merge$f_460 + update_merge$f_461 + update_merge$f_462 + update_merge$f_463
update_merge$drinksT[update_merge$drinksT == -30] <- NA
update_merge$alc_best_3 <- NA
update_merge$alc_best_3[update_merge$drinksT == 0] <- 0
update_merge$alc_best_3[update_merge$drinksT > 0 & update_merge$drinksT <= 14] <- 1
update_merge$alc_best_3[update_merge$drinksT > 14] <- 2
table(update_merge$alc_best_3, useNA = "ifany")

# Alcohol consumption: Occasional Drinker: Consumes alcohol occasionally but less than the defined threshold ≤ 14 alcohol units/week for women and ≤ 21 alcohol units/week for men (alternative way to define is Less than 28 grams/day (just for model fit purposes maybe). Habitual Drinker: Consumes alcohol regularly or more than the defined threshold
update_merge$sex <- ifelse(update_merge$gender == 0, "M", "F")
occ_female <- update_merge$sex == "F" & update_merge$drinksT > 0 & update_merge$drinksT <= 14
occ_male <- update_merge$sex == "M" & update_merge$drinksT > 0 & update_merge$drinksT <= 21
update_merge$alc_code <- NA
update_merge$alc_code[update_merge$drinksT == 0] <- "Never Drinker"
update_merge$alc_code[occ_female | occ_male] <- "Occasional"
update_merge$alc_code[update_merge$drinksT > 21] <- "Habitual"
update_merge$alc_code[is.na(update_merge$drinksT)] <- NA
table(update_merge$alc_code, useNA = "ifany")


# education
update_merge$edu_code[update_merge$f_003=="GSCE/O-Level/CSE"|update_merge$f_003=="Left school before taking O levels / GCSEs"|update_merge$f_003=="Vocational qualifications (NVQ1+2)"] <- "Low"
update_merge$edu_code[update_merge$f_003=="A levels / Highers or equivalent (NVQ3)"] <- "Medium"
update_merge$edu_code[update_merge$f_003=="Bachelor Degree or equivalent (NVQ4)"|update_merge$f_003=="Postgraduate qualifications"] <- "High"
table(update_merge$edu_code)

# BMI
update_merge$body_mass_index[update_merge$body_mass_index<0] <- NA
mean(update_merge$body_mass_index)

# Household income: Low: Income less than 38000, Medium: Income between 38000 and 77999, High: Income more than 78000
update_merge$hhincome_code[update_merge$income_household=="A) Less than 25999"|update_merge$income_household=="B) 26000 - 37999"] <- "Low"
update_merge$hhincome_code[update_merge$income_household=="C) 38000 - 57999"|update_merge$income_household=="D) 58000 - 77999"] <- "Medium"
update_merge$hhincome_code[update_merge$income_household=="E) More than  78000"] <- "High"
table(update_merge$hhincome_code)



# Physical activity: Low: Does not meet criteria for low or high activity; Moderate: 5 or more days of moderate-intensity activity or walking of at least 30 minutes per day; High: 7 or more days of any combination of walking, moderate-intensity, or vigorous-intensity activities achieving a minimum of at least 3000 MET-minutes/week.
# define moderate
con1_m <- update_merge$f_245>=30|update_merge$f_247>=600
update_merge$phy_code[update_merge$f_243>5 & con1_m] <- "Moderate"
# define high
## Define MET values
vigorous_met <- 8
moderate_met <- 4
walking_met <- 3.3
## calculate MET per week
vigorous_met_min <- ifelse(update_merge$f_241 == 0, update_merge$f_237 * update_merge$f_239 * vigorous_met, update_merge$f_241 * vigorous_met)
moderate_met_min <- ifelse(update_merge$f_247 == 0, update_merge$f_243 * update_merge$f_245 * moderate_met, update_merge$f_247 * moderate_met)
walking_met_min <- update_merge$f_251 * walking_met
vigorous_met_min[vigorous_met_min<0] <- NA
moderate_met_min[moderate_met_min<0] <- NA
walking_met_min[walking_met_min<0] <- NA
total_met_min_week <- vigorous_met_min + moderate_met_min + walking_met_min
update_merge$phy_code[vigorous_met_min>7000|moderate_met_min>7000|walking_met_min>7000|total_met_min_week>7000] <- "High"
all_na_rows <- apply(update_merge[,c("f_243", "f_247", "f_245", "f_237", "f_239", "f_241", "f_251")], 1, function(x) all(is.na(x)))
summary(all_na_rows) # no all na rows
update_merge$phy_code[is.na(update_merge$phy_code)] <- "Low"
table(update_merge$phy_code)

## *** CHUNGHO's METHOD *** ##
update_merge$totmets <- 8 *update_merge$f_237 * update_merge$f_239  + 4 *update_merge$f_243 * update_merge$f_245 + 3.3 * update_merge$f_251
update_merge$phy_code_chung<- "Low"
update_merge$phy_code_chung[(update_merge$f_237 >= 3 & update_merge$f_239  >= 20) | (update_merge$f_243 >= 5 & update_merge$f_245  >= 30) | ((update_merge$f_237 + update_merge$f_243 +update_merge$f_249) >= 5 &  update_merge$totmets >= 600)]<- "Moderate"
update_merge$phy_code_chung[(update_merge$f_237 >= 3 & update_merge$totmets >= 1500) | ((update_merge$f_237 + update_merge$f_243 +update_merge$f_249) >= 7 &  update_merge$totmets >= 3000)]<- "High"
update_merge$phy_code_chung<-factor(update_merge$phy_code_chung, c("Low", "Moderate", "High"))
table(update_merge$phy_code, update_merge$phy_code_chung)
table(update_merge$phy_code_chung)
table(update_merge$phy_code)

# anxiety: sum(f_163:f_169), based on reference code
table(update_merge$f_163)
update_merge$f_163 <- tolower(update_merge$f_163)
update_merge$HADS_1[update_merge$f_163 == "most of the time"]<-3
update_merge$HADS_1[update_merge$f_163 == "a lot of the time"]<-2
update_merge$HADS_1[update_merge$f_163 == "occasionally"]<-1
update_merge$HADS_1[update_merge$f_163 == "not at all"]<-0

table(update_merge$f_164)
update_merge$f_164 <- tolower(update_merge$f_164)
update_merge$HADS_2[update_merge$f_164 == "very definitely and quite badly"]<-3
update_merge$HADS_2[update_merge$f_164 == "yes but not too badly"]<-2
update_merge$HADS_2[update_merge$f_164 == "a little but it does not worry me"]<-1
update_merge$HADS_2[update_merge$f_164 == "not at all"]<-0
table(update_merge$HADS_2)

table(update_merge$f_165)
update_merge$f_165 <- tolower(update_merge$f_165)
update_merge$HADS_3[update_merge$f_165 == "a great deal of the time"]<-3
update_merge$HADS_3[update_merge$f_165 == "a lot of the time"]<-2
update_merge$HADS_3[update_merge$f_165 == "not too often"]<-1
update_merge$HADS_3[update_merge$f_165 == "very little"]<-0
table(update_merge$HADS_3)

table(update_merge$f_166)
update_merge$f_166 <- tolower(update_merge$f_166)
update_merge$HADS_4[update_merge$f_166 == "not at all"]<-3
update_merge$HADS_4[update_merge$f_166 == "not often"]<-2
update_merge$HADS_4[update_merge$f_166 == "usually"]<-1
update_merge$HADS_4[update_merge$f_166 == "definitely"]<-0
table(update_merge$HADS_4)

table(update_merge$f_167)
update_merge$f_167 <- tolower(update_merge$f_167)
update_merge$HADS_5[update_merge$f_167 == "very often"]<-3
update_merge$HADS_5[update_merge$f_167 == "quite often"]<-2
update_merge$HADS_5[update_merge$f_167 == "occasionally"]<-1
update_merge$HADS_5[update_merge$f_167 == "not at all"]<-0
table(update_merge$HADS_5)

table(update_merge$f_168)
update_merge$f_168 <- tolower(update_merge$f_168)
update_merge$HADS_6[update_merge$f_168 == "very much indeed"]<-3
update_merge$HADS_6[update_merge$f_168 == "quite a lot"]<-2
update_merge$HADS_6[update_merge$f_168 == "not very much"]<-1
update_merge$HADS_6[update_merge$f_168 == "not at all"]<-0
table(update_merge$HADS_6)

table(update_merge$f_169)
update_merge$f_169 <- tolower(update_merge$f_169)
update_merge$HADS_7[update_merge$f_169 == "very often indeed"]<-3
update_merge$HADS_7[update_merge$f_169 == "quite often"]<-2
update_merge$HADS_7[update_merge$f_169 == "not very often"]<-1
update_merge$HADS_7[update_merge$f_169 == "not at all"]<-0
table(update_merge$HADS_7)

update_merge$HADS_score <- update_merge$HADS_1 + update_merge$HADS_2 + update_merge$HADS_3 + update_merge$HADS_4 + update_merge$HADS_5 + update_merge$HADS_6 + update_merge$HADS_7

update_merge$anx_code[update_merge$HADS_score<8] <- "Normal"
update_merge$anx_code[update_merge$HADS_score>=8 & update_merge$HADS_score<11] <- "Borderline Abnormal"
update_merge$anx_code[update_merge$HADS_score>=11] <- "Abnormal"
table(update_merge$anx_code)

# fruit and veg
update_merge$veg_portions_per_day <- update_merge$f_276 / 3
update_merge$veg_portions_per_week <- update_merge$f_275 * update_merge$veg_portions_per_day
update_merge$fruit_portions_per_week <- update_merge$f_277 * update_merge$f_278
update_merge$total_fruit_veg_portions_per_week <- update_merge$veg_portions_per_week + update_merge$fruit_portions_per_week
update_merge$avg_daily_portions <- update_merge$total_fruit_veg_portions_per_week / 7

write.table(meth_1114, file = "new_meth_1114.txt", sep = "\t", quote = FALSE, row.names = FALSE)
