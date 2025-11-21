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
