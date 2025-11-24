# BKMR model for individual and mixed effect of air pollution exposure

# BKMR model only allow 1 column foreach pollutant, which is lagged exposure, i.e., average exposure during a certain period

setwd("~/OneDrive - Imperial College London/PhD/Project work/AIRWAVE/1129_meth")
#install.packages("bkmr")
library(bkmr)
library(tidyr)

meth_1114 <- read.delim("new_meth_1114.txt")
#colnames(meth_1114)

## this screen time is accurate to year-month-day h:m:s. here I just use the year and month and ignore day and time
whenscreen_year_month <- separate(
  meth_1114["first_clinic_date"],
  first_clinic_date,
  into = c("year", "month", "day", "hour", "min", "second"),
  sep = "[^[:alnum:]]+"
)
whenscreen_year_month <- as.data.frame(sapply(whenscreen_year_month ,as.numeric))
#head(whenscreen_year_month)

#summary(whenscreen_year_month$year<2008) ## 180 people had their blood sample collected in 2007, which is not covered by air pollution exposure estimates
meth_1114_trimmed <- meth_1114[whenscreen_year_month$year>=2008,]
pollution <- meth_1114_trimmed[,4:243]
whenscreen_year_month <- subset(whenscreen_year_month, year >= 2008)
## the starting month for air pollution estimates is Jan-2008
## here we need to know for each person, how many months of exposure data they have, so we can decide on lag intervals
avai_mon_exposure <- (whenscreen_year_month$year-2008)*12+whenscreen_year_month$month
hist(avai_mon_exposure, main = "Month of Prior Exposure", xlab = "month")


## to define the lag time period, we check the statistical power, and for a correlations r/OR >= 0.10 between the exposure index and continuous outcome (α = 0.05, power = 0.80):
## the minimum sample size is 782
freq <- table(avai_mon_exposure)
cum_freq <- cumsum(freq)
print(cum_freq[cum_freq<948-782]) ## showing the maxium lag intervel is 10M

clocks <- c("AgeAccelPheno", "AgeAccelGrim", "bAgeAccel", "epiTOC1", "DunedinPACE") 
y <- meth_1114_trimmed[,clocks]

## BKMR requires a full numeric covariate matrix
covariates <- c("age_when_screened", "gender", "smoking3","edu_code", "hhincome_code", "anx_code", "total_fruit_veg_portions_per_week", "body_mass_index", "drinker_category")
covs <- meth_1114_trimmed[,covariates]
covs$gender <- ifelse(covs$gender=="M", 1, 0)
covs$smoking3 <- as.numeric(factor(covs$smoking3, levels = c("never smoker", "former smoker", "current smoker")))
covs$edu_code <- as.numeric(factor(covs$edu_code))
covs$hhincome_code <- as.numeric(factor(covs$hhincome_code))
covs$drinker_category <- as.numeric(factor(covs$drinker_category))
covs$anx_code <- as.numeric(factor(covs$anx_code))
cell <- c("CD8.naive", "CD8pCD28nCD45RAn","Gran", "NK", "PlasmaBlast", "Mono", "CD4T")
cellp <- meth_1114_trimmed[,cell]


pollution_index <- data.frame("index" = colnames(pollution))
pollution_index <- separate(pollution_index, index, into = c("M", "Pollutant", "year"), sep = "[_.]")

fit_bkmr <- function(clock, lag = 1, adj.cell = "Y"){ # lag is a integer 1-10
  set.seed(1234)
  # prepare exposure matrix 
  temp_pollution <- pollution[avai_mon_exposure>=lag,]
  Z <- rep(1, nrow(temp_pollution))
  for (pollutant in c("NO2", "O3", "PM10", "PM25")){
    tempdf <- pollution[avai_mon_exposure>=lag,pollution_index$Pollutant==pollutant]
    temp_whenscreen_year_month <- whenscreen_year_month[avai_mon_exposure>=lag,]
    temp_whenscreen_year_month$month_index <- (temp_whenscreen_year_month$year-2008)*12+temp_whenscreen_year_month$month
    start <- temp_whenscreen_year_month$month_index-lag+1
    end <- temp_whenscreen_year_month$month_index
    lag_exposure <- integer()
    for (i in 1:nrow(tempdf)){
      start_i <- start[i]
      end_i <- end[i]
      lag_exposure[i] <- rowMeans(tempdf[i,start_i:end_i])
    }
    Z <- cbind(Z, lag_exposure)
  }
  Z <- Z[,-1]
  colnames(Z) <- c("NO2", "O3", "PM10", "PM25")
  
  # covariates
  if (adj.cell == "N"){
    X <- covs[avai_mon_exposure>=lag,]
  } else {
    X <- cbind(covs[avai_mon_exposure>=lag,], cellp[avai_mon_exposure>=lag,])
  }
  
  # age acc
  Y <- y[avai_mon_exposure>=lag,match(clock, clocks)]
  
  outfile <- paste(clock, paste0(lag, "M"), paste0("Cell", adj.cell), sep = "_")
  res <- kmbayes(y = Y,  X = X, Z = Z, iter = 10000, varsel = T) # MCMC iterations x 10000 (for parameter adjustment), select variables via PIP = T
  print("#########################################################################")
  print("#########################################################################")
  print("#########################################################################")
  print("#########################################################################")
  print("#########################################################################")
  print(outfile)
  print(summary(res))
  print("#########################################################################")
  print("#########################################################################")
  print("#########################################################################")
  print("#########################################################################")
  print("#########################################################################")
  
  PIP <- ExtractPIPs(res)
  write.table(PIP, paste(outfile, "PIP.txt", sep = "_"), quote = F, row.names = F, sep = "\t")
  
  SVSum <- SingVarRiskSummaries(res)
  write.table(SVSum, paste(outfile, "SVSum.txt", sep = "_"), quote = F, row.names = F, sep = "\t")
  
  AllSum <- OverallRiskSummaries(res)
  write.table(AllSum, paste(outfile, "AllSum.txt", sep = "_"), quote = F, row.names = F, sep = "\t")
 }

fit_bkmr("AgeAccelPheno", lag = 5)

X
?kmbayes


summary(res)
