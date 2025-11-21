# BKMR model for individual and mixed effect of air pollution exposure

# BKMR model only allow 1 column foreach pollutant, which is lagged exposure, i.e., average exposure during a certain period

setwd("~/OneDrive - Imperial College London/PhD/Project work/AIRWAVE/1129_meth")
#install.packages("bkmr")
library(bkmr)
library(tidyr)

meth_1114 <- read.delim("new_meth_1114.txt")
colnames(meth_1114)

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
covariates <- c("age_when_screened", "gender", "smoking3","alc_code", "edu_code", "hhincome_code", "anx_code", "total_fruit_veg_portions_per_week", "body_mass_index")
covs <- meth_1114_trimmed[,covariates]
cell <- c("CD8.naive", "CD8pCD28nCD45RAn","Gran", "NK", "PlasmaBlast", "Mono", "CD4T")
cellp <- meth_1114_trimmed[,cell]


pollution_index <- data.frame("index" = colnames(pollution))
pollution_index <- separate(pollution_index, index, into = c("M", "Pollutant", "year"), sep = "[_.]")

fit_bkmr <- function(clock, lag, adj.cell = F){ # lag is a integer 1-10
  # prepare exposure matrix 
  X <- rep(1, nrow(pollution))
  for (pollutant in c("NO2", "O3", "PM10", "PM25")){
    tempdf <- pollution[c(1:lag),pollution_index$Pollutant==pollutant]
    lag_exposure <- rowMeans(tempdf)
    X <- cbind(X, lag_exposure)
  }
  X <- X[,-1]
  colnames(X) <- c("NO2", "O3", "PM10", "PM25")
  
  # covariates
  if (adj.cell == FALSE){
    Z <- covs
  } else {
    Z <- cbind(covs, cellp)
  }
 
  
  Y <- y[,match(clock, clocks)]
}

