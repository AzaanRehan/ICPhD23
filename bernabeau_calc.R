# bernabeu
library(tidyverse)
library(survival)


## Loading in model coefficients. Make sure these files are present in the current working directory or change the paths to the correct directory.
coefficients <- read.delim("data/bage_coefficients.tsv")

## Loading in CpG coefficients for episcore projection
cpgs <- read.delim("data/cpg_episcore_weights.tsv")

data <- read.csv("~/928_meth/trimmed_betas.csv", stringsAsFactors = F)
# data <- as.data.frame(t(data))
coef <- data[intersect(rownames(data), cpgs$CpG_Site),]

## Check if Beta or M Values
message("2.5 Checking of Beta or M-values are present") 
m_to_beta <- function(val) {
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef <- if((range(coef, na.rm = T) > 1)[[2]] == "TRUE") {
  message("Suspect that M Values are present. Converting to Beta Values")
  m_to_beta(coef)
} else {
  message("Suspect that Beta Values are present");
  coef
}

## Scale data if needed
message("2.6 Scaling data (if needed)") 
ids <- colnames(coef)
scaled <- apply(coef, 1, function(x) sd(x, na.rm = T)) 

coef <-  if(range(scaled)[1] == 1 & range(scaled)[2] == 1) { 
  coef
} else { 
  coef_scale <- apply(coef, 1, scale)
  coef_scale <- t(coef_scale)
  coef_scale <- as.data.frame(coef_scale)
  colnames(coef_scale) <- ids
  coef_scale
}

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site
message("2.7 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from training sample") 
coef <- if (nrow(coef)==length(unique(cpgs$CpG_Site))) { 
  message("All sites present")
  coef
} else if (nrow(coef)==0) { 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)), c("CpG_Site", "Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - adding to dataset with mean Beta Value from training sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)), ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if (length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {
    missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat = mat*missing_cpgs1$Mean_Beta_Value
  coef = rbind(coef,mat)
} 

message("2.8 Convert NA Values to mean for each probe") 
## Convert NAs to Mean Value for all individuals across each probe 
na_to_mean <-function(methyl) {
  methyl[is.na(methyl)] <- mean(methyl, na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))


###### (3) Calculate Episcores
#########################################################################################################
#########################################################################################################

message("3. Calculating Episcores") 
loop <- unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 

## Save file
message("3.1. Exporting Episcores")  
write.table(out, "~/928_meth/episcore_projections.tsv", sep = "\t", quote = FALSE)

###### (5) bAge prediction
#########################################################################################################
#########################################################################################################
out <- read.delim("~/928_meth/episcore_projections.tsv", stringsAsFactors = F)
grim <- read.csv("~/928_meth/GrimAge_Results.csv")
grim <- grim[, c("DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1")]

message("4.1. Scale GrimAge components (if needed)") 
scaled_grim <- apply(grim, 2, function(x) sd(x, na.rm = T)) 
ids <- colnames(grim)

grim <-  if(range(scaled_grim)[1] == 1 & range(scaled_grim)[2] == 1) { 
  grim
} else { 
  grim_scale <- scale(grim)
  grim_scale <- as.data.frame(grim_scale)
  grim_scale
}

features <- read.csv("~/928_meth/update_meth.csv", stringsAsFactors = F)
df <- as.data.frame(cbind(features, out))
df <- df[!duplicated(df$barcode), ]
scores <- cbind(df, grim)

##Filter to elements in predictor
coefficients <- read.delim("data/bage_coefficients.tsv")
scores <- scores[,coefficients$Variable]


## Calculate bAge
message("5.1. Calculating bAge") 
scores <- t(scores)
pred <- scores * coefficients[,"Coefficient"]
pred_pp <- colSums(pred)

## Scale to same scale as age in testing
message("5.2. Scaling bAge") 
scale_pred <- function(x, mean_pred, sd_pred, mean_test, sd_test) { 
  scaled <- mean_test + (x - mean_pred)*(sd_test/sd_pred)
  return(scaled)
}

# Scale to same Z scale
scale_Z <- function(x, mean_pred, sd_pred) { 
  scaled <- (x - mean_pred)/sd_pred
  return(scaled)
}

mean_pred <- mean(pred_pp)
age <- scores[1,]
mean_test <- mean(age) # Mean age in testing data
sd_pred <- sd(pred_pp)
sd_test <- sd(age) # SD age in testing data

pred_pp_Z <- scale_Z(pred_pp, mean_pred, sd_pred)
pred_pp_scaled <- scale_pred(pred_pp, mean_pred, sd_pred, mean_test, sd_test)

## Make df with everything

pred_df <- data.frame(pred_pp_Z, pred_pp_scaled,  age)
names(pred_df) <- c("bAge", "bAge_Years", "Age")

## Obtain bAgeAccel
message("5.3. Obtaining bAgeAccel") 
pred_df$bAgeAccel <- resid(lm(bAge ~ Age, data=pred_df, na.action=na.exclude))

## Export
message("5.4. Exporting predictions") 
write.table(pred_df, "~/928_meth/bage_predictions.tsv", quote = FALSE, sep = "\t", row.names = T)
