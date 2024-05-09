# Load packages
library(tidyverse)
library(tidyr)
library(dplyr)
library(readxl)
library(magrittr)
library(nnet)

### Import Gulf of Mexico data
GOMdata <- read_xlsx("GOMdata.xlsx", col_names = TRUE) %>%
  select(-Region, -Percent_ID, -Taxonomy)

###### CRITICAL: all likelihood odds are based on thought process of...
####      this combo of ISOLATE + MEDIA + METHOD + TEMP is the most common in data
###           how likely are other genera to be recovered in comparison to the "most successful" approach?
##              ~indirect measurement of how general these isolates' recoveries are

# Convert to factors and relevel based on best ref
GOMdata$Isolate <- factor(GOMdata$Isolate)
GOMdata$Isolate <- relevel(GOMdata$Isolate, ref = "Bacillus")

GOMdata$Media <- factor(GOMdata$Media)
GOMdata$Media <- relevel(GOMdata$Media, ref = "M1low")

GOMdata$Method <- factor(GOMdata$Method)
GOMdata$Method <- relevel(GOMdata$Method, ref = "Stamp")

GOMdata$Temp <- factor(GOMdata$Temp)
GOMdata$Temp <- relevel(GOMdata$Temp, ref = "RT")

# Perform multinomial logistic regression
GOMmultinom_model <- nnet::multinom(Isolate ~ Media + Method + Temp,
                                    data = GOMdata)
# Summarize the model
GOMmodel_summary <-summary(GOMmultinom_model)

# Extract coefficients, standard errors
GOMcoefficients <- GOMmodel_summary$coefficients
GOMstanderr <- GOMmodel_summary$standard.errors

# Check the Z-score for the model (wald Z)
zGOM <- GOMcoefficients/GOMstanderr
zGOM

# 2-tailed z test
pGOM <- (1 - pnorm(abs(zGOM), 0, 1)) * 2
pGOM

# Define a function to convert p-values to asterisks
pval_to_asterisks <- function(pval) {
  if (is.na(pval)) {
    return("NA")  # For missing or NA values
  } else if (pval <= 0.001) {
    return("***")  # Highly significant
  } else if (pval <= 0.01) {
    return("**")   # Very significant
  } else if (pval <= 0.05) {
    return("*")    # Significant
  } else {
    return("NA")     # Not significant
  }
}

# Convert p-values to asterisks
pGOMasterisks <- sapply(pGOM, pval_to_asterisks)

# Creating a matrix of coefficients with p-values
GOMpsummar <- paste(
  round(GOMcoefficients, 4),
  " (", pGOMasterisks, ")",
  sep = ""
) %>%
  matrix(nrow = nrow(GOMcoefficients))
# Assigning row and column names from GOMcoefficients
dimnames(GOMpsummar) <- list(rownames(GOMcoefficients), colnames(GOMcoefficients))

# Export the table to a csv text file
write.csv(GOMpsummar, file = "pGOMsummar.csv", row.names = T)





##### Import Antarctica only data
ARCdata <- read_xlsx("ARCdata.xlsx", col_names = TRUE) %>%
  select(-PercentID, -Taxonomy)

# Convert to factors and relevel based on best ref
ARCdata$Isolate <- factor(ARCdata$Isolate)
ARCdata$Isolate <- relevel(ARCdata$Isolate, ref = "Sporosarcina")

ARCdata$Media <- factor(ARCdata$Media)
ARCdata$Media <- relevel(ARCdata$Media, ref = "AMM")

ARCdata$Method <- factor(ARCdata$Method)
ARCdata$Method <- relevel(ARCdata$Method, ref = "D")

ARCdata$Temp <- factor(ARCdata$Temp)
ARCdata$Temp <- relevel(ARCdata$Temp, ref = "RT")

ARCdata$Depth <- factor(ARCdata$Depth)
ARCdata$Depth <- relevel(ARCdata$Depth, ref = "60ft")


# Perform multinomial logistic regression
ARCmultinom_model <- nnet::multinom(Isolate ~ Depth + Media + Method + Temp, 
                                    data = ARCdata)

# Summarize the model
ARCmodel_summary <-summary(ARCmultinom_model)

# Extract coefficients, standard errors, and then combine
ARCcoefficients <- ARCmodel_summary$coefficients
ARCstanderr <- ARCmodel_summary$standard.errors

# Check the Z-score for the model (wald Z)
zARC <- ARCcoefficients/ARCstanderr
zARC

# 2-tailed z test
pARC <- (1 - pnorm(abs(zARC), 0, 1)) * 2
pARC

# Convert p-values to asterisks
pARCasterisks <- sapply(pARC, pval_to_asterisks)

# Creating a matrix of coefficients with p-values
ARCpsummar <- paste(
  round(ARCcoefficients, 4),
  " (", pARCasterisks, ")",
  sep = ""
) %>%
  matrix(nrow = nrow(ARCcoefficients))
# Assigning row and column names from ARCcoefficients
dimnames(ARCpsummar) <- list(rownames(ARCcoefficients), colnames(ARCcoefficients))

# Export the table to a csv text file
write.csv(ARCpsummar, file = "pARCsummar.csv", row.names = T)
