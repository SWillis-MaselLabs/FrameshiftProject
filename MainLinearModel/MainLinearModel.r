# Main Linear Model 
# This is the linear model that was used to calculate the relative effect sizes of the various hypotheses 
# The data this linear model was run on is included in the directory. To run it, specify the path to the 
# directory on your computer to read in the dataframe.

library(MASS)
library(nlme)
library(lme4)
require(lmerTest)
library(optimx)

# The file that is needed to run this is located in this directory as FileForMainModel.csv
df <- read.csv(file = '', header = T)
head(df)

# below is used to determine the optimal value of lambda for a box cox transformation
bc <- boxcox(df$ISD~1, lambda = seq(.1,.7,0.01))
bc$x[which.max(bc$y)]

# lambda has been rounded to 0.4 for our analyses
lambda = 0.4

# We define the box cox transformation (note: we only included the power and did not include the
# scalars for the transformation. This isn't important since they are only a scaling factor.)
bc.transform <- function(x,L){
  x^L
}

# The transformed ISD values are then saved as a new column in the dataframe
ISD.transform <- bc.transform(df$ISD,lambda)
df$ISD.transform<-ISD.transform

# A mixed linear model is then created using two random effects and two fixed effects. 
# Designation (artificially-frameshifted non-overlapping controls, ancestral genes, novel genes) & Frame (+1 vs. +2) are the fixed effects
# Species and Homology group are the random effects
data.two.random.effects <- lmer(df$ISD.transform ~ df$Designation + df$Frame + (1|df$Species)+(1|df$HomologyGroup))

# The output from the linear model is then found using Summary
summary(data.two.random.effects)
