##Ordinal scale data statistical analysis and spineplot representations
#February 23, 2023 - Murat Yaman

#Here, we are going to generate a dummy data for ordinal scale analyses. Accordingly, a grade range between 0 and 3 will be utilized. In this regard, you can consider 
# - 0: Samples that do not show any signs of effect, 
# - 1: Samples with minor signs, 
# - 2: Samples with mild levels of effect
# - 3: samples with severe levels of the phenotype of interest
# We will assign uniform distribution for our control samples, while assigning skewed distributions to capture differences versus the control group
rm(list=ls())
gc()

suppressPackageStartupMessages(source('OrdinalScaleStatistics_PowerAnalyses_Functions.R'))

mingrade=0; maxgrade=3
grades=mingrade:maxgrade

set.seed(1994)
controls=sample(grades, size=50,replace = T, prob = rep(x = 1,length(grades)))
exp1=sample(grades, size=30,replace = T, prob = c(rep(x = 1,floor(length(grades)/2)), rep(x = 5,length(grades)-floor(length(grades)/2))))
exp2=sample(grades, size=40,replace = T, prob = rep(x = 1,length(grades)))
exp3=sample(grades, size=35,replace = T, prob = c(rep(x = 5,floor(length(grades)/2)), rep(x = 1,length(grades)-floor(length(grades)/2))))

ln1=max(length(controls), length(exp1)); ln2=max(length(exp2), length(exp3)); ln=max(ln1, ln2)

#Sample sizes do not have to be even between the groups for the follow-up analyses
dat=data.frame(controls=c(controls, rep(NA, ln-length(controls))), 
               exp1=c(exp1, rep(NA, ln-length(exp1))), 
               exp2=c(exp2, rep(NA, ln-length(exp2))), 
               exp3=c(exp3, rep(NA, ln-length(exp3))))


#Test the function
m=ordinalScaleStats(dataset = dat)

#Test the function for correct reference group assignment
n=ordinalScaleStats(dataset = dat,refgroup = 'control')

#Instead of significance signs, we can also present the p-values
o=ordinalScaleStats(dataset = dat, padjval = T)

