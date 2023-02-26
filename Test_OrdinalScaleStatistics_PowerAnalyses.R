##Ordinal scale data statistical analysis and spineplot representations
#February 26, 2023 - Murat Yaman

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
dat=ordinalScaleStats(dataset = dat, padjval = T)

dataset_long= dat %>% group_by(Exposures, Scores, Counts) |> uncount(Counts) %>% dplyr::select(Exposures, Scores) %>% as.data.frame()

#Let's estimate the power levels and obtain representative figures. Note1: Since there is some computational power needed to generate simulations, you may need to wait for a while for the function to execute. Note2: you may be prompted to modify some of the arguments, in case there is some misspelling or missing data. For example, the command below will prompt you with correcting your control, and will take EXP2 and exp3 into account in making the analyses, since there is no group called Experiment1 in our samples, yet the function can handle case-sensitivity.
a=powerEstimator(dataset = dataset_long,refgroup = 'controllls', powergroups = c('Experiment1', 'EXP2','exp3'))

#No need to provide a powergroups argument if you want to test for all groups in your experiment. It will automatically assign the groups for you/ Please, just be careful with the reference argument. If not provided, the function will automatically select the first group from the list of groups in your data.
b=powerEstimator(dataset = dataset_long,refgroup = 'controls')

#You can also change the power level (default 80), range of sample sizes (default 10-250) to be tested and the set.seed (default 1994) arguments
c=powerEstimator(dataset = dataset_long,refgroup = 'controls', powerlevel = 90, seedlevel = 2, nsample_min = 30, nsample_max = 500, nsample_step = 50)

#If you prefer, you can have a look at the total experiment size (N) needed for the whole experiment, instead of group sizes (n). In addition, you can also generate static plots rather than interactive plots and html documents
d=powerEstimator(dataset = dataset_long,refgroup = 'controls', GroupSizeBasis = F, interactivePlot = F, seedlevel = 2)

#Note: you may wonder what the Cumulative is, It mainly represents the difference in distribution of our data (model) versus a null model. Hence, it can also show on average how much sample size would be needed to indicate a difference in our settings against a null hypothesis. Nonetheless, group wise levels are more definitive towards drawing implications based on individual contrasts between the respective groups against the controls. Therefore, the variable cumulative would give you some holistic view, while each individual group requiring specific sample size thresholds to represent significant differences.
