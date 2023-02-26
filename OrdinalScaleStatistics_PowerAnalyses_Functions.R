##Functions to do statistical analyses on ordinal scale data and apply power analyses to determine minimum numbers of samples needed to detect significant differences

##February 26, 2023-Murat Yaman

##ordinalScaleStats: analyses ordinal scale data, produces spineplots and statistical descriptions if needed (showSign argument is set TRUE by default to obtain these results). Currently it accepts data in a tabular format where grades of ordinal measures for each sample is given across columns for each respective group. Users can also define a control group of their choice. Otherwise, sample from the first column will serve as the reference group for the follow-up analyses. Consequently, along with .csv format statistical description data, users can obtain detailed spineplots indicating p-values (numerically or by significance signs - padjval argument can be set TRUE/FALSE interchangeably, for these purposes). In addition, ordinal data analyses utilize log-odd ratio based evaluations. Therefore, samples that are dominated by one big proportion of a grade from the ordinal scale can result in +/- infinite values, which can scrutinize the analyses. To avoid that, you can use dummy variable where 1 dummy sample is assigned for each grade to each group. This can become quite useful, especially when your sample size is relatively large. To do so, you can set the dummie argument to TRUE. Besides that, decreasingGrades argument is set TRUE by default, assuming that higher levels of the grading system indicates severe cases while the lower ends represent milder situations. Please, make sure that you are representing the grading system numerically. If you have a different naming system for the ordinal scale, you can consider assigning numeric values for them. 

ordinalScaleStats=function(dataset,refgroup=NULL,decreasingGrades=T,dummie=F, showSign=T,padjval=F,pcorrection='dunnettx',colorpalette="BuGn",spineplot_x='Exposures',spineplot_y='Fraction of phenotypes',width_plot=16, height_plot=12, res_plot=300){
  require(tidyverse)
  require(reshape2)
  require(emmeans)
  suppressWarnings({
    if (is.null(refgroup)) {refgroup=colnames(dataset)[1]
    }  else {refgroup=refgroup}
    
    while (sum(tolower(refgroup) %in% tolower(colnames(dataset)))==0) { 
      refgroup=as.character(readline(prompt='Name of the reference group is not found in your data. Please, type it again:  '))}
    refgroup=colnames(dataset)[which(tolower(colnames(dataset))%in%tolower(refgroup))]
    
    dataset=melt(data = dataset,na.rm = T) ##
    dataset=table(dataset); dataset=melt(dataset)
    colnames(dataset)=c('Exposures','Scores','Counts')
    
    dataset$Scores=factor(dataset$Scores,levels = sort(unique(dataset$Scores),decreasing = decreasingGrades), ordered = T)
    dataset$Exposures=relevel(dataset$Exposures, refgroup)
    
    if (dummie) {dataset$Counts=dataset$Counts+1}
    dataset=dataset[dataset$Counts>0,]
    
    if(showSign){
      require(ordinal)
      dataset_long= dataset %>% group_by(Exposures, Scores, Counts) |> uncount(Counts) %>% as.data.frame()
      model = clm(Scores~Exposures , data=dataset_long)
      summary(model);
      
      statdetails=emmeans(model, trt.vs.ctrl ~ Exposures, adjust=pcorrection)
      statdetails=as.data.frame(broom::tidy(statdetails$contrasts))
      
      filenames=paste0(gsub('-','',Sys.Date()),'_SpineplotStastics_1.csv')
      while (file.exists(filenames)) {
        filenames=paste0(gsub('_[^_]*$', '', filenames),'_', as.numeric(gsub('\\.csv','',gsub("^.+_", "", filenames)))+1,'.csv')
      }
      
      write.csv(statdetails,filenames, row.names = F)
      
      rm(filenames)
      
      signs=statdetails[,c('contrast', 'adj.p.value')]
      signs$Exposures=gsub(paste0(' - ', refgroup), '',signs$contrast)
      tmp_sign=as.data.frame(matrix(nrow = 1, ncol = 3,data = c(refgroup,NA,refgroup))); colnames(tmp_sign)=colnames(signs)
      signs=rbind(tmp_sign, signs);signs$adj.p.value=as.numeric(signs$adj.p.value)
      signs$adj.p.value=round(signs$adj.p.value,3)
      
      dataset=merge(dataset,signs,all.x=T, by='Exposures')
      dataset=dataset%>%mutate(signif=case_when(is.na(adj.p.value)~'',
                                                adj.p.value>0.05~'ns',
                                                adj.p.value<=0.05&adj.p.value>0.01~'*',
                                                adj.p.value<=0.01&adj.p.value>0.001~'**',
                                                adj.p.value<=0.001~'***'))
      
      dataset$pval_label='p-val=';dataset$pval_label[is.na(dataset$adj.p.value)]=''
      dataset$adj.p.value[is.na(dataset$adj.p.value)]=''
      
      dataset=dataset %>%
        group_by(Exposures, Scores, contrast, adj.p.value, signif, pval_label) %>%
        summarise(Counts = sum(Counts)) %>%
        group_by(Exposures, contrast, adj.p.value, signif, pval_label) %>%
        mutate(Prop = Counts / sum(Counts)) %>%
        mutate(tot = sum(Counts))
      label_df = dataset %>% group_by(Exposures, tot, contrast, adj.p.value, signif, pval_label) %>% dplyr::summarise(n=n())}else{
        dataset=dataset %>%
          group_by(Exposures, Scores) %>%
          summarise(Counts = sum(Counts)) %>%
          group_by(Exposures) %>%
          mutate(Prop = Counts / sum(Counts)) %>%
          mutate(tot = sum(Counts))
        label_df = dataset %>% group_by(Exposures, tot) %>% dplyr::summarise(n=n())}
    
    
    label_df$score<-NA
    
    filenames=paste0(gsub('-','',Sys.Date()),ifelse(showSign,'','_Annotated'),'_Spineplot_1.jpeg')
    while (file.exists(filenames)) {
      filenames=paste0(gsub('_[^_]*$', '', filenames),'_', as.numeric(gsub('\\.jpeg','',gsub("^.+_", "", filenames)))+1,'.jpeg')
    }
    
    p=ggplot(dataset, aes(x = Exposures, y = Prop, fill = Scores, label=Counts)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = scales::percent(Prop)),position = position_stack(vjust = 0.5))+
      labs(x=spineplot_x, y=spineplot_y)+
      theme(axis.title = element_text(size = 17, face = 'bold'),axis.text=element_text(size=15),axis.title.x =element_text(margin = margin(t=15)),axis.text.x=element_text(face = 'bold',angle=60, vjust = 1, hjust = 1),axis.title.y = element_text(margin =  margin(r= 15)), strip.text.y = element_text(size=12, face="bold.italic"),legend.text = element_text(size=12),legend.title = element_text(size=15))+ scale_fill_brewer(guide = guide_legend(reverse = TRUE),palette = colorpalette)+ scale_y_continuous(limits = c(0,1.1), breaks=seq(0,1,len=5))
    if(showSign & padjval){
      p=p+geom_text(data=label_df, aes(fill = score, x = Exposures, y=1,label=paste(pval_label,adj.p.value,'\n','n = ',tot)), position = position_dodge(0.9),vjust=-0.5, size = 15/.pt, inherit.aes = FALSE)}else if(showSign & !padjval) {p=p+geom_text(data=label_df, aes(fill = score, x = Exposures, y=1,label=paste(signif,'\n','n = ',tot)), position = position_dodge(0.9),vjust=-0.5, size = 15/.pt, inherit.aes = FALSE)}else{p=p+geom_text(data=label_df, aes(fill = score, x = Exposures, y=1,label=paste('n = ',tot)), position = position_dodge(0.9),vjust=-0.5, size = 15/.pt, inherit.aes = FALSE)}
    
    jpeg(filenames, units="in", width=width_plot, height=height_plot, res=res_plot)
    print(p)
    dev.off()})
  rm(filenames)
  return(dataset)
  
}

##Power analyses - simulation based - with a preliminary data
##Many thanks to Corey Powell @UMICH CSCAR for his valuable suggestions, support and guidance

##ordprobs() function allows calculating a vector of probabilities based on results from an ordinal #logistic regression. This is based on the parametrization of ordinal
#logistic regression used in the polr() function from the MASS package.
#See https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/
#intercepts - vector of intercepts from the ordinal logistic regression
#betas - vector of betas from the ordinal logistic regression
#Vals - vector of values of the variables in the ordinal logistic regression
ordprobs =  function(intercepts,betas,vals=NULL) {
  require(MASS)
  require(Hmisc)
  #The invariant proportional log-odds
  prop_odds = -1*betas %*% vals
  log_odds = intercepts + c(prop_odds)
  #Cumulative probabilities p(x <= j)
  cum_probs = c(sapply(log_odds,function(x){exp(x)/(1 + exp(x))}),1)
  #Calculates p(x = j)  = p(x <= j) - p(x <= j-1)
  probs = rep(0,length(cum_probs))
  probs[1]=cum_probs[1]
  for (i in 2:length(probs)) {
    probs[i] = cum_probs[i]-cum_probs[i-1]
  }
  return(probs)
}


##newPower2 was adopted and further modified from the suggestions of Greg Snow and gwern on #https://stats.stackexchange.com/questions/22406/power-analysis-for-ordinal-logistic-regression. Below function can estimate the power levels for detecting significant differences for a given sample size, probabilities of the ordinal scale data distribution for each sample, and one-hot encoding type values matrix as a reference map. To do so, it samples data from our distribution, makes models and estimates the powers. The function will be later utilized across a range of sample sizes which will be repeated multiple times to simulate the data.
newPower2 = function(n,vals_mat, probs) {
  #Due to random sampling we may run into NaN errors. To handle this, we can run consecutive attempts until we do not run into NaNs.
  require(MASS)
  nan_err=NULL
  while(is.null(nan_err)){
    lmodel=NULL
    attempts=1
    while( is.null(lmodel) && attempts <= 100 ) {
      attempts = attempts + 1
      ds = apply(probs, 2, sample, x=0:(nrow(probs)-1), size=100, replace=T)
      ds = apply(ds, 2, sample, size = round(n/ncol(ds)), replace=TRUE)
      ds = melt(ds); ds$Var1=NULL; colnames(ds) = c('OriginalNames', 'FDS')
      try({lmodel = polr(as.ordered(FDS)~OriginalNames , data=ds, Hess=TRUE)},silent = T)}
    if (sum(is.nan(summary(lmodel)$coefficients[1:length(unique(ds$OriginalNames))-1,3]))>0) {nan_err=NULL}  else {nan_err=summary(lmodel)$coefficients[1:length(unique(ds$OriginalNames))-1,3]
    #Collect the t-values
    ts=summary(lmodel)$coefficients[1:length(unique(ds$OriginalNames))-1,3]
    names(ts)=gsub('OriginalNames','',names(ts))
    #Let's see how much different our lmodel is in comparison to a null model. Accordingly, we can assign a variable named 'cumulative' that represents the difference between lmodel and null model. If the probability of them being similar is less than a certain threshold (0.05 or an adjusted p -value level - as arbitrarily chosen), the variable 'cumulative' will be updated to a certain t-value threshold - similar to the lmodel's. So that, in addition to sample-size estimations for each individual exposure comparisons, we can also consider them as a whole, and estimate total sample-size needed to reach a significance threshold in this type of experimental settings. 
    null = polr(as.ordered(FDS) ~ 1, data = ds, Hess = T)
    cumulative=0; names(cumulative)='Cumulative'
    # if (anova(lmodel,null)[2,7] < 0.05) {cumulative = cumulative+abs(qt(p = 0.05/2, df =lmodel$edf))} #edf: the (effective) number of degrees of freedom used by the polr() model
    # For a more strict p-value threshold, you can divide 0.05 by the numbers of tests (like p-value adjustment)
    if (anova(lmodel,null)[2,7] < 0.05/length(ts)) {cumulative = cumulative +length(unique(ds$OriginalNames))} 
    ts=append(ts, cumulative)
    }}
  
  return(ts)
}

# Sample sizes in our approach are estimated for the whole experiment itself. Meaning that, if our sample size is 200 and if we have 4 groups to study, we would need on average 50 samples per group. Users have the option to choose whether they would like to obtain figures for the whole experiment size (N) or approximately for each group on average (n). To do so, they can use the GroupSizeBasis argument T/F interchangeably. We can further play around with some of the aesthetics, too, and indicate whether we would like to get interactive plots (by plotly package) or static plots (ggplot2 package alone) - for static plots users can also define size and resolution aspects of the figures (width_heat and height_heat for width and heat, respectively, and res_heat for dpi resolution). In addition, sample size intervals to assess the power levels can be provided inside the nsample_min (minimum sample sizes), nsample_max (maximum sample sizes) and nsample_step (step sizes to advance from minimum  to maximum). Note: this values are assigned for total experiment size (N). For reproducibility, seed values to the set.seed() function can be assigned by the seedlevel argument
# dataset=dataset_long; powerlevel=80; refgroup=NULL; powergroups=NULL; GroupSizeBasis=T; interactivePlot=T; nsample_min=10; nsample_max=250; nsample_step=20; seedlevel=1994; width_heat=16; height_heat=10; res_heat=300
# refgroup='controllls'; powergroups=c('Experiment1', 'EXP2', 'exp3')

powerEstimator=function(dataset, powerlevel=80, refgroup=NULL, powergroups=NULL, GroupSizeBasis=T, interactivePlot=T, nsample_min=10, nsample_max=250, nsample_step=20, seedlevel=1994, width_heat=16, height_heat=10, res_heat=300){
  require(tidyverse); require(ordinal); require(MASS); require(Hmisc)
  require(boot)
  require(rms)
  require(mgcv) 
  suppressWarnings({
  if (is.null(refgroup)) {refgroup=as.vector(unique(dataset$Exposures))[1]
  }  else {refgroup=refgroup}
  
  while (sum(tolower(refgroup) %in% tolower(unique(dataset$Exposures)))==0) { 
    refgroup=as.character(readline(prompt='Name of the reference group is not found in your data. Please, type it again:  '))
    refgroup=as.vector(unique(dataset$Exposures)[which(tolower(unique(dataset$Exposures))%in%tolower(refgroup))])}
  
  if (is.null(powergroups)) {powergroups=as.vector(unique(dataset$Exposures)[!unique(dataset$Exposures)%in%refgroup])
  } else {powergroups=as.vector(unique(dataset$Exposures)[tolower(unique(dataset$Exposures))%in%tolower(powergroups)])
  powergroups=powergroups[!powergroups%in%refgroup]}
  
  resp='y'
  while (sum(tolower(powergroups) %in% as.vector(unique(dataset$Exposures)[!tolower(unique(dataset$Exposures))%in%tolower(refgroup)]))==0) { 
    powergroups=as.character(readline(prompt='Names of the test groups for power analyses are not found in your data. Please, provide a valid name:  '))
    resp=tolower(as.character(readline(prompt = 'Wanna add more groups: yes (Y) / no (N)   ')))
    while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more groups: yes (Y) / no (N)   ')))}
    while (resp!='n') {nextcomp=tolower(as.character(readline(prompt = 'Name of the next group: ')))
    while (sum(nextcomp %in% as.vector(unique(dataset$Exposures)[!tolower(unique(dataset$Exposures))%in%tolower(refgroup)]))==0) { nextcomp=as.character(readline(prompt = 'Please recheck the compound name, and provide it once again: '))}
    powergroups=unique(append(powergroups, nextcomp))
    powergroups=as.vector(unique(dataset$Exposures)[which(tolower(unique(dataset$Exposures))%in%tolower(powergroups))])
    resp=tolower(as.character(readline(prompt = 'Wanna add more groups: yes (Y) / no (N)   ')))
    while (!resp %in% c('y','n')) { resp=tolower(as.character(readline(prompt = 'Wanna add more groups: yes (Y) / no (N)   ')))}}
  }
  dataset=dataset[dataset$Exposures%in%c(refgroup,powergroups),]
  dataset$Exposures=as.character(dataset$Exposures)
  dataset$Exposures=factor(dataset$Exposures, levels = c(refgroup,powergroups))
  dataset=cbind(dataset,as.data.frame(model.matrix(~Exposures,data = dataset)))
  colnames(dataset)[colnames(dataset)=='(Intercept)']=paste0('Exposures',refgroup)  
  
  
  m1 = polr(Scores~Exposures , data=dataset, Hess=TRUE)
  
  exps1=as.data.frame(cbind(as.vector(levels(dataset$Exposures)),0))
  
  vals_mat1=as.data.frame(model.matrix(~0+V1,data = exps1))
  colnames(vals_mat1)=gsub('^V1','', colnames(vals_mat1))
  vals_mat1=vals_mat1[-1,]
  col_orders=levels(dataset$Exposures)[which(levels(dataset$Exposures)%in%colnames(vals_mat1))]
  vals_mat1=vals_mat1[,col_orders]
  
  probs1=apply(vals_mat1, 2, ordprobs,intercepts = m1$zeta,betas = coefficients(m1))
  
  n=seq(from=nsample_min, to=nsample_max, by=nsample_step)
  
  tmpreptable=as.data.frame(matrix(nrow = length(n), ncol = ncol(probs1)+1))
  tmpreptable$V1=n
  colnames(tmpreptable)=c('SampleSize', colnames(vals_mat1)[-1], 'Cumulative')
  tval=qt(p = 0.05/2, df = m1$edf)
  
  set.seed(seedlevel)
  ##Repeat the simulations 1000 times based on the probability levels that we obtained from our original data
  for (sampsize in 1:length(n)) {
    out=replicate(1000, newPower2(n = n[sampsize], probs= probs1, vals_mat = vals_mat1))
    #Obtain the power values for each time samples reach signficance (here we have used abs(t-value) score of 2 as a threshold)
    tmpreptable[sampsize, 2:ncol(tmpreptable)]=rowMeans( abs(out) >= abs(tval) ,na.rm = T)*100}
  
  tmpreptable_m=melt(tmpreptable, id.vars = 'SampleSize')
  colnames(tmpreptable_m)=c('SampleSize', 'Groups', 'Powers')
  
  #Let's create a table where we can store the data points corresponding to the assigned power levels for each group. We can later update the table with some additional columns. So that, it can become a reference table which can help us toimprove our plots' aesthetics
  int80=as.data.frame(matrix(nrow = 1, ncol = length(unique(tmpreptable_m$Groups))))
  colnames(int80)=unique(tmpreptable_m$Groups)
  
  #Extrapolate sample sizes from the curves
  for (powers in unique(tmpreptable_m$Groups)) { 
    dp.model <- mgcv::gam(Powers~ s(log(SampleSize)), data = tmpreptable_m[tmpreptable_m$Groups==powers,])
    int80[1,which(colnames(int80)==powers)] <- approx(x = dp.model$fitted.values, y = n, xout=powerlevel)$y
  }
  
  int80=melt(int80)
  int80$GroupSize=ceiling(int80$value/length(unique(tmpreptable_m$Groups)))
  colnames(int80)[colnames(int80)=='value']='TotalExperimentSize'
  int80$TotalExperimentSize=ceiling(int80$TotalExperimentSize)
  # int80=na.omit(int80)
  int80$TargetPower=paste0('%',powerlevel)
  
  #Modify our table with some additional columns which can help us to improve some of the aesthetics of the power plots that we will generate
  tmpreptable_m=merge(tmpreptable_m, int80, by.x='Groups',by.y='variable')
  tmpreptable_m$For80=paste0('To achieve ', tmpreptable_m$TargetPower,' power level: \nTotal experiment size (N) is: ', tmpreptable_m$TotalExperimentSize, '\nAverage size (n) needed for the ', tmpreptable_m$Groups,' group is: ', tmpreptable_m$GroupSize)
  tmpreptable_m$SampleSizeByGroupsAlone=tmpreptable_m$SampleSize/ncol(vals_mat1)
  
  int80$powers = paste0(int80$variable, ', Power level: ',int80$TargetPower,', Experiment size (N): ',ceiling(int80$TotalExperimentSize),', Group size (n): ', int80$GroupSize)
  
  
  ##Let's make our power analysis plots based on the size of total experiments or average group sizes in each experiment. For average group wise sample sizes, we included an additional column (GroupSize) by averaging the total sample size to the numbers of groups tested. Unfortunately, current version of the ggplot2 and plotly did not allow ifelse statements inside the aesthetics lists. Instead, we will go old school to ensure proper representations.
  require(ggplot2)
  require(grid)
  if(GroupSizeBasis){xint=unique(tmpreptable_m$GroupSize)
  }else{xint=unique(tmpreptable_m$TotalExperimentSize)}
  
  if(interactivePlot&GroupSizeBasis){p=ggplot(tmpreptable_m, aes(x=SampleSizeByGroupsAlone, y=Powers, col=Groups, text = paste0("\n",For80)))
  }else if(interactivePlot&!GroupSizeBasis){p=ggplot(tmpreptable_m, aes(x=SampleSize, y=Powers, col=Groups, text = paste0("\n",For80)))
  }else if(!interactivePlot&GroupSizeBasis){p=ggplot(tmpreptable_m, aes(x=SampleSizeByGroupsAlone, y=Powers, col=Groups))
  }else{p=ggplot(tmpreptable_m, aes(x=SampleSize, y=Powers, col=Groups))}
  
  p=p+geom_point()+ 
    geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)))+theme_bw()+theme(plot.title = element_text( size=rel(1.4), face="bold",hjust = 0.5), strip.text.x = element_text(size = rel(1.2), face = "bold.italic"), strip.background = element_rect(color="black", fill="gray66", size=1.5, linetype="solid"),axis.text.x=element_text(size=rel(1),face = 'bold',angle=0),axis.title.x = element_text(size=rel(1.2), face="bold"),axis.title.y = element_text(size=rel(1.2), face="bold"))+
    labs(x= ifelse(GroupSizeBasis,'Group sizes (n)','Total experiment size (N)'), y='Power (%)')+ggtitle('Simulations on the original data')+
    scale_color_brewer(palette="Dark2")+ scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))+
    scale_x_continuous(limits = c(0, ifelse(GroupSizeBasis,max(tmpreptable_m$SampleSizeByGroupsAlone)+10,nsample_max)), breaks = seq(0, ifelse(GroupSizeBasis,max(max(tmpreptable_m$SampleSizeByGroupsAlone)+10,100),max(max(tmpreptable_m$SampleSize)+10,nsample_max)), by = ifelse(GroupSizeBasis,10,nsample_step)))+
    geom_hline(yintercept=powerlevel, linetype="dashed", color = "gray", size=1, alpha=.9)+geom_vline(xintercept = xint, linetype="dashed", color = "gray", size=1, alpha=.9)
  
  filenames=paste0(gsub('-','',Sys.Date()),'_PowerAnalyses_PowerLevel',powerlevel,ifelse(interactivePlot,'_InteractivePlot_1.html','_StaticPlot_1.jpeg'))
  while (file.exists(filenames)) {if(!interactivePlot){
    filenames=paste0(gsub('_[^_]*$', '', filenames),'_', as.numeric(gsub('\\.jpeg','',gsub("^.+_", "", filenames)))+1,'.jpeg')
  }else{filenames=paste0(gsub('_[^_]*$', '', filenames),'_', as.numeric(gsub('\\.html','',gsub("^.+_", "", filenames)))+1,'.html')}
  }
  
  if (!interactivePlot) {
    cnts=0
    jpeg(filenames, units="in", width=width_heat, height=height_heat, res=res_heat) 
    print(p)
    grid.text(expression(underline('Sample sizes (% power levels)')),x = 0.05, y=(0.95), gp=gpar(fontsize=13, fontface='bold'),just = "left")
    for (i in 1:(nrow(int80))) {grid.text(int80$powers[i],x = 0.05, y=(0.95-i*0.022), gp=gpar(fontsize=12),just = "left")
      cnts=cnts+1}
    grid.polygon(x=c(0.04,0.45,0.45,0.04), y=c((0.95-cnts*0.03),(0.95-cnts*0.03),0.96,0.96), gp = gpar(fill="gray", alpha=0.2))
    dev.off()} else{ 
      require(htmlwidgets)
      require(plotly)
      saveWidget(ggplotly(p), file=filenames,selfcontained = T)}
  return(tmpreptable_m)
  })}
