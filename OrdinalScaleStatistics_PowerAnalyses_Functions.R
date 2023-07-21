##Functions to do statistical analyses on ordinal scale data and apply power analyses to determine minimum numbers of samples needed to detect significant differences

##February 26, 2023-Murat Yaman
##Updated July 21, 2023-Murat Yaman
##ordinalScaleStats: analyses ordinal scale data, produces spineplots and statistical descriptions if needed (showSign argument is set TRUE by default to obtain these results). Currently it accepts data in a tabular format where grades of ordinal measures for each sample is given across columns for each respective group. Users can also define a control group of their choice. Otherwise, sample from the first column will serve as the reference group for the follow-up analyses. Consequently, along with .csv format statistical description data, users can obtain detailed spineplots indicating p-values (numerically or by significance signs - padjval argument can be set TRUE/FALSE interchangeably, for these purposes). In addition, ordinal data analyses utilize log-odd ratio based evaluations. Therefore, samples that are dominated by one big proportion of a grade from the ordinal scale can result in +/- infinite values, which can scrutinize the analyses. To avoid that, you can use dummy variable where 1 dummy sample is assigned for each grade to each group. This can become quite useful, especially when your sample size is relatively large. To do so, you can set the dummie argument to TRUE. Besides that, decreasingGrades argument is set TRUE by default, assuming that higher levels of the grading system indicates severe cases while the lower ends represent milder situations. Please, make sure that you are representing the grading system numerically. If you have a different naming system for the ordinal scale, you can consider assigning numeric values for them. 

ordinalScaleStats=function(dataset,refgroup=NULL,decreasingGrades=T,dummie=F, dummynumber=1, showSign=T,padjval=F,pcorrection='dunnettx',ordinalmethod='clm', sheetname=NULL, colorpalette="BuGn",spineplot_x='Exposures',spineplot_y='Fraction of phenotypes',width_plot=16, height_plot=12, res_plot=300){
  require(tidyverse)
  require(reshape2)
  require(emmeans)
  suppressWarnings({
    if (is.null(refgroup)) {refgroup=colnames(dataset)[1]
    }  else {refgroup=refgroup}
    
    while (sum(tolower(refgroup) %in% tolower(colnames(dataset)))==0) { 
      refgroup=as.character(readline(prompt='Name of the reference group is not found in your data. Please, type it again:  '))}
    refgroup=colnames(dataset)[which(tolower(colnames(dataset))%in%tolower(refgroup))]
    
    if (is.null(ordinalmethod)) {ordinalmethod='clm'
    }  else {ordinalmethod=ordinalmethod}
    
    while (sum(tolower(ordinalmethod) %in% tolower(c('clm','polr')))==0) { 
      ordinalmethod=as.character(readline(prompt='Only, clm and polr methods are supported for the ordinal regresion. Please, type it again:  '))}
    ordinalmethod=tolower(ordinalmethod)
    if (ordinalmethod=='polr') {dummie=F}
    
    dataset=melt(data = dataset,na.rm = T) ##
    dataset=table(dataset); dataset=melt(dataset)
    colnames(dataset)=c('Exposures','Scores','Counts')
    
    dataset$Scores=factor(dataset$Scores,levels = sort(unique(dataset$Scores),decreasing = decreasingGrades), ordered = T)
    dataset$Exposures=relevel(dataset$Exposures, refgroup)
    
    if (dummie) {dataset$Counts=dataset$Counts+dummynumber}
    dataset=dataset[dataset$Counts>0,]
    
    dataset_long <- dataset %>% group_by(Exposures, Scores, Counts) |> uncount(Counts) %>% as.data.frame()
    
    if(showSign){
      
      if(ordinalmethod=='polr'){require(MASS)
        model = polr(formula = Scores~Exposures, Hess=T, data = dataset_long)
      }else{require(ordinal)
        model = clm(formula = Scores~Exposures , data=dataset_long)}
      
      summary(model);
      
      statdetails=emmeans(model, trt.vs.ctrl ~ Exposures, adjust=pcorrection)
      statdetails=as.data.frame(broom::tidy(statdetails$contrasts))
      
      filenames=paste0(gsub('-','',Sys.Date()),ifelse(is.null(sheetname),'', paste0('_',sheetname)),'_SpineplotStastics_',ordinalmethod,ifelse(dummie, paste0('_dummies',as.character(dummynumber)),'_withoutDummies'),'_1.csv')
      while (file.exists(filenames)) {
        filenames=paste0(gsub('_[^_]*$', '', filenames),'_', as.numeric(gsub('\\.csv','',gsub("^.+_", "", filenames)))+1,'.csv')
      }
      
      write.csv(statdetails,filenames, row.names = F)
      
      rm(filenames)
      
      if (length(levels(dataset_long$Exposures)<=2)) {signs=statdetails[,c('contrast', 'p.value')]
      colnames(signs)=c('contrast', 'adj.p.value')
      
      } else{signs=statdetails[,c('contrast', 'adj.p.value')]}
      
      signs$Exposures=gsub(paste0(' - ', refgroup), '',signs$contrast)
      require(stringr)
      
      signs$Exposures[startsWith(signs$Exposures,'(')&endsWith(signs$Exposures,')')]=str_sub(signs$Exposures[startsWith(signs$Exposures,'(')&endsWith(signs$Exposures,')')],2,-2)
      
      tmp_sign=as.data.frame(matrix(nrow = 1, ncol = 3,data = c(refgroup,NA,refgroup))); colnames(tmp_sign)=colnames(signs)
      
      signs=rbind(tmp_sign, signs);signs$adj.p.value=as.numeric(signs$adj.p.value)
      if (!0%in%round(signs$adj.p.value,3)) {
        signs$adj.p.value=round(signs$adj.p.value,3)
      } 
      
      
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
      label_df = dataset %>% group_by(Exposures, tot, contrast, adj.p.value, signif, pval_label) %>% dplyr::summarise(n=n())
      }else{
        dataset=dataset %>%
          group_by(Exposures, Scores) %>%
          summarise(Counts = sum(Counts)) %>%
          group_by(Exposures) %>%
          mutate(Prop = Counts / sum(Counts)) %>%
          mutate(tot = sum(Counts))
        label_df = dataset %>% group_by(Exposures, tot) %>% dplyr::summarise(n=n())}
    
    
    label_df$score<-NA
    
    filenames=paste0(gsub('-','',Sys.Date()),ifelse(is.null(sheetname),'', paste0('_',sheetname)),ifelse(showSign,'','_Annotated'),'_Spineplot_',ordinalmethod,ifelse(dummie, paste0('_dummies',as.character(dummynumber)),'_withoutDummies'),'_1.jpeg')
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
