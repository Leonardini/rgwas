
## setwd("~/INI/ReverseGWAS")
library(ggplot2)
library(gridExtra)

### example for DNF data file but would be the same for the CNF file 
UKBB_Immuno_DNF_feasible_discovery_annotated=read.table('UKBB_Immuno_DNF_feasible_discovery_annotated.txt',sep='\t',header=T)

datDNF=as.data.frame(cbind(as.character(UKBB_Immuno_DNF_feasible_discovery_annotated$SNP_Discovery),
                           UKBB_Immuno_DNF_feasible_discovery_annotated$K,UKBB_Immuno_DNF_feasible_discovery_annotated$L,
                           -(UKBB_Immuno_DNF_feasible_discovery_annotated$bestSingleStat_Discovery),
                           UKBB_Immuno_DNF_feasible_discovery_annotated$ratio_discovery,UKBB_Immuno_DNF_feasible_discovery_annotated$ratio_validation,
                           -(UKBB_Immuno_DNF_feasible_discovery_annotated$bestSingleStat_Validation)))

colnames(datDNF)[1] ='SNP_discovery'
colnames(datDNF)[2] ='K'
colnames(datDNF)[3] ='L'
colnames(datDNF)[4] ='Single_SNP_discovery'
colnames(datDNF)[5] ='Ratio_discovery'
colnames(datDNF)[6] ='Ratio_validation'
colnames(datDNF)[7] ='Single_SNP_validation'
datDNF[,4] = as.numeric(as.character(datDNF[,4]))
datDNF[,5] = as.numeric(as.character(datDNF[,5]))
datDNF[,6] = as.numeric(as.character(datDNF[,6]))
datDNF[,7] = as.numeric(as.character(datDNF[,7]))

KL=interaction(datDNF$K,datDNF$L)
datDNF=cbind(datDNF,KL)
datDNF=datDNF[!datDNF$KL == 1.1,]

# final plot facet by complexity K,L across all feasible combinations
ggplot(datDNF,aes(Ratio_discovery,Ratio_validation,size=Single_SNP_discovery, fill=Single_SNP_validation)) + 
  geom_point(shape=21,color='black',alpha=0.8) + scale_size_area(max_size = 15) + facet_wrap(vars(K,L), ncol=3, labeller = "label_both",scales='free')+ 
  geom_hline(yintercept = 0, colour = 'red') + theme(strip.background = element_blank(), strip.placement = "outside") 


#ggplot(datDNF,aes(Ratio_discovery,Ratio_validation,size=Single_SNP_discovery, fill=Single_SNP_validation)) + geom_point(shape=21,color='black',alpha=0.8) + scale_size_area(max_size = 20) + facet_wrap(vars(KL))+ geom_hline(yintercept = 0, colour = 'red')
#ggplot(datDNF,aes(Ratio_discovery,Ratio_validation,size=Single_SNP_discovery, fill=Single_SNP_validation)) + geom_point(shape=21,color='black',alpha=0.8) + scale_size_area(max_size = 20) + facet_grid(rows=vars(K),cols=vars(L)) + geom_hline(yintercept = 0, colour = 'red'))

#### plot across complexities per SNP
uSNP=unique(datDNF$SNP_discovery)
d=split(uSNP, ceiling(seq_along(uSNP)/10))
p_list=list()
for(i in 1:length(d)){
  datDNF_tmp=datDNF[datDNF$SNP_discovery %in% as.character(uSNP[d[[i]]]),]
  p=ggplot(datDNF_tmp,aes(K,Ratio_discovery,fill=L)) + 
    geom_point(shape=21,color='black', size=6) + facet_wrap(vars(SNP_discovery), scales='free')+ 
    theme(strip.background = element_blank(), strip.placement = "outside") + 
    guides(fill = guide_legend(override.aes = list(size = 5))) + scale_fill_brewer()
  p_list[[i]] = p
}

ml <- marrangeGrob(p_list, nrow=1, ncol=1)
ggsave("ReverseGWAS_DNF_perSNP.pdf", ml)










