cfpr = read.csv('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Analysis of all data/final_combined_results_CFPR.csv',header=T)
pyr = read.csv('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Analysis of all data/final_combined_results_PYR.csv',header=T)
comb = read.csv('/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Analysis of all data/final_combined_results_Combined.csv',header=T)

merged = merge(cfpr,comb,by.x = 'ID',by.y='ID',all=T)

merged = merge(merged,pyr,by.x = 'ID',by.y='ID',all=T)


write.csv(merged,'/Users/vickyingham/Seafile/Wasim & Vicky/RNAseq/Analysis of all data/final_combined_results_fixed.csv',row.names=F)
