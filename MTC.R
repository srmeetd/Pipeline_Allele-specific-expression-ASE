setwd("../")
system(paste0("mkdir Quasar_stats_results"))

setwd("Quasar_results_in_bed_format")

library(plyr)


files  <- list.files(pattern = '.bed$')
tables <- lapply(files, read.table,sep = "\t", header = F,stringsAsFactors = FALSE)
combined.df <- do.call(rbind , tables)
#value = split(combined.df,list(combined.df$V1,combined.df$V2,combined.df$V3,combined.df$V4),drop=TRUE)

value = dlply(combined.df, c("combined.df$V1","combined.df$V2","combined.df$V3","combined.df$V4"))

chip = function(x)
{
  y = as.data.frame(x)
  names(y) = c("chr","strt","end","rsID","Ref","Alt","Ref_count","Alt_count","betas","betas.se","Pval")
  chisq <- -2*sum(log(y$Pval))
  dof  = 2* length(y$Pval)
  combined_p = pchisq(chisq, df=dof, lower.tail=FALSE)
  beta_bion = mean(y$betas)
  file = y[!duplicated(y[,c('chr','strt',"end","rsID")]),]
  file$combined_pvalue = combined_p
  file$avg_betas = beta_bion
  file$n = nrow(y)
  return(file)
  
}

tmp = lapply(value,FUN=chip)
df = do.call(rbind,tmp)

df = subset (df, select = c("chr","strt","end","rsID","Ref","Alt","betas","betas.se","Pval","combined_pvalue", "avg_betas","n"))
row.names(df) = NULL
p_value = df[df[,"Pval"] <=0.1,]
p_value = unname(p_value)
write.table (p_value,"../Quasar_stats_results/p_value1_avg_betas.results",quote=F,row.names=F, sep = "\t")


df$FDR =  p.adjust(df$combined_pvalue, method="BH")
signif = df[df[,"FDR"] <=0.1,]
bed_format = unname(signif)
write.table (bed_format,"../Quasar_stats_results/FDR_0.05.avg_betas.results",quote=F,row.names=F, sep = "\t")


