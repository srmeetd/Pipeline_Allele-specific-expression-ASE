
f = list.files (path = "RNA_editingsites.dir/", 
           pattern = "SPRINT_identified_all.res",
           recursive  = TRUE,full.names = TRUE)

convert_bed = function(x)
{
  a = read.table(x, header = FALSE, sep = "\t")
  names = gsub ("RNA_editingsites.dir//|/SPRINT_identified_all.res|.dir","",x)
  rna_editing_sites = subset (a,select = c("V1","V2","V3"))
  colnames(rna_editing_sites) = NULL
  write.table (rna_editing_sites,paste0("RNA_editing_sites/",names,"_RES.bed"), sep = "\t",
               quote = FALSE, row.names = FALSE)
}
lapply(f,convert_bed)
