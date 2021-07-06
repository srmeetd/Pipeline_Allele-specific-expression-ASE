#setwd("../")
#system(paste0("mkdir Quasar_results_in_bed_format"))

#setwd("QUSAR_results.dir/")

library (optparse)
#file = list.files(pattern = ".results$")

bed = function (x,output_root)
{
  name = gsub (".results","",x)
  print(name)
  data = read.table (x,sep = "\t", header=T)
  bed = subset (data, select = c("chr","start","End","Rsid","Ref","alt","Ref_cout","Alt_count","betas","betas.se","pval2.het.ind."))
  bed= unname(bed)
  #path = "../Quasar_results_in_bed_format/"
  #print (paste("completed ",name))
  #write.table (bed, paste0(path,name,".bed"),quote = F, row.names =F,sep = "\t")
  write.table(bed, paste0(output_root, ".bed"), sep="\t", row.names=FALSE, quote=FALSE)
}

#lapply(file,FUN=bed)

options_list <- list(
  make_option(c("-i", "--input-file"), action="store", dest="input_file"),
  make_option(c("-o", "--ouput-root"), action="store", dest="output_root"))

opt <- parse_args(OptionParser(option_list=options_list))

bed(opt$input_file, opt$output_root)
