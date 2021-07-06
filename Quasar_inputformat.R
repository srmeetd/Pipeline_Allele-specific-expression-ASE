####to convert frequency column into numeric

library (optparse)

#system(paste0("mkdir Input_to_quasar_format_corrected.dir"))
#setwd("Input_Quasar.dir/")
#bed = list.files(pattern  = ".bed$")

process = function(x,output_root)
{
  file = read.table (x,header=F,sep = "\t")
  file$V7 = as.numeric(as.character(file$V7))
  final = na.omit(file)
  final = unname(final)
  name = gsub(".bed","",x)
  #path = "../Input_to_quasar_format_corrected.dir/"
  write.table(final, paste0(output_root, ".bed"), sep="\t", row.names=FALSE, quote=FALSE)
  #write.table (final,paste0(path,name),sep = "\t",quote= F, row.names=F)
  system(paste0("gzip -c ", output_root,".bed", ">",  output_root,".gz" ))
}
#lapply (bed,process)

options_list <- list(
  make_option(c("-i", "--input-file"), action="store", dest="input_file"),
  make_option(c("-o", "--ouput-root"), action="store", dest="output_root"))

opt <- parse_args(OptionParser(option_list=options_list))

process(opt$input_file, opt$output_root)
