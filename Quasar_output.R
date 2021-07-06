
#setwd("../")
#system(paste0("mkdir QUSAR_results.dir"))

#setwd("Input_to_quasar_format_corrected.dir//")


library(QuASAR)
library (optparse)

#file = list.files (pattern =  ".gz$")

quasar = function(x,output_root)
{
  indiviual_file = read.table (x,sep = "\t")
  colnames(indiviual_file) =  c("chr","start","End","Ref","alt","Rsid","MAF","Ref_cout","Alt_count","other_count")
  name = gsub (".gz",".results",x)
  ase.dat <- UnionExtractFields(x, combine=TRUE)
  ase.dat.gt <- PrepForGenotyping(ase.dat = ase.dat, min.coverage=5)
  
  if (length(x) == 1) {
    ase.joint <- fitAseNull(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))
    
  } else {
    ase.joint <- fitAseNullMulti(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))
  }
  
  ourInferenceData <- aseInference(gts=ase.joint$gt, 
                                   eps.vect=ase.joint$eps, 
                                   priors=ase.dat.gt$gmat, 
                                   ref.mat=matrix(ase.dat.gt$ref, ncol = ncol(ase.dat$ref)),
                                   alt.mat=matrix(ase.dat.gt$alt, ncol = ncol(ase.dat$alt)),
                                   min.cov=10, 
                                   sample.names=c(x), 
                                   annos=ase.dat.gt$annotations)
  
  for (sample in length(ourInferenceData)) {
    m = merge (indiviual_file,ourInferenceData[[sample]]$dat, by.x = c("chr","start","Rsid"),
           by.y = c("annotations.chr","annotations.pos0","annotations.rsID"),all.y= T)
    #path = "QUSAR_results.dir/"
    write.table(m, file = paste0(output_root, ".results"), 
		sep="\t", row.names=FALSE, quote=FALSE)
    #write.table(m,file = paste0(path,name[sample]),
    #            quote = FALSE,
    #            sep="\t",
    #            row.names=FALSE)
  }
}

#tmp = lapply(file,FUN=quasar)

options_list <- list(
  make_option(c("-i", "--input-file"), action="store", dest="input_file"),
  make_option(c("-o", "--ouput-root"), action="store", dest="output_root"))

opt <- parse_args(OptionParser(option_list=options_list))

quasar(opt$input_file, opt$output_root)

 
