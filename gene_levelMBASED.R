#! /usr/bin/Rscript

library (optparse)
library (GenomicRanges)
library (MBASED)
#############
summarizeASEResults_1s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency=assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE=assays(MBASEDOutput)$pValueASE[,1],
    pValueHeterogeneity=assays(MBASEDOutput)$pValueHeterogeneity[,1]
  )
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID)))
  return(
    list(
      geneOutput=geneOutputDF,
      locusOutput=lociOutputList
    )
  )
}
########################

	args <- commandArgs(trailingOnly = TRUE)
	test = read.table(args[1], header = FALSE,sep= "\t")
	test_sample1 = subset (test, select = c("V1","V6","V7","V8","V9","V10","V12","V13","V4"))
	#t = test_sample1[!duplicated(test_sample1),]
	#test_sample1 = t
	test_sample1$sum = rowSums(test_sample1[7:8])
	test_sample1 = as.data.frame(test_sample1[test_sample1$sum > 30,])
#	t = test_sample1[!duplicated(test_sample1),]
#       test_sample1 = t
        test_sample1$V4 = toupper (test_sample1$V4)
        test_sample1 =  test_sample1[!grepl("(1 OF MANY)", test_sample1$V4),]
        test_sample1 =  test_sample1[!grepl("UNKOWN", test_sample1$V4),]
        t = test_sample1[!duplicated(test_sample1),]
        test_sample1 = t
        test_sample1$V4 = as.character(test_sample1$V4)
	test_sample1$V13 = as.integer(test_sample1$V13)
	test_sample1$V4 = gsub("^NA$","unknown",test_sample1$V4)
	#test_sample1$sum = NULL
	test_sample1$gene_name = paste0(test_sample1$V4,":",test_sample1$V10)
	mySNVs_2s <- GRanges(
	seqnames=test_sample1$V1,
	ranges=IRanges(start=test_sample1$V6, width=1),
	aseID=test_sample1$V4,
	allele1=test_sample1$V8,
	allele2=test_sample1$V9
	)
	names(mySNVs_2s) = test_sample1$gene_name
	
	m = SummarizedExperiment(assays =  list(lociAllele1Counts = matrix(test_sample1$V12,
                              ncol = 1, dimnames = list(names(mySNVs_2s),"sample1"))
                              ,lociAllele2Counts = matrix(test_sample1$V13,ncol = 1 ,
                              dimnames = list(names(mySNVs_2s),"sample1"))),
                              rowRanges=mySNVs_2s)
							  
	ASEresults_2s <- runMBASED(
	ASESummarizedExperiment=m,
	isPhased=FALSE,
	numSim=10^6,
	BPPARAM = SerialParam()
	)

	results = summarizeASEResults_1s(ASEresults_2s)$geneOutput
	write.table(results, file = paste0(args[2],"_unphased.tsv"), sep="\t", row.names=TRUE, quote=FALSE)

	
	phased <- runMBASED(
        ASESummarizedExperiment=m,
        isPhased=TRUE,
        numSim=10^6,
        BPPARAM = SerialParam()
        )

       results = summarizeASEResults_1s(phased)$geneOutput
       write.table(results, file = paste0(args[2],"_phased.tsv"), sep="\t", row.names=TRUE, quote=FALSE)

