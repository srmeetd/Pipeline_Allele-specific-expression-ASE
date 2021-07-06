"""
rements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import shutil
import CGATCore.Experiment as E
from CGATCore import Pipeline as P
import re
import glob
import collections
import CGAT.GTF as GTF
import CGATCore.IOTools as IOTools
import CGAT.BamTools.bamtools as BamTools
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC

# Pipeline configuration
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS

# Add parameters from the annotation pipeline, but
# only the interface
#PARAMS.update(P.peek_parameters(
#    PARAMS["annotations_dir"],
#    "genesets",
#    prefix="annotations_",
#    update_interface=True,
#    restrict_interface=True))

####merged file#####

@follows(mkdir("merged_vcf"))
@merge ("FilterRNA_editingsites.dir/*.gz","merged_vcf/merged.vcf")
def merge_vcf(infile,outfile):
    job_memory = "8G"
    statement = '''bcftools merge FilterRNA_editingsites.dir/*.vcf.gz > %(outfile)s '''
    P.run(statement)
    

####MAF calculation#####

@follows(merge,mkdir("FilterHWE.dir"))
@transform("merge.dir/*.vcf",regex(r"merge.dir/merged.vcf"),r"FilterHWE.dir/HWE_filter.recode.vcf")

def HWE (infile,outfile):
    output = P.snip(outfile, ".recode.vcf")
    statement= '''vcftools --vcf %(infile)s --hwe 0.0001 --recode --out %(output)s'''
    P.run(statement)


@follows(HWE,mkdir("MAF_score.dir"))
@transform("FilterHWE.dir/*",regex(r"FilterHWE.dir/HWE_filter.recode.vcf"),r"MAF_score.dir/merged_variants_corrected_MAF.vcf")
def MAF (infile,outfile):
    job_memory="6G"
    #infile="merged_vcf/merged_variants.vcf"
    #outfile="MAF_score.dir/merged_variants_corrected_MAF.vcf"
    statement = ''' bcftools +fill-tags %(infile)s -Ov -o %(outfile)s -- -t AN,AC,AC_Hemi,AC_Hom,AC_Het,AN,HWE,MAF,NS '''
    P.run(statement)

####filter bam files######

@follows(mkdir("sorted_bams.dir"))
@transform("star.dir/*.bam",regex(".+/(.+).bam"),r"sorted_bams.dir/\1.sorted.bam")

def sort(infiles,outfile):
    job_memory = "16G"
    infile=infiles
    statement = '''samtools sort -n %(infile)s > %(outfile)s'''
    P.run(statement)


@follows(mkdir("Filter_bam.dir"))
@transform(sort, regex(".+/(.+).star.genome1.star.sorted.bam"),r"Filter_bam.dir/\1.star.filtered.bam")
def filter_bam (infiles,outfile):
    job_memory = "16G"
    infile=infiles
    bam1="sorted_bams.dir/" + os.path.basename(infile).split(".",3)[0] + ".star.genome1.star.sorted.bam"
    bam2="sorted_bams.dir/" + os.path.basename(infile).split(".",3)[0] + ".star.genome2.star.sorted.bam"
   # bam1="sorted_bams.dir/" + os.path.basename(infile).split(".star")[0] + ".star.sorted.bam"
   # bam2="sorted_bams.dir/" + os.path.basename(infile).split(".star")[0] + ".star.sorted.bam"
    outfile1 = "Filter_bam.dir/" + os.path.basename(infile).split(".",3)[0] + ".star.filtered.bam"

    outfile="Filter_bam.dir/" + os.path.basename(infile).split(".",3)[0] + ".filtered.bam"
    statement = '''python /data/md1srd/Software/pipeline_allele_imbalance/pipeline_split_genome/modified_split_genome_pipeline/filter_bam.py %(bam1)s %(bam2)s %(outfile1)s '''
    P.run(statement)


@follows(mkdir("sorted_filter_bam.dir"))
@transform("Filter_bam.dir/*.filtered.bam",
           regex(".+/(.+).star.filtered.bam"),
           r"sorted_filter_bam.dir/\1.sorted.bam")
def sort_filter (infile, outfile):
    job_memory = "8G"
    
    output = "sorted_filter_bam.dir/" + os.path.basename(infile).split(".")[0] + ".sorted.bam"
    statement = '''samtools sort %(infile)s > %(output)s && 
                   samtools index %(output)s '''
    P.run(statement)


@follows(sort_filter,mkdir("Input_Quasar.dir"))
@transform(sort_filter, 
           regex(".+/(.+).sorted.bam"),
           add_inputs(r"Filter_vcf.dir/\1.star.vcf"),
           r"Input_Quasar.dir/\1.bed")
           

def input_quasar(infiles,outfile):
    job_memory = "16G"
    bam,vcf = infiles
    path_to_script = os.path.dirname(__file__)
    fafile="patient_genomes.dir/" + os.path.basename(bam).split(".")[0] + ".star.genome1.fasta"
    output = "Input_Quasar.dir/" + os.path.basename(bam).split(".")[0] + ".bed"
    statement = '''python %(path_to_script)s/base_count.py %(bam)s %(output)s %(fafile)s %(vcf)s'''
    P.run(statement)

@follows(merge_vcf, MAF, sort, filter_bam, sort_filter,input_quasar)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


