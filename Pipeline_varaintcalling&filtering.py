"""===========================
Pipeline Germline Variant Calling
===========================

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See ref`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_rnaseqmismatches.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
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
import CGATCore.Experiment as E
from CGATCore import Pipeline as P
import re

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS["projectsrc"] = os.path.dirname(__file__)
#for key, value in PARAMS.iteritems():
#    print "%s:\t%s" % (key,value)


# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.

PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


# ---------------------------------------------------
# Specific pipeline tasks

@follows(mkdir("readgroups.dir"))
@transform("input_files.dir/*.bam",formatter(),r"readgroups.dir/{basename[0]}.readgroups.bam")
def add_read_groups(infile, outfile):
    platform = PARAMS["platform"]
    groupsample = P.snip(os.path.basename(infile), ".bam")
    statement = '''java -Xmx8G -jar /shared/sudlab1/General/apps/bio/picard-tools-1.135/picard.jar
                   AddOrReplaceReadGroups
                   I=%(infile)s
                   O=%(outfile)s
                   RGLB=lib1
                   RGPL=%(platform)s
                   RGPU=unit1
                   RGSM=%(groupsample)s'''

    job_memory = "16G"
    P.run(statement)



@follows(mkdir("deduped.dir"))
@transform(add_read_groups,
           regex(r"readgroups.dir/(.+).readgroups.bam"),
           r"deduped.dir/\1.bam")
def dedup_bams(infile, outfile):
    '''Use MarkDuplicates to mark dupliceate reads'''
    job_memory = "16G"

    tempfile=P.snip(outfile, ".bam") + ".temp.bam"   
    metrics=P.snip(outfile, ".bam") + ".metrics.tsv"
    temporary = PARAMS["tmpdir"]
    statement = '''MarkDuplicates I=%(infile)s
                                  O=%(tempfile)s
                                  M=%(metrics)s
                                  TMP_DIR=%(temporary)s > %(outfile)s.log;

                                samtools view 
                                -F 1024
                                -b
                                %(tempfile)s
                                > %(outfile)s;
                  
                                rm -r %(tempfile)s;

                                samtools index %(outfile)s'''
    P.run(statement)


@follows(mkdir("split.dir"))
@transform(dedup_bams,regex(r"deduped.dir/(.+).bam"),r"split.dir/\1.split.bam")
def splitbams(infile,outfile):
    '''use GATK splitNcigar to split reads into exon segements'''
    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar /shared/sudlab1/General/git_repositories/GATK_file/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar 
                   -T SplitNCigarReads 
                   -R %(fastamap)s
                   -I %(infile)s 
                   -o %(outfile)s 
                   -rf ReassignOneMappingQuality 
                   -RMQF 255 
                   -RMQT 60 
                   -U ALLOW_N_CIGAR_READS ''' 
                                      

   
    job_memory = "32G"
    P.run(statement)

#@follows(mkdir("BaseRecalibration.dir"))
#@transform(splitbams,regex(r"split.dir/(.+).split.bam"),r"BaseRecalibration.dir/\1.recal.bam")
#def baserecal(infile,outfile):
#    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
#    fastamap = PARAMS["mapfasta"]
#    drctry= PARAMS["tmpdir"]
#    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar ~/Downloads/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar 
#                   -T BaseRecalibrator
#                   -R %(fastamap)s
#                   -I %(infile)s 
#                   -o %(outfile)s
#                   '''
#  
#    job_memory = "12G"
#    P.run(statement)
             

#fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"

@follows(mkdir("Variantcalls.dir"))
@transform(splitbams,regex(r"split.dir/(.+).split.bam"),r"Variantcalls.dir/\1.vcf.gz")
def variantcalling(infile,outfile):
    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    tempfile=P.snip(outfile,".gz")
    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar /shared/sudlab1/General/git_repositories/GATK_file/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar 
                   -T HaplotypeCaller
                   -R %(fastamap)s 
                   -I %(infile)s
                   -dontUseSoftClippedBases 
                   -stand_call_conf 10.0
                   --dbsnp /shared/sudlab1/General/projects/Sumeet/dbSNP/All_20180418_chr.vcf.gz
                   -o %(tempfile)s &&
                   bgzip -c %(tempfile)s > %(outfile)s &&
                   tabix -p vcf %(outfile)s
                   '''  

    job_memory = "16G"
    P.run(statement)         

###################################

@follows(variantcalling,mkdir("phased.dir"))
@transform("Variantcalls.dir/*.vcf.gz", regex(r"Variantcalls.dir/(.+).vcf.gz"),
            add_inputs(r"split.dir/\1.split.bam"),
            r"phased.dir/\1.vcf.gz")
def phasevariants(infiles, outfile):
    pass
    vcf, bam = infiles
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    tempfile=P.snip(outfile,".gz")
    statement = ''' java -Xmx10g -Djava.io.tmpdir=%(drctry)s -jar /shared/sudlab1/General/git_repositories/GATK_file/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
                          -T ReadBackedPhasing
                          -R %(fastamap)s
                          -I %(bam)s
                          --variant %(vcf)s
                          -L %(vcf)s
                          -o %(outfile)s
                          --phaseQualityThresh 20.0
                           '''

    job_memory = "16G"
    P.run(statement)

#######################################################
@follows(phasevariants,mkdir("Genotype_vcf.dir"))
@transform("phased.dir/*.vcf.gz",regex(r"phased.dir/(.+).vcf.gz"),
           add_inputs(r"split.dir/\1.split.bam"),r"Genotype_vcf.dir/\1.vcf.gz")


def gvcf(infiles,outfile):
    vcf,bam = infiles
    geno = P.snip(outfile,".gz")
    fastamap = PARAMS["mapfasta"]
    statement = '''freebayes
                   -f %(fastamap)s
                   -@ %(vcf)s %(bam)s > %(geno)s &&
                    bgzip -c %(geno)s > %(outfile)s &&
                    tabix -p vcf %(outfile)s '''
    job_memory = "16G"
    P.run(statement)


#-------------------------------------------------
@follows (mkdir("dbsnpid.dir"))
@transform(gvcf,regex(r"Genotype_vcf.dir/(.+).vcf.gz"),r"dbsnpid.dir/\1.vcf.gz")

def annotate(infile,outfile):    
    vcf = P.snip(outfile,".gz")
    statement = ''' bcftools 
                    annotate
                    -a /shared/sudlab1/General/projects/Sumeet/dbSNP/All_20180418_chr.vcf.gz
                    -c ID %(infile)s > %(vcf)s &&
                    bgzip -c %(vcf)s > %(outfile)s &&
                    tabix -p vcf %(outfile)s '''
    P.run(statement)
#---------------------------------------------------

#------------------------------------------------

@follows (annotate,mkdir("BiallicSNPs.dir"))
@transform(annotate,regex(r"dbsnpid.dir/(.+).vcf.gz"),r"BiallicSNPs.dir/\1.vcf.gz")

def readquality (infile,outfile):
    vcf  = P.snip (outfile,".gz")
    statement = ''' bcftools view
                    --max-alleles 2
                    --exclude-types indels
                    %(infile)s
                    -o%(outfile)s '''
    job_memory = "16G"
    P.run(statement)


#######################################

@follows (readquality,mkdir ("merge.dir"))
@merge ("BiallicSNPs.dir/*.gz","merge.dir/merged.vcf")

def merge (infiles,outfile):
    infile = infiles
    inputs = "BiallicSNPs.dir/*.gz"
    statement = '''bcftools merge %(inputs)s > %(outfile)s'''
    P.run (statement)

####################################

@follows(merge, mkdir("genome.dir"))
@transform(PARAMS["mapfasta"],formatter(),r"genome.dir/reference.fasta")
def reference_creation(infile, outfile):

    dirs = "genome.dir"
    name = PARAMS["name"]
    statement = '''ln -s %(infile)s %(dirs)s &&
                  mv genome.dir/%(name)s genome.dir/reference.fasta &&
                  sprint prepare genome.dir/reference.fasta bwa'''
    job_threads = 3
    job_memory = "16G"
    P.run(statement,job_condaenv="py2")


# --------------------------------------------------
@follows(reference_creation,mkdir("RNA_editingsites.dir"))
@transform("input.dir/*.fastq.gz", regex(r"input.dir/(.+).fastq.gz"),
           add_inputs("genome.dir/reference.fasta"),
           r"RNA_editingsites.dir/\1.dir")

def RES(infiles,outfile):
    samtools = PARAMS["samtool"]
    fastq,reference=infiles
    fas = os.path.basename(fastq).split(".gz")[0]
    path_to_script = os.path.dirname(__file__)
    statement = ''' gunzip -c %(fastq)s > input.dir/%(fas)s  &&
                   sprint main -1 input.dir/%(fas)s %(reference)s %(outfile)s bwa %(samtools)s &&
                   rm input_fastq.dir/%(fas)s &&
                   cd %(outfile)s && Rscript %(path_to_script)s/convert_bed_rnaediting.R'''

    job_threads = 3
    job_memory = "16G"
    P.run(statement,job_condaenv="py2")


###############################################

@follows(RES,mkdir("FilterRNA_editingsites.dir"))
@transform("RNA_editingsites.dir/*.RES.bed", regex(r"RNA_editingsites.dir/(.+).RES.bed"),
           add_inputs("BiallicSNPs.dir/*.star.vcf.gz"),
           r"FilterRNA_editingsites.dir/\1.vcf")

def filterRES(infiles,outfile):
    RES,variant=infiles
    vcf = P.snip(outfile,".gz")
  
    statement = ''' bedtools subtract -a %(variant)s -b %(RES)s -header -s > %(vcf)s && 
                    bgzip -c %(vcf)s > %(outfile)s &&
                    tabix -p vcf %(outfile)s '''

    job_threads = 3
    job_memory = "16G"
    P.run(statement,job_condaenv="py2")


#############################################


@follows(add_read_groups,dedup_bams,splitbams,variantcalling,phasevariants,gvcf,annotate,readquality,merge,reference_creation,RES,filterRES)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))




