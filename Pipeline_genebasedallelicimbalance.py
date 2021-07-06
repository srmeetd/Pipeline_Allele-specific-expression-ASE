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
import cgatcore.experiment as E
from cgatcore import pipeline as P
import re
import glob
import collections
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
import cgat.BamTools.bamtools as BamTools
import cgatpipelines.tasks.geneset as PipelineGeneset
import cgatpipelines.tasks.mapping as PipelineMapping
import cgatpipelines.tasks.mappingqc as PipelineMappingQC


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

###################################
@follows (mkdir ("Genome_cordinates.dir"))
@transform (PARAMS["gtf"],formatter(),r"Genome_cordinates.dir/genome_cordinate.txt")

def genomcord (infile,outfile):
    job_memory = "16G"
    job_threads = 4
    statement = '''gtf2bed < %(infile)s | cut -f1-4 > %(outfile)s '''
    P.run(statement)


@follows(genomcord,mkdir("MABSED_input.dir"))
@transform ("Input_Quasar.dir/*.bed", regex(r"Input_Quasar.dir/(.+).bed"),
            add_inputs(r"Genome_cordinates.dir/*.txt"),
            r"MBASED_input.dir/\1.txt")

def bedintersect(infiles,outfile):
    bed,genome_cordinates =infiles
    statement = '''bedtools intersect -a %(genome_cordinates)s -b %(bed)s -wa -wb > %(outfile)s '''
    P.run(statement)

#############################################
@follows(bedintersect,mkdir("MBASED_results.dir"))
@transform ("MBASED_input.dir/*",formatter("(.+).txt"), ["MBASED_results.dir/{basename[0]}_phased.tsv", "MBASED_results.dir/{basename[0]}_unphased.tsv"])
#@transform (bedintersect,regex("MBASED_input.dir/(.+).txt"),r"MBASED_results.dir/\1.results.txt")

def MBASED (infiles,outfiles):
    job_memory = "16G"
    job_threads = 4
    infile = infiles
    phased,unphased = outfiles
    output =  P.snip(unphased, "_unphased.tsv")
    statement = ''' Rscript /data/md1srd/Software/AI_pipeline/gene_levelMBASED.R %(infile)s  %(output)s '''
    P.run (statement,job_condaenv="R-rstudio")

#############################################
@follows(bedintersect,MBASED)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


