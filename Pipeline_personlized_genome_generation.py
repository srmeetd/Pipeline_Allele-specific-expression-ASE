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
import cgat.BamTools as bamtools
#import CGAT.BamTools.bamtools as BamTools
import cgatpipelines.tasks.geneset as PipelineGeneset
#import cgatpipelines.pipeline_geneset as PipelineGeneset
import cgatpipelines.tasks.mapping as PipelineMapping
import cgatpipelines.tasks.mappingqc as PipelineMappingQC
import cgatpipelines.PipelineGeneset as PipelineGeneset
import cgatpipelines.tasks.rnaseq as PipelineRnaseq
import cgatpipelines.tasks.tracks as PipelineTracks
from cgatpipelines.report import run_report
import cgat.expression as Expression


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
PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    "genesets",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

PipelineGeneset.PARAMS = PARAMS
PipelineMappingQC.PARAMS = PARAMS

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except NameError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise.

def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])

    if not os.path.exists(PARAMS["annotations_database"]):
        raise ValueError(
            "can't find database '%s'" %
            PARAMS["annotations_database"])

    statement = '''ATTACH DATABASE '%s' as annotations''' % \
                (PARAMS["annotations_database"])

    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()
    return dbh
################################
#add phasing variant

############################################
# This function should be run if we found split bam files has bad cigar
##########################################


@follows(filtering,mkdir("patient_genomes.dir"))
@transform(filtering,
           regex("FilterRNA_editingsites.dir/(.+).vcf"),
           add_inputs(PARAMS["mapfasta"]),
           r"patient_genomes.dir/\1.fasta")
def generate_patient_genome(infiles, outfiles):
    '''This task uses generate_diploid_genome to make a
    diploid genome specific to each patient'''

    vcf, genome = infiles
    #sample = os.path.basename(vcf).split(".",2)[0]
    sample = os.path.basename(vcf).split(".vcf",1)[0]

    #genome2 = os.path.basename(vcf).split(".",2)[0] + ".genome2.fasta" 
    outfile = P.snip(os.path.basename(vcf), ".vcf")

    if not os.path.isfile("patient_genomes.dir/" + sample + ".genome1.fasta"):
        path_to_script = os.path.dirname(__file__)
        statement = ''' python %(path_to_script)s/generate_split_genomes.py
                          --stdin=%(vcf)s
                          --genome-fasta=%(genome)s
                          --log=%(outfile)s.log
                          --stdout=%(sample)s
                           '''
        job_memory = "8G"
        P.run(statement)

    if not os.path.isfile("patient_genomes.dir/" + sample + ".genome2.fasta"):
        path_to_script = os.path.dirname(__file__)
        statement = ''' python %(path_to_script)s/generate_split_genomes.py
                          --stdin=%(vcf)s
                          --genome-fasta=%(genome)s
                          --log=%(outfile)s.log
                          --stdout=%(sample)s
                           '''
        job_memory = "8G"
        P.run(statement)


@follows(generate_patient_genome)
@transform("patient_genomes.dir/*.fasta",
           formatter(),
           "patient_genomes.dir/{basename[0]}.dir")
def patient_genome_index(infiles, outfile):
    '''This task uses STAR to generate indeces for
       every diploid genome specific to each patient'''

    job_threads = 4
    job_memory = "16G"
    tempdir = P.get_temp_filename(dir=".")
    infile = infiles
    gtf =  PARAMS["gtf"]
    gf = "<(zcat %(gtf)s )" % locals()


    if not os.path.exists(outfile):
#        shutil.rmtree(outdir)
        statement = '''mkdir %(outfile)s &&
                      STAR  --runMode genomeGenerate
                      --genomeDir %(outfile)s
                      --genomeFastaFiles %(infile)s
                      --sjdbGTFfile %(gtf)s
                      --genomeSAsparseD 2
                      --sjdbOverhang 100
                      --outTmpDir %(tempdir)s 2&> %(outfile)s/index.log
                       '''



#        statement = ''' mkdir %(outfile)s;
#                    STAR  --runMode genomeGenerate
#                     --genomeDir %(outfile)s
#                     --genomeFastaFiles %(infile)s
#                     --sjdbGTFfile <(zcat %(gtf)s)
#                     --genomeSAsparseD 2
#                     --sjdbOverhang 100 
#                     --outTmpDir %(tempdir)s
#                      '''
        P.run(statement)
    else:
        pass
# gf = zcat %(gtf)s > %(outfile)s/geneset.gtf  rm %(outfile)s/geneset.gtf
#############################

#########################################################################
#########################################################################
#########################################################################
# Read mapping
#########################################################################

SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    "*.remote",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r".*/(\S+).(fastq.1.gz|fastq.gz|fa.gz|sra|csfasta.gz|csfasta.F3.gz|export.txt.gz|remote)")


@follows(mkdir("nreads.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"nreads.dir/\1.nreads")
def countReads(infile, outfile):
    job_threads = 2
    job_memory = "12G"
    '''Count number of reads in input files.'''
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run(statement)

#### need to make change so that it should consider each new genome for each fasqt files.

@follows(patient_genome_index, countReads, mkdir("star.dir"))
@subdivide(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           [r"star.dir/\1.genome1.star.bam",
            r"star.dir/\1.genome2.star.bam"])
def mapReadsWithSTAR(infile,outfiles):

    job_threads = PARAMS["star_threads"]
    job_memory = PARAMS["star_memory"]
    
   # name = os.path.basename(infile).split(os.extsep, 1)[0]
    name = os.path.basename(infile).split(".fastq")[0]


    star_mapping_genome1 =  name + ".star.genome1"
    star_mapping_genome2 =  name + ".star.genome2"

    for genome in [star_mapping_genome1, star_mapping_genome2]:
        
        star_mapping_genome = genome
        outfile = os.path.join("star.dir", os.path.basename(genome) + ".star.bam")
        m = PipelineMapping.STAR(
            executable=P.substitute_parameters(**locals())["star_executable"],
            strip_sequence=PARAMS["strip_sequence"])

        statement = m.build((infile,), outfile)
        P.run(statement)

@follows(mapReadsWithSTAR)
@merge("star.dir/*.bam", "star_stats.tsv")
def buildSTARStats(infiles, outfile):
    '''Compile statistics from STAR run
    Concatenates log files from STAR runs and reformats them as a tab
    delimited table.
    Parameters
    ----------
    infile: list
        :term:`bam` files generated with STAR.
    outfile: str
        :term: `tsv` file containing statistics about STAR run
    '''

    data = collections.defaultdict(list)
    for infile in infiles:
        fn =  re.sub(r'.star.bam', '', infile) + "Log.final.out"
        if not os.path.exists(fn):
            raise ValueError("incomplete run: %s" % infile)

        for line in IOTools.open_file(fn):
            if "|" not in line:
                continue
            header, value = line.split("|")
            header = re.sub("%", "percent", header)
            data[header.strip()].append(value.strip())

    keys = list(data.keys())
    outf = IOTools.open_file(outfile, "w")
    outf.write("track\t%s\n" % "\t".join(keys))
    for x, infile in enumerate(infiles):
        track = P.snip(os.path.basename(infile), ".bam")
        outf.write("%s\t%s\n" %
                   (track, "\t".join([data[key][x] for key in keys])))
    outf.close()

@follows(buildSTARStats)
@transform(buildSTARStats, suffix(".tsv"), ".load")
def loadSTARStats(infile, outfile):
    '''
    Loads statistics about a star run from a tsv file to a database table -
    star_stats.
    Parameters
    ----------
    infile: term:`tsv` file containing a table of tophat statistics.
    outfile: .load file logging database loading
    '''
    P.load(infile, outfile)


@follows(countReads)
@merge(countReads, "reads_summary.load")
def loadReadCounts(infiles, outfile):
    ''' load the read counts
    individual read counts are merged and loaded into a table called
    reads_summary
    Parameters
    ----------
    infile : str
       Input filename in :term:`tsv` format
    outfile : str
       Output filename, the table name is derived from `outfile`
    '''

    outf = P.get_temp_file(".")
    outf.write("track\ttotal_reads\n")
    for infile in infiles:
        track = P.snip(infile, ".nreads")
        lines = IOTools.open_file(infile).readlines()
        nreads = int(lines[0][:-1].split("\t")[1])
        outf.write("%s\t%i\n" % (track, nreads))
    outf.close()

    P.load(outf.name, outfile)

    os.unlink(outf.name)

@follows(mapReadsWithSTAR)
@transform(mapReadsWithSTAR,
           regex(".bam"),
           ".bed.gz")
def buildBed(infile, outfile):
    ''' Generate :term:`bed` format file from :term:`bam` alignment file
    Parameters
    ----------
    infile : str
       Input filename in :term:`bam` format
    outfile : str
       Output filename in :term:`bed` format
    '''

    statement = '''
    cat %(infile)s
    | cgat bam2bed
          %(bed_options)s
          --log=%(outfile)s.log
          -
    | sort -k1,1 -k2,2n
    | bgzip
    > %(outfile)s;
    tabix -p bed %(outfile)s
    '''
    P.run(statement)


@follows(generate_patient_genome, patient_genome_index,countReads, mapReadsWithSTAR, buildSTARStats, loadSTARStats, loadReadCounts)
def full():
    pass



@follows(mkdir("MultiQC_report.dir"))
@originate("MultiQC_report.dir/multiqc_report.html")
def renderMultiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LC_ALL=en_GB.UTF-8 && "
        "export LANG=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")
    P.run(statement)


@follows(renderMultiqc)
def build_report():
    '''report dummy task to build reports'''
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
