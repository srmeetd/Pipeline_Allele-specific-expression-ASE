'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''




from ruffus import *
import glob

from cgat import FastaIterator
import os
import sys
import re
from cgat import IndexedFasta
import cgatcore.iotools as IOTools
import cgatcore.experiment as E

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-fasta", dest="fasta", type="string",
                      help="Location of indexed fasta sequence for genome")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)


    last_contig = None
    index_fasta = IndexedFasta.IndexedFasta(options.fasta)
    contig_sizes = index_fasta.getContigSizes()
    last_contig = None
    gen1 = None
    gen2 = None
     
#        outfile_gen1 = open("genome_files.dir/" + os.path.basename(filename) + ".genome1.fasta", "w")
 #       outfile_gen2 = open("genome_files.dir/" + os.path.basename(filename) + ".genome2.fasta", "w")
	
    #stdout1 = options.stdout + ".genome1.fasta"
    #stdout2 = options.stdout + ".genome2.fasta" 
    stdout1 = "patient_genomes.dir/"+os.path.basename(options.stdout.name) + ".genome1.fasta"
    stdout2 = "patient_genomes.dir/"+os.path.basename(options.stdout.name) + ".genome2.fasta"
    #stdout1 = os.path.basename(options.stdout) + ".genome1.fasta"
    #stdout2 = os.path.basename(options.stdout) + ".genome2.fasta"    
    stdout1 = open(stdout1, "w")
    stdout2 = open(stdout2, "w")
    
    for line in options.stdin:
        if(line[0].lstrip().startswith('#')):
            continue
        line=line.split('\t')
        if (line[9].lstrip().startswith('./.')): 
            continue
        genotype=line[9].split(':') 
        #print(genotype[4])
        reference = genotype[4].split(",")
           # print(reference)
        if line[0] != last_contig:
        # some code here to output current chr
            if gen1 is not None:
                fasta1 = FastaIterator.FastaRecord(last_contig,"".join(gen1), fold=60)
                fasta2 = FastaIterator.FastaRecord(last_contig,"".join(gen2), fold=60)
                #options.stdout1.write(str(fasta1) + "\n")
                #options.stdout2.write(str(fasta2) + "\n")
#                    outfile_gen2.write(str(fasta2) +"\n")

                stdout1.write(str(fasta1) + "\n")
                stdout2.write(str(fasta2) + "\n")
                del fasta1
                del fasta2
            last_contig = line[0]
            #print(last_contig)
            gen1 = list(index_fasta.getSequence(last_contig, start=0, end = contig_sizes[line[0]]))
            gen2 = list(index_fasta.getSequence(last_contig, start=0, end = contig_sizes[line[0]]))
        
        if (genotype[0] =="0/1" and reference[1][-1] == "1"):
            gen1[int(line[1])-1] = line[4]
        
        elif(genotype[0] =="0/1" and reference[1][-1] == "2"):
            gen2[int(line[1])-1] = line[4]
	
        elif(genotype[0] =="1/1"):
            gen1[int(line[1])-1] = line[4]
            gen2[int(line[1])-1] = line[4]

    fasta1 = FastaIterator.FastaRecord(last_contig, "".join(gen1), fold=60)
    fasta2 = FastaIterator.FastaRecord(last_contig, "".join(gen2), fold=60)
			
   # options.stdout1.write(str(fasta1) + "\n")
   # options.stdout2.write(str(fasta2) + "\n")


    stdout1.write(str(fasta1) + "\n")
    stdout2.write(str(fasta2) + "\n")
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

