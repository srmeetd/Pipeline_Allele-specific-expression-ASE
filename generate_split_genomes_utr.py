
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
    parser.add_option("-f", "--input-format", dest="format", type="choice",
                      choices = ["bed", "vcf"],
                      default = "bed",
                      help="Format of varients input file, default is output from"
                      "intersection quasar input and bed")
    parser.add_option("-o", "--one-sequence", dest="one", action="store_true",
                      default=False,
                      help="Output all mutations in a single sequence")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)


    last_seq = None
    last_trans = None
    index_fasta = IndexedFasta.IndexedFasta(options.fasta)
    contig_sizes = index_fasta.getContigSizes()
    varients = options.stdin

    for line in varients:
        if(line[0].lstrip().startswith('#')):
            continue
        
        # parse input line
        line=line.strip().split('\t')
        start = int(line[1])
        ref_allele = line[3]
        alt_allele = line[4]
        UTR_start = int(line[11])
        UTR_end = int(line[12])
        contig = line[10]

        
        if options.format == "bed":
            rsID = line[5]
        elif optins.format == "vcf":
            rsID = line[2]

        strand = line[15]
 
        # getFasta puts only sequence location as sequence name, not transcript_id
        sequence_id = "%s:%i-%i(%s)" % (contig,UTR_start,UTR_end, strand)
        
        if last_seq is None:
            last_seq = sequence_id
            last_transcript = line[13]
            utr_sequence = list(index_fasta.getSequence(sequence_id, start=0, end = contig_sizes[sequence_id]))

        if options.one and \
           last_seq is not None and \
           sequence_id != last_seq:
            # Moved on to new UTR

            output = last_transcript
            fasta1 = FastaIterator.FastaRecord(output, "".join(utr_sequence))
            options.stdout.write(str(fasta1) + "\n")

            last_seq = sequence_id
            last_transcript = line[13]
            utr_sequence = list(index_fasta.getSequence(sequence_id, start=0, end = contig_sizes[sequence_id]))

        elif not options.one:
             utr_sequence = list(index_fasta.getSequence(sequence_id, start=0, end = contig_sizes[sequence_id]))


        if strand == "+":
            mutant = int(start-UTR_start)

        elif strand == "-":
            # need rev comp of bases
            ref_allele = Genomics.reverse_complement(ref_allele)
            alt_allele = Genomics.reverse_complement(alt_allele)
            
            # deal with 0-based, half-closed on reverse strand
            mutant = UTR_end - (start + 1)
           
        
        else:
            raise ValueError("Unstranded record in input")
      
        assert utr_sequence[mutant] == ref_allele, \
            "Base %i of sequence %s should be %s but is actaully %s" % \
            (mutant, sequence_id, ref_allele, utr_sequence[mutant])
    
        
        utr_sequence[mutant] = alt_allele

        if not options.one:
            # output name is >UTR_identifier:base_change:rsID"
            # where UTR_identifier should include transcript ID
            output = "%s:%s:%s" % (line[13],
                                   ref_allele + str(mutant) + alt_allele,
                                   rsID)
            fasta1 = FastaIterator.FastaRecord(output, "".join(utr_sequence))
            options.stdout.write(str(fasta1) + "\n")


    if options.one:
        output = last_transcript
        fasta1 = FastaIterator.FastaRecord(output, "".join(utr_sequence))
        options.stdout.write(str(fasta1) + "\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

