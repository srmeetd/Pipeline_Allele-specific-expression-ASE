import pysam
import sys

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

bam1 = pysam.AlignmentFile(sys.argv[1], "rb")
bam2 = pysam.AlignmentFile(sys.argv[2], "rb")

bam_out= pysam.AlignmentFile(sys.argv[3], "wb", header = bam1.header)
for read1,read2 in zip (bam1.fetch (until_eof=True),bam2.fetch (until_eof=True)):
    assert read1.query_name == read2.query_name
    #    print (read1.query_name == read2.query_name)
    if read1.has_tag('AS'):
        r1_AS = int(read1.get_tag('AS'))
    else:
        r1_AS = int(0)
    if read2.has_tag('AS'):
        r2_AS = int(read2.get_tag('AS'))
    else:
        r2_AS = int(0)
    if r1_AS >  r2_AS:
        bam_out.write(read1)
    elif r1_AS < r2_AS:
        bam_out.write(read2)
    elif r1_AS == r2_AS:
        if read1.mapping_quality > read2.mapping_quality:
            bam_out.write(read1)
        elif read1.mapping_quality < read2.mapping_quality:
            bam_out.write(read2)
        else:
            pass
       


bam1.close()
bam2.close()
bam_out.close()






