import pysam
import pysamstats
import re
import sys

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv


bam1 = pysam.AlignmentFile(sys.argv[1])
outfile = open(sys.argv[2],"w") 
mflines=open('MAF_score.dir/merged_variants_corrected_MAF.vcf','r').readlines()

freq_dict={}
for line in mflines:
	if(line[0].lstrip().startswith('#')):
		continue
	L=line.split('\t')
	mafchr=L[0]
	m= int(L[1])
	mafstart=str(m-1)
	r=re.compile("MAF=.*")
	freq_dict[mafchr+mafstart]=list(filter(r.match, L[7].split(';')))[0][4:]



def genome(vcf_file=sys.argv[4]):
	
	file_lines=open(vcf_file,"r").readlines()
	
	
	for line in file_lines:
		if(line[0].lstrip().startswith('#')):
			continue
		line=line.split('\t')
		chr = line[0]
		start_pos = int(line[1])
		start = int(start_pos-1)
		end = int(start+1)
		ref = line[3]
		alt = line[4]
		db_id = line[2]
		
		strangealt_count=0
		
		for rec in pysamstats.stat_variation(bam1,chrom = chr, start = start, end = end,truncate=True, pad=True, fafile = sys.argv[3]):
			#print(rec)
			bases=['T','C','G','A']
			bases.remove(ref)
			if alt in "TCGA":
				bases.remove(alt)
				others_count=rec[bases[0]] +rec[bases[1]]
			else:
				others_count=0
			
			if(rec['chrom']+str(start) in freq_dict.keys()):
				frequency=freq_dict[rec['chrom']+str(start)]
			else:
				frequency=0
				
			if alt in "TCGA":
				outfile.write(rec['chrom'] + "\t" + str(start) + "\t" + str(end) + "\t" +  ref + "\t" + alt + '\t' + db_id + '\t' + str(frequency) + "\t" + str(rec[ref])+ "\t" + str(rec[alt]) + "\t" + str(others_count) + "\n")
			else:
				strangealt_count=+1
		#print(strangealt_count)

genome()
outfile.close()
