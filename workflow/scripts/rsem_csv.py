import csv
import pysam

samfile = pysam.AlignmentFile("WT_1.transcript.bam", "rb")

myfile = open('WT_1_rsem.csv', 'w')
wr = csv.writer(myfile)

header = ['qname', 'TF', 'pt']
wr.writerow(header)

for r in samfile:
    r_info = [r.query_name, r.reference_name, r.get_tag('ZF')]
    wr.writerow(r_info)
