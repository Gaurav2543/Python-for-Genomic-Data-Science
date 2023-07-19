# Program to compare sequence in a fasta file with the NCBI Database and fine th closest matches from the database with it  

import Bio
from Bio.Blast import NCBIWWW      #importing a large collection of NCBI built programs that we can call over the web
from Bio.Blast import NCBIXML    # This module formats things more nicely call over the web
fasta_string = open("/home/gaurav/Documents/Research/Python_for_Genomic_Data_Science/FASTA_files/myseq.fasta").read()
# print(fasta_string)
result_handle = NCBIWWW.qblast("blastn","nt",fasta_string) 

blast_record = NCBIXML.read(result_handle)
len(blast_record.alignments)       # this length is the number of closest matching sequences found, its max value can be 50
# blast_record.alignments contains the closest matching sequences
E_VALUE_THRESH = 0.01
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)    # hsp is high scoring pair (mean to measure distance between two sequences)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
            
# fasta_string.close    # close the file
