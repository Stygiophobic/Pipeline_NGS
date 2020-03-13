#!/usr/bin/env python3
import os
import sys

#function for store fasta in a dictionary
def read_fasta(file_path):
    fasta_file=open(file_path,"r")
    fastas={}
    for line in fasta_file:
        if (line[0]=='>'):
            header=line[1:]
            fastas[header]=''
        else:
            fastas[header]+=line
    fasta_file.close()
    return fastas

fasta_sample=read_fasta("/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/BAM_FINAL/S20121_S52.fasta")

#work in progress
for sequence in fasta_sample:
    header=sequence.rstrip('\n')
    #print(header+"\n")
    DNA=fasta_sample[sequence].rstrip('\n')
    #print(DNA+"\n")

header="BVIC_S1"
print(header)
DNA=fasta_sample["BVIC_S1"+"\n"]
print(DNA)
newseq=""
