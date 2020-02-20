#!usr/bin/en python3
import os
import pandas as pd

configfile:"config.yaml"

#Get information from config file
result_repository=config['Result_Repository']
fastq_repository_name=config['Project_folder']
samplefile=config['Samplesheet_Location']

#Get sample name list from the samplesheet
table_samplefile=pd.read_csv(samplefile,sep=";",header=0,)  
SAMPLE_LIST=list(table_samplefile['SAMPLE'])



rule all:
    input:
        raw_fastq_R1 = expand(result_repository + fastq_repository_name + "/{sample}_R1.fastq.gz",sample=SAMPLE_LIST),
        raw_fastq_R2 = expand(result_repository + fastq_repository_name + "/{sample}_R2.fastq.gz",sample=SAMPLE_LIST)
    output:
    
    
    message:
        "lolilolz, ca parse les fastq"  

