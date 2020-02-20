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
    message:
        "lolilolz, ca parse les fastq"  
    input:
        raw_fastq_R1 = expand(result_repository + fastq_repository_name + "/{sample}_R1.fastq.gz",sample=SAMPLE_LIST),
        raw_fastq_R2 = expand(result_repository + fastq_repository_name + "/{sample}_R2.fastq.gz",sample=SAMPLE_LIST)
    #output:
    run:
        #Check number of sample in the fastq folder
        SAMPLE_NUMBER=len(os.listdir(result_repository + fastq_repository_name))/2
        SAMPLEFILE_NUMBER=len(SAMPLE_LIST)
        if (SAMPLEFILE_NUMBER<SAMPLE_NUMBER):
            print("WARNING:More samples in fastq directory than samplefile. Please verify the samplefile.\n")
        if (SAMPLEFILE_NUMBER>SAMPLE_NUMBER):
            print("WARNING:More samples in samplefile than fastq directory.Snakemake will stop.Please verify the samplefile.\n")

        shell("echo 'LELZ'")



