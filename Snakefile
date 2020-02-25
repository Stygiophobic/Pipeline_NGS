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
        hg19_ref = "REF_HG19/hg19.fa" ,
        raw_fastq_R1 = expand(result_repository + fastq_repository_name + "/{sample}_R1.fastq.gz",sample=SAMPLE_LIST),
        raw_fastq_R2 = expand(result_repository + fastq_repository_name + "/{sample}_R2.fastq.gz",sample=SAMPLE_LIST) ,
        database = expand("REF_HG19/hg19.fa.{ext}", ext=["nhr", "nin", "nsq"]),
        bitmask = "REF_HG19/hg19.bitmask",
        srprism = expand("REF_HG19/hg19.srprism.{subfile}",subfile=['amp','idx','imp','map','mpm','rmp','ss','ssa','ssd'])



rule check_samplefile:
    message:
        "Checking samplefile, and decompress fastq."  
    input:
        raw_fastq_R1 = expand(result_repository + fastq_repository_name + "/{sample}_R1.fastq.gz",sample=SAMPLE_LIST),
        raw_fastq_R2 = expand(result_repository + fastq_repository_name + "/{sample}_R2.fastq.gz",sample=SAMPLE_LIST)
    output:
        fastq_R1 = result_repository + "FASTQ/{sample}_R1.fastq.gz",
        fastq_R2 = result_repository + "FASTQ/{sample}_R2.fastq.gz"    
    run:
        #Check number of sample in the fastq folder
        SAMPLE_NUMBER=len(os.listdir(result_repository + fastq_repository_name))/2
        SAMPLEFILE_NUMBER=len(SAMPLE_LIST)
        if (SAMPLEFILE_NUMBER<SAMPLE_NUMBER):
            print("WARNING:More samples in fastq directory than samplefile. Please verify the samplefile.\n")
        if (SAMPLEFILE_NUMBER>SAMPLE_NUMBER):
            print("WARNING:More samples in samplefile than fastq directory.Snakemake will stop.Please verify the samplefile.\n")
        shell("cp {input.raw_fastq_R1} > {output.fastq_R1} ")
        shell("cp {input.raw_fastq_R2} > {output.fastq_R2} ")

rule get_hg19:
    message:
        "download if necessary the hg19 reference genome."
    output:
        hg19_ref = "REF_HG19/hg19.fa"
    shell:
        """
        if [ ! -d REF_HG19 ] ;then 
            mkdir -p REF_HG19 
        fi 
        wget -P REF_HG19/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
        gunzip REF_HG19/hg19.fa.gz       
        #rm -r REF_HG19 
        """

rule prepare_reference:
    message:
        "prepare reference genome before using BMtagger. Executed once per reference genome."
    input:
        HG19 = rules.get_hg19.output.hg19_ref
    output:
        database = expand("REF_HG19/hg19.fa.{ext}", ext=["nhr", "nin", "nsq"]),
        bitmask = "REF_HG19/hg19.bitmask",
        srprism = expand("REF_HG19/hg19.srprism.{subfile}",subfile=['amp','idx','imp','map','mpm','rmp','ss','ssa','ssd'])
    shell:
        """
        bmtool -d {input} -o {output.bitmask} -A 0 -w 18 -z
        srprism mkindex -i {input} -o {output.srprism} -M 7168
        makeblastdb -in {input} -dbtype nucl
        """        

     