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

#Check number of sample in the fastq folder
SAMPLE_NUMBER=len(os.listdir(result_repository + fastq_repository_name))/2
SAMPLEFILE_NUMBER=len(SAMPLE_LIST)



print("\nPIPELINE INFORMATION:\n")
if (SAMPLEFILE_NUMBER<SAMPLE_NUMBER):
    print("WARNING:More samples in fastq directory than samplefile. Please verify the samplefile.\n")
if (SAMPLEFILE_NUMBER>SAMPLE_NUMBER):
    print("WARNING:More samples in samplefile than fastq directory.Snakemake will stop.Please verify the samplefile.\n")


#Expected files at the end of the analysis.
rule output_pipeline:
    input:
        #hg19_ref = "REF_HG19/hg19.fa" ,
        #fastq_R1 = expand(result_repository + "FASTQ/{sample}_R1.fastq.gz",sample=SAMPLE_LIST),
        #fastq_R2 = expand(result_repository + "FASTQ/{sample}_R2.fastq.gz",sample=SAMPLE_LIST),
        unzip_fastq_R1 = expand(result_repository + "FASTQ/{sample}_R1.fastq",sample=SAMPLE_LIST),
        unzip_fastq_R2 = expand(result_repository + "FASTQ/{sample}_R2.fastq",sample=SAMPLE_LIST), 
        #human_read_list = expand(result_repository + "BMtagger/{sample}.blacklist" ,sample=SAMPLE_LIST), 
        R1_cleaned =  expand(result_repository + "FASTQ_CLEANED/{sample}_R1_cleaned.fastq",sample=SAMPLE_LIST),
        R2_cleaned =  expand(result_repository + "FASTQ_CLEANED/{sample}_R2_cleaned.fastq",sample=SAMPLE_LIST),
        #HG19_filter =  expand(result_repository + "FASTQ_CLEANED/{sample}_HG19_filter.fastq",sample=SAMPLE_LIST),
        BBMAP = "tool/bbmap/bbsplit.sh",       
        #database = expand("REF_HG19/hg19.fa.{ext}", ext=["nhr", "nin", "nsq"]),
        #bitmask = "REF_HG19/hg19.bitmask",
        #srprism = expand("REF_HG19/hg19.srprism.{subfile}",subfile=['amp','idx','imp','map','rmp','ss','ssa','ssd'])


#control samplefile and fastq's repository contents. Copy procceed fastq.
rule fastq_copy:
    message:
        "Checking samplefile, and decompress fastq."  
    input:
        raw_fastq_R1 = result_repository + fastq_repository_name + "/{sample}_R1.fastq.gz",
        raw_fastq_R2 = result_repository + fastq_repository_name + "/{sample}_R2.fastq.gz",
    params:
        sample=SAMPLE_LIST
    output:
        fastq_R1 = result_repository + "FASTQ/{sample}_R1.fastq.gz",
        fastq_R2 = result_repository + "FASTQ/{sample}_R2.fastq.gz"    
    shell:
        """
        cp {input.raw_fastq_R1}  {output.fastq_R1}
        cp {input.raw_fastq_R2}  {output.fastq_R2}
        """

#Download hg19 if any REF_HG19 repository can be found.
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

#Unzip fastq
rule fastq_unzip:
    message:
        "Unziping file"  
    input:
        copy_fastq_R1 = rules.fastq_copy.output.fastq_R1,
        copy_fastq_R2 = rules.fastq_copy.output.fastq_R2,
    output:
        unzip_fastq_R1 = result_repository + "FASTQ/{sample}_R1.fastq",
        unzip_fastq_R2 = result_repository + "FASTQ/{sample}_R2.fastq"    
    shell:
        """
        gunzip -k {input.copy_fastq_R1} 
        gunzip -k {input.copy_fastq_R2}  
        """
#Download BBmap tool
rule get_BBmap:
    message:"Download tool if necessary."
    output:
        BBMAP = "tool/bbmap/bbsplit.sh"
    shell:
        """
        if [ ! -d tool ] ;then 
            mkdir -p tool 
        fi 
        wget -P tool/ http://downloads.sourceforge.net/project/bbmap/BBMap_38.79.tar.gz 
        tar -C tool/ -xzvf tool/BBMap_38.79.tar.gz 
        chmod +x tool/bbmap/bbsplit.sh
        """        
#Removing human reads from fastq
rule clean_fastq:
    message:
        "Removing human reads from fastq."
    input:
        unzip_fastq_R1 = rules.fastq_unzip.output.unzip_fastq_R1,
        unzip_fastq_R2 = rules.fastq_unzip.output.unzip_fastq_R2,
        HG19 = rules.get_hg19.output.hg19_ref,
        BBMAP = rules.get_BBmap.output.BBMAP
    output:
        R1_cleaned = result_repository + "FASTQ_CLEANED/{sample}_R1_cleaned.fastq",
        R2_cleaned = result_repository + "FASTQ_CLEANED/{sample}_R2_cleaned.fastq", 
    shell:
        """
        {input.BBMAP} in1={input.unzip_fastq_R1} in2={input.unzip_fastq_R2} ref=REF_HG19/hg19.fa \
        basename={wildcards.sample}_%.fastq outu1={output.R1_cleaned} outu2={output.R2_cleaned} path=temp/
        """        

rule trim_fastq:
    message:
        "Use of the trimgalore tool for filtering fastq and get optimal reads."
    input:
        R1_cleaned = rules.clean_fastq.output.R1_cleaned , 
        R2_cleaned = rules.clean_fastq.output.R2_cleaned ,          

            