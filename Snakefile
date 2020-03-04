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

#WILCARD SUBTYPE:
SUBTYPE=['BVIC_Malaysia2506','BYAM_Florida4','pH1N1_California07','H3N2_Perth16']
dict_subtype={}

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
        #unzip_fastq_R1 = expand(result_repository + "FASTQ/{sample}_R1.fastq",sample=SAMPLE_LIST),
        #unzip_fastq_R2 = expand(result_repository + "FASTQ/{sample}_R2.fastq",sample=SAMPLE_LIST), 
        #R1_trimmed =  expand(result_repository + "FASTQ_TRIMM/{sample}_R1_trimmed.fastq",sample=SAMPLE_LIST) ,
        #R2_trimmed =  expand(result_repository + "FASTQ_TRIMM/{sample}_R2_trimmed.fastq",sample=SAMPLE_LIST) ,
        #SAM = expand(result_repository + "SAM/{sample}.sam",sample=SAMPLE_LIST),
        count_premapping = expand(result_repository + "COUNT_MAPPING/{sample}_premapping.csv",sample=SAMPLE_LIST) ,
        #R1_cleaned =  expand(result_repository + "FASTQ_CLEANED/{sample}_R1_cleaned.fastq",sample=SAMPLE_LIST),
        #R2_cleaned =  expand(result_repository + "FASTQ_CLEANED/{sample}_R2_cleaned.fastq",sample=SAMPLE_LIST),
        #T_G = "tool/TrimGalore-0.6.5/trim_galore" ,
        #BBMAP = "tool/bbmap/bbmap.sh",       



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
    message:
        "Download tool if necessary."
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
        rm tool/BBMap_38.79.tar.gz
        """        
rule create_index:
    message:
        "Produce index on hg19 reference."
    params:
        path_human= result_repository 
    input:
        BBMAP = rules.get_BBmap.output.BBMAP,
        HG19 = rules.get_hg19.output.hg19_ref
    output:
        index = "temp/ref/genome/1/summary.txt"
    shell:
        """
        tool/bbmap/bbsplit.sh  ref=REF_HG19/hg19.fa path=temp/
        #tool/bbmap/bbmap.sh ref=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz  \
        #-Xmx16g -usemodulo=true path=temp/
        """

#Removing human reads from fastq
rule clean_fastq:
    message:
        "Removing human reads from fastq."   
    #threads:1        
    resources: mem_gb= 20
    input:
        unzip_fastq_R1 = rules.fastq_unzip.output.unzip_fastq_R1,
        unzip_fastq_R2 = rules.fastq_unzip.output.unzip_fastq_R2,
        BBMAP = rules.get_BBmap.output.BBMAP    ,
        index = rules.create_index.output.index    
    output:
        R1_cleaned = result_repository + "FASTQ_CLEANED/{sample}_R1_cleaned.fastq",
        R2_cleaned = result_repository + "FASTQ_CLEANED/{sample}_R2_cleaned.fastq",
        HUMAN = result_repository + "FASTQ_CLEANED/{sample}_hg19.fastq"
    params:
        path_human= result_repository + "FASTQ_CLEANED/"   
    shell:
        """
        {input.BBMAP} in1={input.unzip_fastq_R1} in2={input.unzip_fastq_R2}  \
        basename={rules.clean_fastq.params.path_human}{wildcards.sample}_%.fastq outu1={output.R1_cleaned} outu2={output.R2_cleaned} \
        path=temp/ 
        """        

#Downloading the trimgalore binary if needed
rule get_trimgalore:
    message:
        "download Trimgalore binary. cutadapt and fastqc are already installed into the singularity IMG."
    output:
        T_G = "tool/TrimGalore-0.6.5/trim_galore"
    shell:
        """
        if [ ! -d tool ] ;then 
            mkdir -p tool 
        fi 
        curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o tool/trim_galore.tar.gz
        tar -C tool/ -xvzf tool/trim_galore.tar.gz
        #chmod +x tool/TrimGalore-0.6.5/trim_galore
        rm tool/trim_galore.tar.gz
        """        

#Perform trimming on all cleaned fastq using TrimGalore
rule trim_fastq:
    message:
        "Use of the trimgalore tool for filtering fastq and get optimal reads."
    input:
        R1_cleaned = rules.clean_fastq.output.R1_cleaned , 
        R2_cleaned = rules.clean_fastq.output.R2_cleaned ,       
        T_G = rules.get_trimgalore.output.T_G
    output:
        R1_trimmed =  result_repository + "FASTQ_TRIMM/{sample}_R1_trimmed.fastq" ,
        R2_trimmed =  result_repository + "FASTQ_TRIMM/{sample}_R2_trimmed.fastq" ,
    params:
        TG_output = result_repository + "FASTQ_TRIMM/"        
    shell:
        """
        {input.T_G} --dont_gzip --no_report_file --trim-n --quality 20 --length 50 --basename {wildcards.sample} \
        --paired {input.R1_cleaned} {input.R2_cleaned} --output_dir {rules.trim_fastq.params.TG_output}
        mv {rules.trim_fastq.params.TG_output}{wildcards.sample}_val_1.fq {rules.trim_fastq.params.TG_output}{wildcards.sample}_R1_trimmed.fastq
        mv {rules.trim_fastq.params.TG_output}{wildcards.sample}_val_2.fq {rules.trim_fastq.params.TG_output}{wildcards.sample}_R2_trimmed.fastq
        """                   

rule premapping_align:
    message:
        "Alignment on each influenza subtype to the sample using bwa."
    threads:4    
    input:
        R1_trimmed = rules.trim_fastq.output.R1_trimmed ,
        R2_trimmed = rules.trim_fastq.output.R2_trimmed ,
    output:
        SAM = result_repository + "SAM/{sample}.sam"
    shell:
        """
        bwa index -p temp/influenza_subtype mapping/pre_mapping/influenza_subtype.fasta
        bwa mem -t 4 -O 10 -E 2 temp/influenza_subtype {input.R1_trimmed} {input.R2_trimmed} > {output}
        """

rule premapping_count:  
    message:
        "Counting and sum up the mapping for each samples using samtools and unix commands."
    input: 
        SAM = rules.premapping_align.output.SAM
    output:  
        count_premapping = result_repository + "COUNT_MAPPING/{sample}_premapping.csv"
    run:
        shell("samtools view -S {input} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | tr '\t' ';'> temp/count_{wildcards.sample}.txt")
        table=pd.read_csv("temp/count_"+wildcards.sample+".txt",sep=";",names=['COUNT','MAPPING'])  
        table = table.loc[table['MAPPING'] != "*"]
        SubType = table['MAPPING'].values[0]
        result=open(result_repository + "COUNT_MAPPING/"+wildcards.sample+"_premapping.csv",'w')
        result.write(wildcards.sample + ";" + SubType)
        result.close()

        


            