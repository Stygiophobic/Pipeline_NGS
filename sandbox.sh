#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
heure=$(date +%H%M)
jour=$(date +%Y%m%d) 
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/NGS-PIPELINE.simg"
fastq_repository="200218_NB501048_0691_AHYTKHAFXY_1582192205"
run_name="Test_OLD_DATA"
working_repository="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/"
pathdata=$working_repository$fastq_repository
samplefile="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/Samplefile_test.csv"
################################################################################

################################################################################
#########################       PRE-PROCESS      ###############################
################################################################################ 
ln -s /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/${fastq_repository}/ViroEst-Routine $pathdata



################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################ 
k5start -U -f /home/chu-lyon.fr/regueex/login.kt -- nohup 
singularity exec $singularity_img snakemake \
    --resources mem_gb=32 \
    --config Result_Repository=$working_repository \
             Project_folder=$fastq_repository \
             Samplesheet_Location=$samplefile > $working_repository${rep_report}report_${jour}_${heure}.txt    
                        
################################################################################
###option:
# --dryrun => fait tourner le pipeline à vide pour controler

#Create symlink on data
ln -s /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/${fastq_repository}/ViroEst-Routine $pathdata

echo "kek"

#Remove symlink on data
unlink $pathdata


#test BBtool
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/NGS-PIPELINE.simg"
singularity shell $singularity_img

#full job
tool/bbmap/bbsplit.sh in1=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ/S20105_S15_R1.fastq in2=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ/S20105_S15_R2.fastq ref=REF_HG19/hg19.fa         basename=S20105_S15_%.fastq outu1=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ_CLEANED/S20105_S15_R1_cleaned.fastq outu2=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ_CLEANED/S20105_S15_R2_cleaned.fastq path=temp/ -Xmx16g
#index:
tool/bbmap/bbsplit.sh  ref=REF_HG19/hg19.fa  basename=S20105_S15_%.fastq outu1=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ_CLEANED/S20105_S15_R1_cleaned.fastq outu2=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ_CLEANED/S20105_S15_R2_cleaned.fastq path=temp/ -Xmx16g
tool/bbmap/bbsplit.sh  ref=REF_HG19/hg19.fa path=temp/

tool/bbmap/bbsplit.sh in1=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ/S20105_S15_R1.fastq in2=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ/S20105_S15_R2.fastq ref=REF_HG19/hg19.fa         basename=S20105_S15_%.fastq outu1=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ_CLEANED/S20105_S15_R1_cleaned.fastq outu2=/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/FASTQ_CLEANED/S20105_S15_R2_cleaned.fastq path=temp/ 

