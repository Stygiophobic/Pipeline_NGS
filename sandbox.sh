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

#BWA
bwa index -p temp/influenza_subtype mapping/pre_mapping/influenza_subtype.fasta
#SAMTOOLS
samtools view -bSc /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/SAM/S20119_S50.sam > /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/SAM/test.bam
samtools view -S /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/SAM/S20119_S50.sam | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | tr '\t'  ';'> RefsReadsCount.txt