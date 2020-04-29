#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
heure=$(date +%H%M)
jour=$(date +%Y%m%d) 
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/pipeline.simg"
fastq_repository="200415_NB501048_0715_AH2GJWAFX2_1587050402"
run_name="FULL_RUN"
working_repository="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/"
pathdata=$working_repository$fastq_repository
samplefile="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/sample_fullrun.csv"
################################################################################

################################################################################
#########################       PRE-PROCESS      ###############################
################################################################################ 
#Create symlink on fastq
ln -s /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/${fastq_repository}/ViroEst-Routine $pathdata

################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################ 
k5start -U -f /home/chu-lyon.fr/regueex/login.kt -- nohup singularity exec $singularity_img snakemake \
    --resources mem_gb=32 \
    --config Result_Repository=$working_repository \
             Project_folder=$fastq_repository \
             Samplesheet_Location=$samplefile > $working_repository${rep_report}report_${jour}_${heure}.txt & 

singularity exec $singularity_img snakemake --unlock --config Result_Repository=$working_repository \
             Project_folder=$fastq_repository \
             Samplesheet_Location=$samplefile


####TEST Gestion TEMOIN#####

k5start -U -f /home/chu-lyon.fr/regueex/login.kt -- nohup singularity exec $singularity_img snakemake -s Trimming_QCcheck \
            --cores 8 \
            --resources mem_gb=32 \
            --config Result_Repository=$working_repository \
                Project_folder=$fastq_repository \
                Samplesheet_Location=$samplefile

k5start -U -f /home/chu-lyon.fr/regueex/login.kt -- nohup singularity exec $singularity_img snakemake -s 2_mappingsubype \
            --cores 8 \
            --resources mem_gb=32 \
            --config Result_Repository=$working_repository \
                Project_folder=$fastq_repository \
                Samplesheet_Location=$samplefile

################################################################################
###option:
# --dryrun => fait tourner le pipeline à vide pour controler

#Create symlink on data
ln -s /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/${fastq_repository}/ViroEst-Routine $pathdata


#Remove symlink on data
unlink $pathdata

#test cmd count

singularity shell "/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/pipeline.simg"
