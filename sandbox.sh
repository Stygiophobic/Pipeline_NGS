#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
heure=$(date +%H%M)
jour=$(date +%Y%m%d) 
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/NGS-PIPELINE.simg"
fastq_repository="200218_NB501480_0458_AH3FM2BGXB_1582206003"
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
    --config Result_Repository=$working_repository \
             Project_folder=$fastq_repository \
             Samplesheet_Location=$samplefile > $working_repository${rep_report}report_${jour}_${heure}.txt    
                        
################################################################################
###option:
# --dryrun => fait tourner le pipeline à vide pour controler

#Create symlink on data
ln -s /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/${fastq_repository}/ViroEst-Routine $pathdata


#Remove symlink on data
unlink $pathdata

