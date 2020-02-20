#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/IMG_SINGULARITY/PIPELINE-NGS.simg"
fastq_repository="191217_NB501048_0665_AHWM32AFXY_1576689602"
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
singularity exec $singularity_img snakemake \
    --dryrun \
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