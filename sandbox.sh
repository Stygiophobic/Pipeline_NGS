#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
heure=$(date +%H%M)
jour=$(date +%Y%m%d) 
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/pipeline.simg"
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

#Tri Par R des fichiers VCF

cp temp/S20121_S52.fasta /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/VARCALL/

samtools mpileup -u -d 1000 -f mapping/subtype_mapping/BVIC_Malaysia2506.fasta /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/BAM_SUBTYPE/S20119_S50.bam \
| bcftools call --ploidy 1 -c | vcfutils.pl vcf2fq | seqtk seq -a -  > mdr.txt

#Bedtools coverage
bedtools genomecov -ibam S20121_S52_sorted.bam -d -strand -  > COV_S20121_S52_reverse.txt
bedtools genomecov -ibam S20121_S52_sorted.bam -d -strand +  > COV_S20121_S52_forward.txt


cd /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/VARCALL/
sed '1,10d' S20121_S52.vcf > varfile.vcf

S20076_S3
S20117_S41
S20106_S16
S20119_S50
S20186_S77
S20107_S17
S20112_S29

