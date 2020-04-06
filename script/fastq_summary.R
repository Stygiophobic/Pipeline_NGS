library (stringr)

argv <- commandArgs(TRUE)
output=argv[1]

removePath<-function(pathfile){
  split<-str_split(pathfile,"/")
  max<-length(split[[1]])
  samplename<-split[[1]][max]
  return(samplename)
}

#PRE
preqc<-read.table("temp/PREQC_count.txt",sep=";")
colnames(preqc)<-c("PREQC","Sample")
preqc$Sample<-sapply(preqc$Sample,removePath)
preqc$Sample<-str_replace(preqc$Sample,"_R1.fastq","")

#POST
postqc<-read.table("temp/POSTQC_count.txt",sep=";")
colnames(postqc)<-c("POSTQC","Sample")
postqc$Sample<-sapply(postqc$Sample,removePath)
postqc$Sample<-str_replace(postqc$Sample,"_R1_cleaned.fastq","")

#HUMAN
human<-read.table("temp/HUMAN_count.txt",sep=";")
colnames(human)<-c("HUMAN","Sample")
human$Sample<-sapply(human$Sample,removePath)
human$Sample<-str_replace(human$Sample,"_hg19.fastq","")

sum_table<-merge(preqc,human,by="Sample")
sum_table<-merge(sum_table,postqc,by="Sample")

#divide per 4 because => fastq
sum_table$PREQC<-sum_table$PREQC/4
sum_table$HUMAN<-sum_table$HUMAN/8 #/(4*2)> human reads in R1 and R2
sum_table$POSTQC<-sum_table$POSTQC/4
sum_table$PERCENT_HUMAN<-sum_table$HUMAN/sum_table$PREQC*100

write.table(sum_table,output,sep="\t",quote=F,row.names=F)