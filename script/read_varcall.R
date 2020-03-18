library(stringr)
library(seqinr)

#setwd("C:/Users/regueex/Desktop/Pipeline_VARCALL")
argv <- commandArgs(TRUE)


#Parsing data
#VCF<-read.table("test.vcf",sep="\t",header=T)
#fasta<-read.table("S20121_S52.fasta")
#samplesheet<-read.csv2("Samplefile_test.csv")
#sample<-"S20121_S52"
#output<-
VCF<-read.table(argv[1],sep="\t",header=T)
fasta<-read.table(argv[2])
samplesheet<-read.csv2(argv[3])
sample<-argv[4]
output<-argv[5]


#Prepare VCF

#getAC
getAC<-function(value_info){
  split<-str_split(value_info,";")
  string_AC<-split[[1]][1]
  split2<-str_split(string_AC,"=")
  if (split2[[1]][2]==".") value_AC_num<-0
  else value_AC_num=as.numeric(split2[[1]][2])
  return(value_AC_num)
}

#getAF
getAF<-function(value_info){
  split<-str_split(value_info,";")
  string_AF<-split[[1]][2]
  split2<-str_split(string_AF,"=")
  if (split2[[1]][2]==".") value_AF_num<-0
  else value_AF_num=as.numeric(split2[[1]][2])
  return(value_AF_num)
}

VCF$AC<-sapply(as.character(VCF$INFO),getAC)
VCF$AF<-sapply(as.character(VCF$INFO),getAF)

VCF<-subset(VCF,QUAL>100000)
VCF<-VCF[,c(1,2,4,5)]

len_fasta<-nrow(fasta)/2
vec_header<-c()
vec_seq<-c()

vec_header[1]<-as.character(fasta[1,])
for (i in 1:(len_fasta-1)){
  row_header<-1+i*2
  vec_header[i+1]<-as.character(fasta[row_header,])  
}

for (i in 1:(len_fasta)){
  row_seq<-0+i*2
  vec_seq[i]<-as.character(fasta[row_seq,])  
}

table_seq<-cbind(vec_header,vec_seq)
table_seq<-as.data.frame(table_seq)

colnames(table_seq)<-c("Header","Seq")
table_seq$Header<-str_replace(table_seq$Header,">","")
table_seq$Seq<-as.character(table_seq$Seq)
table_seq$Header<-as.character(table_seq$Header)

vec_seq_vcf<-c()

for (i in 1:nrow(table_seq)){
  #grab header and sequence
  #i<-4
  header<-as.character(table_seq[i,1])
  sequence<-as.character(table_seq[i,2])
  sequence<-str_split(sequence,"")
  sequence<-unlist(sequence)
  if (!(header %in% VCF$CHROM)) {
    vec_seq_vcf[i]<-paste(sequence,collapse = '')
  }else{
    new_sequence<-c()
    cmpt_new_seq<-1
    position<-0 #0
    while (position<length(sequence)){ #length(sequence)
      position<-position+1
      #test
      #position<-561
      info_VCF<-subset(VCF,CHROM==header & POS==position)
      if (nrow(info_VCF)==1){
        print("Variant Valide")
        REF<-as.character(info_VCF$REF)
        ALT<-as.character(info_VCF$ALT)
        #if deletion
        gap_ref<-nchar(REF)-1
        if (gap_ref>1) position<-position+gap_ref
        #if insertion
        gap_alt<-nchar(ALT)-1
        if (gap_alt>1) position<-position+gap_alt
        new_sequence[cmpt_new_seq]<-ALT
        cmpt_new_seq<-cmpt_new_seq+1
        next
        }
      new_sequence[cmpt_new_seq]<-sequence[position]
      cmpt_new_seq<-cmpt_new_seq+1
    }
    #Concatenate all nucleotides
    new_sequence<-paste(new_sequence,collapse = '')
    vec_seq_vcf[i]<-new_sequence
  }
}

table_seq$new_seq<-vec_seq_vcf
table_seq$new_seq<-as.character(table_seq$new_seq)

#Load
translation_table<-read.csv2("data/translation_data.csv")
table_seq<-merge(table_seq,translation_table,by="Header")

test_sequence<-table_seq[4,3]
header<-table_seq[4,1]

translateSEQ<-function(header){
  seq_to_tr<-table_seq$new_seq[table_seq$Header==header]
  start_trans<-as.numeric(table_seq$START1[table_seq$Header==header])
  stop_trans<-as.numeric(table_seq$STOP1[table_seq$Header==header])
  seq_to_tr<-str_split(seq_to_tr,"")
  seq_to_tr<-unlist(seq_to_tr)
  cut_seq<-seq_to_tr[start_trans:stop_trans]
  prot_seq<-translate(seq=tolower(cut_seq),frame=0,numcode = 1,sens="F")
  stop_pos<-which(prot_seq=="*")-1
  prot_seq<-prot_seq[1:stop_pos]
  final_prot_seq<-paste(prot_seq,collapse = '')
  return(final_prot_seq)
}

table_seq$prot_seq<-sapply(table_seq$Header,translateSEQ)

new_header<-as.character(samplesheet$FASTA_OUTPUT[samplesheet$SAMPLE==sample])


split_seg<-function(header){
  header_split<-str_split(header,"_")
  segment<-as.character(header_split[[1]][2])
  return(segment)
}

table_seq$seg_vec<-sapply(table_seq$Header,split_seg)
table_seq$new_header<-paste0(">",new_header,"_",table_seq$seg_vec)


#Nucleic sequences
if (nrow(table_seq)==8) {
  vec_result<-c(table_seq$new_header[1],table_seq$new_seq[1],
                table_seq$new_header[2],table_seq$new_seq[2],
                table_seq$new_header[3],table_seq$new_seq[3],
                table_seq$new_header[4],table_seq$new_seq[4],
                table_seq$new_header[5],table_seq$new_seq[5],
                table_seq$new_header[6],table_seq$new_seq[6],
                table_seq$new_header[7],table_seq$new_seq[7],
                table_seq$new_header[8],table_seq$new_seq[8])
  write.table(vec_result,file=output,quote = FALSE,row.names = FALSE)
}

if (nrow(table_seq)==7) {
  vec_result<-c(table_seq$new_header[1],table_seq$new_seq[1],
                table_seq$new_header[2],table_seq$new_seq[2],
                table_seq$new_header[3],table_seq$new_seq[3],
                table_seq$new_header[4],table_seq$new_seq[4],
                table_seq$new_header[5],table_seq$new_seq[5],
                table_seq$new_header[6],table_seq$new_seq[6],
                table_seq$new_header[7],table_seq$new_seq[7])
  write.table(vec_result,file=output,quote = FALSE,row.names = FALSE) 
}
