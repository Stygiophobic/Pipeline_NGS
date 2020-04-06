library(stringr)

count_seg<-function(table){
  vec_result<-c()
  #count UNK
  if ("*" %in% table$Segment) vec_result[1]<-table$Count[table$Segment=="*"] else vec_result[1]<-0
  #count BVIC
  if ("BVIC_S1" %in% table$Segment) vec_result[2]<-table$Count[table$Segment=="BVIC_S1"] else vec_result[2]<-0
  if ("BVIC_S2" %in% table$Segment) vec_result[3]<-table$Count[table$Segment=="BVIC_S2"] else vec_result[3]<-0
  if ("BVIC_S3" %in% table$Segment) vec_result[4]<-table$Count[table$Segment=="BVIC_S3"] else vec_result[4]<-0
  if ("BVIC_S4" %in% table$Segment) vec_result[5]<-table$Count[table$Segment=="BVIC_S4"] else vec_result[5]<-0
  if ("BVIC_S5" %in% table$Segment) vec_result[6]<-table$Count[table$Segment=="BVIC_S5"] else vec_result[6]<-0
  if ("BVIC_S6" %in% table$Segment) vec_result[7]<-table$Count[table$Segment=="BVIC_S6"] else vec_result[7]<-0
  if ("BVIC_S7" %in% table$Segment) vec_result[8]<-table$Count[table$Segment=="BVIC_S7"] else vec_result[8]<-0
  if ("BVIC_S8" %in% table$Segment) vec_result[9]<-table$Count[table$Segment=="BVIC_S8"] else vec_result[9]<-0
  #count BYAM
  if ("BYAM_S1" %in% table$Segment) vec_result[10]<-table$Count[table$Segment=="BYAM_S1"] else vec_result[10]<-0
  if ("BYAM_S2" %in% table$Segment) vec_result[11]<-table$Count[table$Segment=="BYAM_S2"] else vec_result[11]<-0
  if ("BYAM_S3" %in% table$Segment) vec_result[12]<-table$Count[table$Segment=="BYAM_S3"] else vec_result[12]<-0
  if ("BYAM_S4" %in% table$Segment) vec_result[13]<-table$Count[table$Segment=="BYAM_S4"] else vec_result[13]<-0
  if ("BYAM_S5" %in% table$Segment) vec_result[14]<-table$Count[table$Segment=="BYAM_S5"] else vec_result[14]<-0
  if ("BYAM_S6" %in% table$Segment) vec_result[15]<-table$Count[table$Segment=="BYAM_S6"] else vec_result[15]<-0
  if ("BYAM_S7" %in% table$Segment) vec_result[16]<-table$Count[table$Segment=="BYAM_S7"] else vec_result[16]<-0
  if ("BYAM_S8" %in% table$Segment) vec_result[17]<-table$Count[table$Segment=="BYAM_S8"] else vec_result[17]<-0
  #count H1N1
  if ("H1N1_S1" %in% table$Segment) vec_result[18]<-table$Count[table$Segment=="H1N1_S1"] else vec_result[18]<-0
  if ("H1N1_S2" %in% table$Segment) vec_result[19]<-table$Count[table$Segment=="H1N1_S2"] else vec_result[19]<-0
  if ("H1N1_S3" %in% table$Segment) vec_result[20]<-table$Count[table$Segment=="H1N1_S3"] else vec_result[20]<-0
  if ("H1N1_S4" %in% table$Segment) vec_result[21]<-table$Count[table$Segment=="H1N1_S4"] else vec_result[21]<-0
  if ("H1N1_S5" %in% table$Segment) vec_result[22]<-table$Count[table$Segment=="H1N1_S5"] else vec_result[22]<-0
  if ("H1N1_S6" %in% table$Segment) vec_result[23]<-table$Count[table$Segment=="H1N1_S6"] else vec_result[23]<-0
  if ("H1N1_S7" %in% table$Segment) vec_result[24]<-table$Count[table$Segment=="H1N1_S7"] else vec_result[24]<-0
  if ("H1N1_S8" %in% table$Segment) vec_result[25]<-table$Count[table$Segment=="H1N1_S8"] else vec_result[25]<-0
  #count H3N2
  if ("H3N2_S1" %in% table$Segment) vec_result[26]<-table$Count[table$Segment=="H3N2_S1"] else vec_result[26]<-0
  if ("H3N2_S2" %in% table$Segment) vec_result[27]<-table$Count[table$Segment=="H3N2_S2"] else vec_result[27]<-0
  if ("H3N2_S3" %in% table$Segment) vec_result[28]<-table$Count[table$Segment=="H3N2_S3"] else vec_result[28]<-0
  if ("H3N2_S4" %in% table$Segment) vec_result[29]<-table$Count[table$Segment=="H3N2_S4"] else vec_result[29]<-0
  if ("H3N2_S5" %in% table$Segment) vec_result[30]<-table$Count[table$Segment=="H3N2_S5"] else vec_result[30]<-0
  if ("H3N2_S6" %in% table$Segment) vec_result[31]<-table$Count[table$Segment=="H3N2_S6"] else vec_result[31]<-0
  if ("H3N2_S7" %in% table$Segment) vec_result[32]<-table$Count[table$Segment=="H3N2_S7"] else vec_result[32]<-0
  if ("H3N2_S8" %in% table$Segment) vec_result[33]<-table$Count[table$Segment=="H3N2_S8"] else vec_result[33]<-0
  #Compute percent
  sum_count<-sum(vec_result)
  vec_result<-vec_result/sum_count*100
  return(vec_result)
}

argv <- commandArgs(TRUE)
path_input=argv[1]
output=argv[2]

label<-c("Segment","Count")

col_mat<-c("UKN","BVIC_S1","BVIC_S2","BVIC_S3","BVIC_S4","BVIC_S5","BVIC_S6","BVIC_S7","BVIC_S8",
           "BYAM_S1","BYAM_S2","BYAM_S3","BYAM_S4","BYAM_S5","BYAM_S6","BYAM_S7","BYAM_S8",
           "H1N1_S1","H1N1_S2","H1N1_S3","H1N1_S4","H1N1_S5","H1N1_S6","H1N1_S7","H1N1_S8",
           "H3N2_S1","H3N2_S2","H3N2_S3","H3N2_S4","H3N2_S5","H3N2_S6","H3N2_S7","H3N2_S8")

files<-list.files(path_input)

sample_lab<-c()
cpt<-1
all_result<-list()

for (file in files) {
  sample<-str_replace(file,"_count.tsv","")
  recomb_table<-read.table(paste0(path_input,file),header=F)
  colnames(recomb_table)<-label
  sample_lab[cpt]<-sample
  result<-count_seg(recomb_table)
  result<-as.numeric(result)
  all_result[[cpt]]<-result
  cpt<-cpt+1
}

recomb_mat<-do.call(rbind.data.frame, all_result)
rownames(recomb_mat)<-sample_lab
colnames(recomb_mat)<-col_mat

recomb_mat<-as.matrix(recomb_mat)

write.table(recomb_mat,output,sep="\t",quote=F)