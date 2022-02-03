
#read in output from cluster, as raw up or dntag counts
dn<-read.csv('His_LP_S25_dntag_rawcounts.csv')
#isolate db names (from column)
colnames(dn)<-gsub('\\.', '-', colnames(dn))


#import summary file to convert dn_db to plate identifier to be used for analysis with old R script 
db_summary<-read.csv('DB_barcodes_V5.csv', sep = '\t')
ad_summary<-read.csv('AD_barcodes_V5.csv', sep = '\t')




#for loop to replace column names with plate identifiers

#this line just fixes the old matrix outputs from guru that still had _ in AVR58 gene, should be required when using latest output from guru on 1-15-19  
#IMPORTANT
#colnames(dn)[853:896]<-sub('_', '',colnames(dn)[853:896])
#colnames(dn)[966]<-'LPG0439-02-349_BC-1'
#-----------------
dn_names<-colnames(dn) 



for (i in 1:length(dn_names)){
  if (dn_names[i] %in% db_summary$Locus){
    dn_names[i] <- as.character(db_summary$identifier[dn_names[i] == db_summary$Locus])
  }
}

#remove X at begining
dn_names<-dn_names[2:length(dn_names)]


#for loop to replace column names with plate identifiers

ad_names<-as.character(dn$X)


for (i in 1:length(ad_names)){
  if (ad_names[i] %in% ad_summary$Locus){
    ad_names[i]<-as.character(ad_summary$Identifier[ad_names[i]==ad_summary$Locus])
  }
}



#--------------now make a new matrix that has the appropriate row and column names
rownames(dn)<-ad_names
#remove old names column
dn$X<-NULL

colnames(dn)<-dn_names
  
#this line just fixes the old matrix outputs from guru that still had _ in AVR58 gene, should be required when using latest output from guru on 1-15-19  
write.csv(dn, 'plus_his_dn.csv')

barcodes.file <- read.csv('res/barcodes_V2.csv', header = F)

#code to find missing barcodes that are in my matrix, but not the barcodes.csv file
#dn_names[!(dn_names %in% as.character(barcodes.file$V6))]
dn_names[!(dn_names %in% as.character(barcodes.file$V6))]
as.character(barcodes.file$V6)[!(as.character(barcodes.file$V6[barcodes.file$V1=='DB']) %in% dn_names)]



DB_barcodes<-as.character(barcodes.file$V6[barcodes.file$V1=='DB'])

dn_names[(!dn_names %in% DB_barcodes)]

DB_barcodes[!(DB_barcodes %in% dn_names)]

'DBPlatinum1A16' %in% dn_names




AD_barcodes<-as.character(barcodes.file$V6[barcodes.file$V1=='AD'])
ad_names[!(ad_names %in% AD_barcodes)]

AD_barcodes[!(AD_barcodes %in% ad_names)]
#----------------------now repeat with up plus his counts
