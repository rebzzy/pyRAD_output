####################################################################################################
#### make structure or geneland input file from allelic data (non-binary) ##########################
insnp<-"~/Desktop/Wagtails/pyrad_filter_results/haplotypes/outfiles/whites_haplo2_cl80c12m12.snps" # unedited output from pyRAD
nspecies<-23
maxallele<-9 #maximum number of alleles at each locus
out.suffix<-"structure"
####################################################################################################
####################################################################################################

a<-read.table(insnp,header=FALSE,skip=1,row.names=1,na.strings="000")
sp.names<-gsub("allele1_","",row.names(a)[1:nspecies])
strux_table<-NULL
for (i in 1:length(a[1,])){
	if (a[1,i]=="_"){ next()}
	rem.N<-a[-(grep("N",a[,i])),i] #vector of snps at locus but without "N"
	if (length(rem.N)==0){ 
		rem.N<-a[,i]
	} 
	list.geno<-unique(rem.N)
	if (length(list.geno)>maxallele){  next()	}
	one_locus<-NULL
	for (r in 1:nspecies){
		if (length(grep("N",a[r,i]))>0){
			geno<-"0 0"
		} else {
			geno<-paste(which(list.geno==a[r,i]),which(list.geno==a[(r+nspecies),i]),sep=" ")
			if (nchar(geno)==1){
				geno<-paste(geno,"N",sep=" ")
			}
			if (geno==""){
				print(i)
			}
		}
		one_locus<-rbind(one_locus,geno)
	}
	strux_table<-cbind(strux_table,one_locus)
}
write.table(cbind(sp.names,strux_table),gsub("snps","test",insnp),row.names=FALSE,quote=FALSE,sep=" ",col.names=FALSE)
