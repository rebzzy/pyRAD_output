########### EDIT BELOW ################
# Assuming you have your pyRAD output file intact and unaltered, set your WD to your pyRAD output directory
wd<-"~/Desktop/Wagtails/pyrad_filter_results/"

# Choose which existing clustering folder you would like to use: 
clust<-"clust.90"

# choose folder name where output files should go
outputdir<-"haplotypes"
consensdir<-"clust.80/" #needs to be named this way to run back through pyRAD; clust.80 corresponds to the number you use in line 10 of the paramaters files

###########
setwd(wd)
dir.create(outputdir)
setwd(paste(wd,outputdir,sep=""))
dir.create(consensdir)
dir.create("outfiles")
dir.create("stats")
setwd(paste(wd,clust,sep=""))
consens_files<-list.files()[grep("consens$",list.files())] 
clustS_files<-list.files()[grep("clustS$",list.files())] 
freqofmins<-NULL
for (i in 1:3){#length(consens_files)){
	consens<-scan(consens_files[i],what="character",sep="\n")
	clustS<-scan(clustS_files[i],what="character",sep="\n")
	output<-paste(paste(wd,outputdir,"/",consensdir,sep=""),consens_files[i],sep="")
	ambig_list<-NULL
	conlength<-NULL # length of consensus that matches ambig_list
	for (j in 1:(length(consens)/2)){ 
		ambig<-grep("[MRWYKSmrwyks]",consens[j*2])
		if (length(ambig)>0){ # if it finds an ambig code
			conlength<-c(conlength,nchar(consens[j*2]))
			ambig_list<-c(ambig_list,(j*2))
		} else { # no ambiguity
			write(paste("allele1_",consens[(j*2)-1],"\n",consens[j*2],sep=""),output,append=TRUE)
			write(paste("allele2_",consens[(j*2)-1],"\n",consens[j*2],sep=""),output,append=TRUE)
		}
	}
	clustEND<-grep("//",clustS) # if ambiguity, search clustS file
	mult_ambig<-NULL
	one_ambig<-NULL
	find_in_clust<-NULL
	for (k in 1:300){#length(ambig_list)){
		multi<-grep("[MRWYKSmrwyks]",strsplit(consens[ambig_list[k]],"")[[1]])
		name<-consens[(ambig_list[k]-1)]
		seq.vec<-strsplit(consens[ambig_list[k]],"")[[1]]
		allele1<-seq.vec
		allele2<-seq.vec
		if (length(multi)>1){ # weird
			lineS<-grep(name,clustS)
			lineE<-clustEND[which(clustEND>lineS)[1]]
			align<-clustS[lineS:lineE] # grab cluster from clustS file
			justseq<-align[seq(2,length(align),2)] # just grabs even elements (aka seq)
			if (nchar(justseq[1])==conlength[k]){ #consens has to be same lenght as clustS
				temp.matrix<-NULL
				for (s in 1:length(justseq)){
					temp.matrix<-rbind(temp.matrix,strsplit(justseq[s],"")[[1]][multi])
					#makes matrix of the variable loci; rows are sequence, col is var site
				}
				for (c in 1:dim(temp.matrix)[2]){
					if (length(which(temp.matrix[,c]=="N"))>0){ #take out rows that have N in them	
						temp.matrix<-as.matrix(temp.matrix[-unique(which(temp.matrix[,c]=="N")),]) 
						}
				}
				if (length(dim(temp.matrix))==0){
					next()
				}
				phase<-as.matrix(unique(temp.matrix[,c(1:length(multi))])) #rows are loci (so read across)
				if (length(phase[,1])==1){
					print("too many Ns") # N counted as SNP
				} else if (length(phase[,1])==2){
					allele1[multi]<-phase[1,]
					allele2[multi]<-phase[2,]
					write(paste("allele1_",name,"\n",paste(allele1,collapse=""),sep=""),output, append=TRUE)
					write(paste("allele2_",name,"\n",paste(allele2,collapse=""),sep=""),output, append=TRUE)
				} else if (length(phase[,1])>2){
					# finds 3 alleles!
						print("Found 3 alleles")
				}			
			} else {
				print("consensus different length from clusters")
			}
		} else if (length(multi)>2) {
				print("three or more ambigs")
		} else if (length(multi)==1){ # ambig occurs once
			if (seq.vec[multi]=="M"|seq.vec[multi]=="m"){
				reps<-c("A","C")
			}
			if (seq.vec[multi]=="R"|seq.vec[multi]=="r"){
				reps<-c("A","G")
			}
			if (seq.vec[multi]=="S"|seq.vec[multi]=="s"){
				reps<-c("G","C")
			}
			if (seq.vec[multi]=="Y"|seq.vec[multi]=="y"){
				reps<-c("T","C")
			}
			if (seq.vec[multi]=="W"|seq.vec[multi]=="w"){
				reps<-c("A","T")
			}
			if (seq.vec[multi]=="K"|seq.vec[multi]=="k"){
				reps<-c("G","T")
			}
			allele1[multi]<-reps[1]
			allele2[multi]<-reps[2]
			write(paste("allele1_",name,"\n",paste(allele1,collapse=""),sep=""),output,append=TRUE)
			write(paste("allele2_",name,"\n",paste(allele2,collapse=""),sep=""),output,append=TRUE)
		}
	}
}
#####################################################
# Run steps 6&7 in pyRAD to get .snps, .fasta etc files

