############################################
###     Helper functions to compute Div  ###
############################################

#Set of functions to compute the diversity of fasta files processed
#where in each fasta contains the haplotypes present in the sample. (Usually
#each sample have more than 1 fasta, each corresponfing to an amplicon)

files_setup<-function(x){
    
    'Function to order fasta files in sample folders'
    # Arguments:
    #   x   Path to the folder that contains the fasta files
    library(tools)
    setwd(x)
    
    files<-list.files(pattern = "*.fna", recursive = TRUE)
    nt_filename<-gsub(".*/","",files)
    
    samplename<-unique(gsub("\\..*","",nt_filename))
    dir.create("data", showWarnings = FALSE)
    
    for (sample in samplename){
        dir.create(paste("data",sample,sep="/"), showWarnings = FALSE)
        nt_path<-paste("data",sample,"nt",sep="/")
        dir.create(nt_path, showWarnings = FALSE)
        tocopy<-grep(paste(sample,".A",sep=""),files,value=TRUE)
        file.copy(tocopy, nt_path)
    }
}

Div_infolder<-function(folder,dS_size){
    
    'Function that computes the diversity indices and plots the profiles of
    the fastas in the folder, downsampling the abundances'
    # Arguments:
    #   folder  Folder that contains the fasta files to process. (Only 1 sample)
    #   dS_size Downsamplig size.
    
    library(RColorBrewer)
	Comp_div<- function(x,type,dS_size){
	    
		'Function that computes the diversity indices of a fasta'
	    # Arguments:
	    #   x       path to the fasta file
	    #   type    "DNA" or "AA" in case of nt sequence or aa.
	    #   dS_size Downsampling size. 
	    
		AA.dist <-function(seqs){
		        ' Function that computes the distances of a aa sequences'
		        # Arguments:
                #   seqs    Sequences to compute the distances. 
				library(Biostrings)
				library(bios2mds)
		        seqs<-as.character(seqs)
				n <- length(seqs)
				D <- matrix(0,nrow=n,ncol=n)
				rownames(D) <- colnames(D) <- names(seqs)
				mseqs <- t(sapply(seqs,function(sq) toupper(strsplit(sq,split="")[[1]])))
				for(i in 1:(n-1))
				  for(j in (i+1):n)
					D[i,j] <- dis(mseqs[i,],mseqs[j,],"BLOSUM80")
				D <- D + t(D)
				return(D)
			}
								  
		sample_E<-ReadAmplSeqs(x,type=type)
		
		if (length(sample_E$hseqs) > 1) {
		    sample_E$hseqs<-sample_E$hseqs[DSFT(sample_E$nr,dS_size)]
		    sample_E$nr<-sample_E$nr[DSFT(sample_E$nr,dS_size)]
			nm <- nmismatch(pairwiseAlignment(sample_E$hseqs,sample_E$hseqs[1]))
			
			if(type=="DNA"){ 
				dst <- DNA.dist(sample_E$hseqs,model="raw")
			}else{ 
				dst<-AA.dist(sample_E$hseqs)}

			res_sample<-c(length(sample_E$hseqs),SegSites(sample_E$hseqs),TotalMutations(sample_E$hseqs),
							  Shannon(sample_E$nr),GiniSimpson(sample_E$nr),Hill(sample_E$nr,q=0),Hill(sample_E$nr,q=1),
							  Hill(sample_E$nr,q=2),Hill(sample_E$nr,q=3),Hill(sample_E$nr,q=4),Hill(sample_E$nr,q=Inf),
							  Renyi(sample_E$nr,q=0),Renyi(sample_E$nr,q=1),Renyi(sample_E$nr,q=2),Renyi(sample_E$nr,q=3),
							  Renyi(sample_E$nr,q=4),Renyi(sample_E$nr,q=Inf),HCq(sample_E$nr, q=0),HCq(sample_E$nr, q= 1),
							  HCq(sample_E$nr, q= 2),HCq(sample_E$nr, q= 3), HCq(sample_E$nr, q= 4),HCq(sample_E$nr, q= Inf),
							  MutationFreq(dst=dst),FAD(dst),NucleotideDiversity(dst),
							  MutationFreq(nm=nm,nr=sample_E$nr,len=width(sample_E$hseqs)[1]),NucleotideDiversity(dst,sample_E$nr),
							  RaoPow(dst,4,sample_E$nr))
		}else{
		    # In case that the quasispecie only have 1 haplotype some indices are 0. 
			res_sample<-c(length(sample_E$hseqs),SegSites(sample_E$hseqs),TotalMutations(sample_E$hseqs),
							  Shannon(sample_E$nr),GiniSimpson(sample_E$nr),Hill(sample_E$nr,q=0),Hill(sample_E$nr,q=1),
							  Hill(sample_E$nr,q=2),Hill(sample_E$nr,q=3),Hill(sample_E$nr,q=4),Hill(sample_E$nr,q=Inf),
							  Renyi(sample_E$nr,q=0),Renyi(sample_E$nr,q=1),Renyi(sample_E$nr,q=2),Renyi(sample_E$nr,q=3),
							  Renyi(sample_E$nr,q=4),Renyi(sample_E$nr,q=Inf),HCq(sample_E$nr, q=0),HCq(sample_E$nr, q= 1),
							  HCq(sample_E$nr, q= 2),HCq(sample_E$nr, q= 3), HCq(sample_E$nr, q= 4),HCq(sample_E$nr, q= Inf),
							  0,0,0,0,0,0)
		}
	}
	
	Comp_plot<- function(x,type,dS_size){
	    ' Function to plot the profiles of some indices given a fasta file'
	    # Arguments:
	    #   x       path to the fasta file
	    #   type    "DNA" or "AA" in case of nt sequence or aa.
	    #   dS_size Downsampling size. 
	    
		par(mfrow=c(3,1))
		sample_E<-ReadAmplSeqs(x,type=type)
		sample_E$hseqs<-sample_E$hseqs[DSFT(sample_E$nr,dS_size)]
		sample_E$nr<-sample_E$nr[DSFT(sample_E$nr,dS_size)]
		Ampl<-regmatches(x, regexpr("\\.A.", x, perl=TRUE))	
		Ampl<-gsub("\\.","",Ampl)
	    plot(HillProfile(sample_E$nr),type ="b", main=paste("Hill numbers obtained with HillProfile",Ampl, sep="_"))
	    plot(RenyiProfile(sample_E$nr),type ="b", main=paste("Rényi entropy obtained with RenyiProfile",Ampl, sep="_"))
	    plot(HCqProfile(sample_E$nr),type ="b", main=paste("Havrda-Charvat entropy obtained with HCqProfile",Ampl, sep="_"))
	}
	
	mont_plot<-function(x,type,dS_size,Sample){
	    ' Function to plot the montserrat plot of a fasta file'
	    # Arguments:
	    #   x       path to the fasta file
	    #   type    "DNA" or "AA" in case of nt sequence or aa.
	    #   dS_size Downsampling size. 
	    
	    my_purples = brewer.pal(n = 9, "BuPu")[3:9]
	    sample_E<-ReadAmplSeqs(x,type=type)
	    if (length(sample_E$nr)>1){
	        sample_E$hseqs<-sample_E$hseqs[DSFT(sample_E$nr,dS_size)]
	        sample_E$nr<-sample_E$nr[DSFT(sample_E$nr,dS_size)]
	    }
	    lst <- SortByMutations(sample_E$hseqs,sample_E$nr) 
	    qs <- lst$bseqs
	    cnm <- cumsum(table(lst$nm))+1
	    nm.pos <- as.vector(cnm)[-length(cnm)]
	    names(nm.pos) <- names(cnm[-1])
	    ampl<-regmatches(x, regexpr("A[0-9]", x, perl=TRUE))
	    name<-paste(Sample,ampl,sep="_")
	    col<-regmatches(ampl, regexpr("[0-9]", ampl, perl=TRUE))
	    bp<-barplot(prop.table(lst$nr),col= my_purples[as.numeric(col)] )
	    axis(1,at=bp[nm.pos],labels=names(nm.pos),cex.axis=0.7)
	    title(main=name, xlab="Number of mutations", ylab="Abundance") 
	}
	
	type<-ifelse(grepl("Aa+", folder[1]),"AA","DNA")
	Sample<- strsplit(Sample_files[1],"/")[[1]][2]
	DivInd_E<-sapply(folder, function(x) Comp_div(x,type,dS_size))
	colnames(DivInd_E)<-sub(".*/","", colnames(DivInd_E))
	
	DivInd_E<-data.frame(DivInd_E,check.names = FALSE)
	rownames(DivInd_E)<-c("Number of haplotypes","The number of polymorphic sites","Number of mutations","Shannon entropy",
	                      "Gini-Simpson", "Hill numbers_0","Hill numbers_1","Hill numbers_2","Hill numbers_3","Hill numbers_4",
	                      "Hill numbers_Inf","Rényi entropy_0","Rényi entropy_1","Rényi entropy_2","Rényi entropy_3",
	                      "Rényi entropy_4","Rényi entropy_Inf","Havrda-Charvat estimator_0","Havrda-Charvat estimator_1",
	                      "Havrda-Charvat estimator_2","Havrda-Charvat estimator_3","Havrda-Charvat estimator_4",
	                      "Havrda-Charvat estimator_Inf","Average mutation frequency by entity","FAD",
	                      "Nucleotide diversity by entity","Average mutation frequency by molecule",
	                      "Nucleotide diversity","Rao entropy")
	
	DivInd_E$Index_mean<-rowMeans(DivInd_E)
	
	if (type=="DNA"){
	    pdf(file=paste(paste("results","nt", paste(Sample,type,sep="_"),sep="/"),".pdf",sep=""), width = 5, height = 7)
	    sapply(folder, function(x) Comp_plot(x,type,dS_size))
	    dev.off()
	    write.csv(DivInd_E,file=paste(paste("results","nt", paste(Sample,type,sep="_"),sep="/"),"DivInx.csv",sep=""))
	    
	    pdf(file=paste(paste("results","nt", paste(Sample,"Montserrat",type,sep="_"),sep="/"),".pdf",sep=""),
	        width = 10, height = 7,paper="a4r" )
	    par(mfrow = c(2,3))
	    sapply(folder, function(x) mont_plot(x,type,dS_size,Sample))
	    dev.off()
	    
	}else{
	    pdf(file=paste(paste("results","aa", paste(Sample,type,sep="_"),sep="/"),".pdf",sep=""), width = 5, height = 7)
	    sapply(folder, function(x) Comp_plot(x,type,dS_size))
	    dev.off()
	    write.csv(DivInd_E,file=paste(paste("results","aa", paste(Sample,type,sep="_"),sep="/"),"DivInx.csv",sep=""))
	}
	
}