###########################################
###         Diversity with DSFT         ###
###########################################

WD<-"/home/mervhir/Escriptori/MEGA/TFM/Diversity"
setwd(WD)

library(QSutils)
source("DivFun.R")

### Set up the data folder
#############################

files_setup(getwd())

NTfiles<-Sys.glob("data/*/nt/*.fna")

### NT Smaller sample size in the comparison
###############################################

NThaplo<-lapply(NTfiles,function(x) ReadAmplSeqs(x))
names(NThaplo)<-NTfiles

len_NT<-lapply(NThaplo, function(x) sum(x$nr))

# Set the downsampling size.
Size_DS<-len_NT[min(unlist(len_NT))==len_NT][[1]][1]

SampleNames<-list.dirs(path="./data",full.names = FALSE,recursive=FALSE)

for (Sample in SampleNames){
    cat(paste("Sample: ",Sample, "is processed \n"))
    Sample_files<-grep(paste(Sample,"/",sep=""),NTfiles,value=TRUE)
    Div_infolder(Sample_files,Size_DS)
}

### Report of each index among samples
##############################################

source("Report.R")


