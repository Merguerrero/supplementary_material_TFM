library(ggpubr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

csv_res<-list.files(path= "results/nt", pattern = "*.csv", recursive = TRUE,
                    full.names = TRUE)

if (file.info("data/SampleOrder.txt")$size == 0) {
    filenames<-sub("_.*","",sub(".*/","",csv_res))
} else {
    filenames<-as.vector(read.delim("data/SampleOrder.txt",header=FALSE)[,1])
}

for(csv in csv_res) { 
    x <- read.csv(csv, check.names=FALSE)
    assign(sub("_.*","",sub(".*/","",csv)), x)
}

all_info <- get(filenames[1])[,1]
for (sampl_file in filenames){
    file<-get(sampl_file)[,-1]
    file<-file[,-ncol(file)]
    all_info<-cbind(all_info,file)
}

all_info_t<-transpose(all_info[,-1])
colnames(all_info_t)<-all_info$all_info

samples<-sub("\\..*$","", colnames(all_info))[-1]


all_info_t$sample<-factor( levels=filenames,ordered=TRUE)
all_info_t$ampl<-regmatches(colnames(all_info), regexpr("A[0-9]",  colnames(all_info), perl=TRUE))

### Line plot per amplicon
###############################################

line_plot<- function(ind){
    ' Function to create a linear plot for a given index'
    # Arguments:
    #   ind Diversity index to plot. Same name than in all_info_t data frame. 
    
    my_purples = brewer.pal(n = length(unique(all_info_t$ampl))+3, "BuPu")[-c(1:3)]
    print(ggplot(all_info_t, aes(x=sample, y=all_info_t[,ind], group=ampl)) + 
          geom_line(aes(color=ampl)) +
          geom_point(aes(color=ampl)) +
          scale_colour_manual(values=my_purples) +
          ggtitle(as.character(ind)) + theme(axis.title.y=element_blank())) 
}

pdf(file="report/All_inx_report_nt_Linear.pdf", paper="a4r")
lapply(colnames(all_info_t)[-length(colnames(all_info_t))], function(col) line_plot(col))
dev.off()

