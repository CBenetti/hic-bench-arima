# Load the libraries
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(tidyr)

# Read arguments
args <- commandArgs(TRUE)
branch <- sprintf("%s", args[1])
output <- sprintf("%s", args[2])
op <- paste("find ",branch," -name 'Arima_QC_deep.txt'",sep="")
file_list <- system(op,intern=T)
tab<-do.call("rbind",lapply(file_list,function(x){read.table(x,sep="\t",header=T)}))
df<-data.frame(Sample=tab[,1],intra=tab[,14],inter=tab[,18],PCR_dup=tab[,5],Optical=tab[,7],Unmapped=(tab[,2]-tab[,3]))
long_df<-df%>%pivot_longer(cols = -Sample,names_to="Class",values_to="Pairs")
pdf(paste(output,"/QC_plots.pdf",sep=""))
plot(ggplot(long_df, aes(fill=Class, y=Pairs, x=Sample)) + geom_bar(position="stack", stat="identity")+scale_fill_brewer(palette="Spectral")+
scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"))+coord_flip()+ggtitle("Pair count - Arima QC")+theme(legend.position = "top",legend.background = element_rect(fill = "#F0F0F0")))
plot(ggplot(long_df, aes(fill=Class, y=Pairs, x=Sample)) + geom_bar(position="fill", stat="identity")+
scale_fill_brewer(palette="Spectral")+coord_flip()+ggtitle("Percentages (%) - Arima QC")+theme(legend.position = "top",legend.background = element_rect(fill = "#F0F0F0")))
dev.off()

