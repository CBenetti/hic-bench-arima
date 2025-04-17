args <- commandArgs(trailingOnly = TRUE)
res <- as.numeric(ifelse(length(args) > 0, args[1], "my_environment.RData"))
#Get all peak files
files <- system("find -L ../../MAPS/results/ -name '*_peaks.broadPeak'",intern=T) ##consider changing this if Arima Maps is not used
#Get and save coordinates for them
total_peaks <- lapply(files,function(x){
        f <- read.table(x,sep="\t",header=F)
        f <- f[order(f[,8],decreasing=T),] ##consider changing this if Arima Maps is not used
        f[c(1:dim(f)[1]),] ##first evaluate minimum peak number
})
if(length(total_peaks)>1){
        total_peaks <- do.call("rbind",total_peaks)
}else{
        total_peaks <- total_peaks[[1]]
}
total_peaks[,4] <- sprintf(paste0("%0", nchar(as.character(dim(total_peaks)[1])), "d"), 1:dim(total_peaks)[1])
write.table(total_peaks,file="MAPS_peaks.broadPeak",sep="\t",row.names=F,col.names=F,quote=F) ##consider changing this if Arima Maps is not us>
#First process the peaks as to bin them in 1Kb bins
bin_peak <- function(x){
        m <- floor((as.numeric(x[2])+as.numeric(x[3]))/2)
        if(as.numeric(x[3])-as.numeric(x[2])>res){
                n <- ceiling((m-as.numeric(x[2])-(res/2))/res)
                tmp <- seq((m-(res/2)-res*n),(m+(res/2)+res*n),res)
		options(scipen = 999)
                f <- lapply(tmp[-length(tmp)],function(y){paste(x[1],y,y+res-1,paste("ENH_",x[4],sep=""),sep="\t")})
                return(f)
        }else{
		options(scipen = 999)
                return(paste(x[1],m-(res/2),m+((res/2)-1),paste("ENH_",x[4],sep=""),sep="\t"))}
}
l <- readLines("MAPS_peaks.broadPeak") ##consider changing this if you have changed the above
li <- strsplit(l,"\t")
a <- lapply(li,bin_peak)
cat(unlist(a),file="peaks.bed",sep="\n")
read.table("peaks.bed") ->b
b[,4] <- make.unique(b[,4],sep="")
#sprintf(paste0("%0", nchar(as.character(dim(b)[1])), "d"), 1:dim(b)[1])
write.table(b,file="peaks.bed",sep="\t",row.names=F,col.names=F,quote=F)

