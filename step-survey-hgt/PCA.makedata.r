library(gmodels)
library(RColorBrewer)
library(getopt)
library(ggplot2)


command=matrix(c( "datafile" , "d" ,1, "character" ,

                  "minsum" , "m" ,1, "integer" ,
				  "weight" , "w",0,"logical",
                  "help" , "h" ,0, "logical" 
                  ),byrow=T,ncol=4)
args=getopt(command)
if  ( ! is.null(args$help) || is.null(args$datafile) || is.null(args$minsum)) {
     cat(paste(getopt(command, usage = T),  "\n" ))
     q()
} 

inname = args$datafile 
minsum = args$minsum
year<-sub("a30_h90gene.mCRASH-","",inname)
year<-sub("-ortho_number_Gene_count","",year)
pcatable<- paste(inname,".pca",sep="")
expr <- read.table(inname, header=T, row.names=1)  

if(is.null(args$weight)){
expr[expr>0]<-1
numcol<-ncol(expr)-1
rowSums(expr[,1:numcol])->expr$Total.
subset(expr,expr$Total. > minsum)->expv
expv<-expv[,!grepl("Total.",colnames(expv))]
data <- t(expv)



}else{

numcol<-ncol(expr)-1
rowSums(expr[,1:numcol])->expr$Total.
subset(expr,expr$Total. > minsum)->expv
expv<-expv[,!grepl("Total.",colnames(expv))]
data <- t(log(expv+1))


}

data.pca <- fast.prcomp(data) 

a <- summary(data.pca)   
tmp <- a[4]$importance 
pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100
xlab <- paste("PC1","(",pro1,"%)",sep="")
ylab <- paste("PC2","(",pro2,"%)",sep="")
xmax <- max(data.pca$x[,1])  
xmin <- min(data.pca$x[,1])
ymax <- max(data.pca$x[,2])
ymin <- min(data.pca$x[,2])
   
  
samples =rownames(data.pca$x)
as.data.frame(data.pca$x[,1:3])->data.fan
as.factor(rownames(data.pca$x))->data.fan$name
c(rep(year,nrow(data.fan)))->data.fan$year
sub("[.]+.*","",rownames(data.fan))->data.fan$term

write.table(data.fan,file=pcatable,sep="\t")
