library(gmodels)
library(RColorBrewer)
library(getopt)
library(ggplot2)


command=matrix(c( 
  "datafile" , "d" ,1, "character" ,
	"statfile" , "s" ,1, "character" ,
  "minsum"  , "m"  ,1, "integer" ,
  "weight"  , "w"  ,0, "logical",
  "help"    , "h" ,0,  "logical" 
                  ),byrow=T,ncol=4)
args=getopt(command)
if  ( ! is.null(args$help) || is.null(args$datafile) || is.null(args$statfile) || is.null(args$minsum)) {
     cat(paste(getopt(command, usage = T),  "\n" ))
     q()
} 

inname = args$datafile 
year<-sub(".*mCRASH-","",inname)
year<-sub("-ortho_number_Gene_count","",year)
minsum = args$minsum
statname = args$statfile
statwname = paste(statname,"w",sep="")

outname <- paste(inname,".pdf",sep="")
table <-  paste(inname,".xls",sep="")
titletext1 <- paste("PCA cluster via eliminateing HGT candidates within ",year," million years",sep="")
titletext2<- "Ratios of HGTs within "
titletext3<- paste(year," million yeas",sep="")

weighttable <-  paste(inname,".weight.xls",sep="")
expr <- read.table(inname, header=T, row.names=1)  

if(is.null(args$weight)){
expr[expr>0]<-1
numcol<-ncol(expr)-1
rowSums(expr[,1:numcol])->expr$Total.
subset(expr,expr$Total. > minsum)->expv
expv<-expv[,!grepl("Total.",colnames(expv))]
data <- t(expv)
datafas <- round(data) 

write.table(datafas,file=table,sep="",col.names= FALSE)
}else{

numcol<-ncol(expr)-1
rowSums(expr[,1:numcol])->expr$Total.
subset(expr,expr$Total. > minsum)->expv
expv<-expv[,!grepl("Total.",colnames(expv))]
data <- t(log(expv+1))
datafas <- round(data) 
datafas[datafas>9]<-9
write.table(datafas,file=weighttable,sep="",col.names= FALSE)
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
as.data.frame(data.pca$x)->data.fan
as.factor(rownames(data.pca$x))->data.fan$name
sub("[.]+.*","",rownames(data.fan))->data.fan$term
factor(data.fan$term,levels=unique(data.fan$term))

rm(data.label)

mycolors <- c(	"black",
				"grey",
				brewer.pal(9,"OrRd")[5],
				brewer.pal(8,"Dark2")[7],
				brewer.pal(8,"Set1")[2],
				brewer.pal(11,"PRGn")[2],
				brewer.pal(9,"Set1")[1],
				brewer.pal(9,"Greens")[8]
				)
sum(data.fan$PC2[which(data.fan$term=="Viridiplantae")])->sumpc2
sum(data.fan$PC1[which(data.fan$term=="Protozoa")])->sumpc1
if(sumpc1<0)data.fan$PC1<-data.fan$PC1*(-1)
if(sumpc2<0)data.fan$PC2<-data.fan$PC2*(-1)

labels<-c(
  "Protozoa..Deuterostomia..Homo_sapiens",
  "Rhodophyta..Rhodelphis..Rhodema",
  "Chromalveolate..Stramenopiles..Phatr",
  "Viridiplantae..Chlamydomonadaceae..Chlamre",
  "Viridiplantae..Tracheophyta..Zeama",
  "Viridiplantae..Characeae..Charabr",
  "Viridiplantae..Bathycoccacea..Ostresp",
  "Protozoa..Salpingoeca..Salpingoeca_rosetta",
  "Glaucophyta..Glaucocystophyceae..Cyanopa",
  "Chromalveolate..Cryptophyta..Guith",
  "Protozoa..Deuterostomia..Oncorhynchus_tshawytscha",
  "Viridiplantae..Tracheophyta..Arachhy",
  "Rhodophyta..Madagascaria..Madager",
  "Chromalveolate--Haptophyceae--Emihu"
  )
as.data.frame(data.fan[which(data.fan$name== "Protozoa..Deuterostomia..Mus_musculus"),]) -> data.label
for (i in labels) {
  print(i)
  rbind(data.label,as.data.frame(data.fan[which(data.fan$name==i),])) -> data.label
  
}

p<-ggplot(data=data.fan,aes(x=PC1,y=PC2,group=term,color=term,size=1))+geom_point(alpha=0.9,size=0.5)+
stat_ellipse(level = 0.8,show.legend = T,size=0.3,alpha=0.5)+
scale_color_manual(values=mycolors,name="Taxon",breaks=c("Archaea","Bacteria","Chromalveolate","Fungi","Glaucophyta","Protozoa","Rhodophyta","Viridiplantae"),
                         labels=c("Archaea","Bacteria","Chromalveolate","Fungi","Glaucophyta","Opisthokonta","Rhodophyta","Viridiplantae"))+
theme_bw()+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
  axis.line= element_line(colour = "black",size = 0.5),
  plot.title=element_text(hjust=0.5,colour="red"),
  plot.margin = unit(c(0.5,0.3,0.5,2.7),"cm"),
  axis.text = element_text(size = 8, colour = "black"),
  rect = element_rect(colour = "black", size = 0.6, linetype = 1),
  legend.text = element_text(size = 6,colour = "black"), 
  legend.title = element_text(size = 8, face = "bold", hjust = 0),
  legend.key.size = unit(1, "lines"),
  legend.key = element_rect(size=2),
  )+
#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"),plot.title=element_text(hjust=0.5,colour="red"))+
#annotate('text', label = 'setosa', x = -0.5, y = -0.25, size = 5, colour = '#f8766d') +
labs(x=xlab,y=ylab,title=titletext1)+
#stat_ellipse(level = 0.95)+
geom_text(data=data.label,aes(x=PC1,y=PC2,colour=term,label=name),size=1.5,check_overlap = TRUE,hjust = 0.5, nudge_x = -0.1,vjust = 1.3, nudge_y = -0.025)


##############################
#make stat histgraph file
####################
x<-read.table(statname,header=T,sep="\t")
stat<-as.data.frame(x)
library(RColorBrewer)
library(ggplot2)
colors<-c(
				"black",
				"grey",
				brewer.pal(8,"Dark2")[7],
				brewer.pal(11,"PRGn")[2],
				brewer.pal(9,"OrRd")[5],
				brewer.pal(8,"Set1")[2],
				brewer.pal(9,"Set1")[1],
				brewer.pal(9,"Greens")[8],
				brewer.pal(11,"Paired")[9],
				brewer.pal(11,"Paired")[1]
	)


stat$Donor<-factor(stat$Donor,levels=c("Archaea","Bacteria","Fungi","Metazoa","CRASH","Glaucophyta","Rhodophyta","Viridiplantae","other-Eukaryotes","unHGT"))
stat$kingdom<-factor(stat$kingdom,levels=c("Glaucophyta","Viridiplantae","Rhodophyta","Chromalveolate","Protozoa","Fungi","Bacteria","Archaea"),ordered=TRUE)
p1<-ggplot(data=stat,aes(x=factor(kingdom),y=num,fill=Donor))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+scale_fill_manual(values=colors)+guides(fill=FALSE)+
labs(x="Reciptor kingdoms",y="Ratio based on gene counts",title=titletext2,color="white")+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	axis.line= element_line(colour = "black",size = 0.5),
	axis.text = element_text(size = 8, colour = "black"),
	rect = element_rect(colour = "black", size = 0.3, linetype = 1),
	plot.title=element_text(hjust=1.0,colour="red"),
	
	plot.margin = unit(c(0,0,0.5,1),"cm"))+
coord_flip()
#######################################################
#make statw histgraph
#######################################################
xw<-read.table(statwname,header=T,sep="\t")
statw<-as.data.frame(xw)
statw$Donor<-factor(statw$Donor,levels=c("Archaea","Bacteria","Fungi","Metazoa","CRASH","Glaucophyta","Rhodophyta","Viridiplantae","other-Eukaryotes","unHGT"))
statw$kingdom<-factor(statw$kingdom,levels=c("Glaucophyta","Viridiplantae","Rhodophyta","Chromalveolate","Protozoa","Fungi","Bacteria","Archaea"),ordered=TRUE)
p2<-ggplot(data=statw,aes(x=factor(kingdom),y=num,fill=Donor))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+scale_fill_manual(values=colors)+
labs(x="",y="Ratio based on OGs counts",title=titletext3,color="red")+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	axis.line= element_line(colour = "black",size = 0.5),
	axis.text.y = element_blank(),
	
	rect = element_rect(colour = "black", size = 0.2, linetype = 1),
	legend.text = element_text(size = 6,colour = "black"), 
	legend.title = element_text(size = 8, face = "bold", hjust = 0),
	legend.key.size = unit(0.6, "lines"),
	plot.title=element_text(hjust=-0.05,colour="red"),
	plot.margin = unit(c(0,0,0.5,0),"cm"))+

coord_flip()

####################################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
####################################################
layout <- matrix(c(1, 1,1, 1, 1, 1, rep(1,6),rep(2, 3),rep(3,3)), nrow = 3, byrow = TRUE)


pdf(outname,width=8,height=6)
multiplot(p,p1,p2,layout=layout)

dev.off()


