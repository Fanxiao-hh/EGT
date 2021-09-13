library(gmodels)
library(RColorBrewer)
library(getopt)
library(ggplot2)


command=matrix(c( "datafile" , "d" ,1, "character" ,
	              "statfile", "s" ,1, "character" ,			  
                  "help" , "h" ,0, "logical" 
                  ),byrow=T,ncol=4)
args=getopt(command)
if  ( ! is.null(args$help) || is.null(args$datafile) || is.null(args$statfile) ) {
     cat(paste(getopt(command, usage = T),  "\n" ))
     q()
} 

inname = args$datafile 
year<-sub("a30_h90gene.mCRASH-","",inname)
year<-sub("-ortho_number_Gene_count.pca","",year)
statname = args$statfile
statwname = paste(statname,"w",sep="")

outname <- paste(inname,".pdf",sep="")
table <-  paste(inname,".xls",sep="")
titletext1 <- paste("PCA cluster via eliminateing HGT candidates within ",year," million years",sep="")
titletext2<- "Ratios of HGTs within "
titletext3<- paste(year," million yeas",sep="")


expr <- read.table(inname, header=T, row.names=1,sep="\t")   
as.data.frame(expr)->data.fan

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
  "Viridiplantae..Tracheophyta..Jatrocu",
  "Viridiplantae..Tracheophyta..Heliaan",
  "Viridiplantae..Tracheophyta..Duriozi")
as.data.frame(data.fan[which(data.fan$name=="Viridiplantae..Tracheophyta..Spinaol"),]) -> data.label
for (i in labels) {
  print(i)
  rbind(data.label,as.data.frame(data.fan[which(data.fan$name==i),])) -> data.label
  
}
data.label

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
geom_text(data=data.label,aes(x=PC1,y=PC2,colour=term,label=name),size=2,check_overlap = TRUE,hjust = 0, nudge_x = -0.1,vjust = 1, nudge_y = -0.025)


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
stat$kingdom<-factor(stat$kingdom,levels=c("Viridiplantae","Rhodophyta","Glaucophyta","Chromalveolate","Protozoa","Fungi","Bacteria","Archaea"),ordered=TRUE)
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


