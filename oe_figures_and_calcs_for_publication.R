#file is:"C:/Users/andre/Downloads/oe_signalptable.txt"
library(data.table)



###let's do this for zhan unmap

oe<-read.table("C:/Users/andre/Documents/Documents/Documents/OE/OE_unmap_Zhan.histo")
#plot the distribution of kmer freqs
plot(oe[50:680,],type="l", main="Kmer distribution for unmapped reads",xlab="Kmer coverage",ylab="Frequency")
#find the maximum for calcs below, I guessed and checked
abline(v=533)

#the estimated genome size is the product of the kmer frequencies and their abundances divided by the max coverage
#note we don't want to include error kmers, so we start at a coverage cutoff of 50
#I divide by 1 million to get the size in megabases 
(sum(as.numeric(oe[50:6888,1]*oe[50:6888,2]))/533)/1000000
#8.108162



apmb<-c(130.2, 61.6, 61.5, 22.9, 14, 9.2, 9, 8.806, 8.887)
apgen<-c(7093, 8789, 7540, 5268, 6375, 3870, 4002, 5270, 2712)

plot(apmb, apgen, xlab = "Genome size (Mb)", ylab = "Number of annotated genes", pch=16, ylim=c(2000, 10000), xlim=c(5, 130),
     main = "Comparative Apicomplexan genomics", cex = 1.44,cex.lab=1.4,cex.main=1.5)
#text( 12, 2700, "OE")




#install.packages('VennDiagram')
library(VennDiagram)



####################retry with BUSCOs because they're better annotated
#choose the appropriate full.tsv
ob<-fread("C:/Users/andre/Documents/Documents/Documents/OE/oe_redux_full.tsv",header=T,stringsAsFactors = F,fill=T)
lb<-fread("C:/Users/andre/Documents/Documents/Documents/OE/OE_L_full_redux.tsv",header=T,stringsAsFactors = F,fill=T)
ab<-fread("C:/Users/andre/Documents/Documents/Documents/OE/PgA_full.tsv",header=T,stringsAsFactors = F,fill=T)
bb<-fread("C:/Users/andre/Documents/Documents/Documents/OE/PgB_full.tsv",header=T,stringsAsFactors = F,fill=T)
gb<-fread("C:/Users/andre/Documents/Documents/Documents/OE/Gn_full.tsv",header=T,stringsAsFactors = F,fill=T)

o2<-ob[,c(1:2,7)]
l2<-lb[,c(1:2,7)]
a2<-ab[,c(1:2,7)]
b2<-bb[,c(1:2,7)]
g2<-gb[,c(1:2,7)]

ol<-merge(o2,l2,by="Buscoid")
ol<-ol[!duplicated(ol$Buscoid),]
olg<-merge(ol,g2,by="Buscoid")
olg<-olg[!duplicated(olg$Buscoid),]
ab<-merge(a2,b2,by="Buscoid")
ab<-ab[!duplicated(ab$Buscoid),]

busco<-merge(olg,ab,by="Buscoid")
busco<-busco[!duplicated(busco$Buscoid),]

#now we can get to work
#how many BUSCOs are conserved in all?
nrow(busco[which(busco$OE!="Missing"&busco$OElike!="Missing"&busco$Gn!="Missing"&busco$PgA!="Missing"&busco$PgB!="Missing"),])
#221 are found in all 5 ....that's just under 50%
nrow(busco[which(busco$OE=="Missing"&busco$OElike=="Missing"&busco$Gn=="Missing"&busco$PgA=="Missing"&busco$PgB=="Missing"),])
#only 25 are missing from all 5 tho, 5.6%

#first get sizes
#
table(busco$OE)
#366
table(busco$OElike)
#268

nrow(busco[which(busco$OE!="Missing"&busco$OElike!="Missing"),])
#264

#let's start with OE vs OE-like
#plot(-4,-4,ylim=c(-1,1),xlim=c(1,1),xaxt='n',yaxt='n',xlab='',ylab='')
dev.new(width=10, height=10)
draw.pairwise.venn(366, 268, 264,fill=c("lightsalmon3","tan1"),alpha=0.3,cex=1.8)
#text(-0.01,0.4,"OE",cex=1.5)


#let's do a quad vd then


table(busco$PgA)
#356
table(busco$PgB)
#355
table(busco$Gn)
#397

#how many BUSCOs are conserved in 2?
nrow(busco[which(busco$OE!="Missing"&busco$PgA!="Missing"),])
nrow(busco[which(busco$OE!="Missing"&busco$PgB!="Missing"),])
nrow(busco[which(busco$OE!="Missing"&busco$Gn!="Missing"),])

nrow(busco[which(busco$PgA!="Missing"&busco$PgB!="Missing"),])
nrow(busco[which(busco$PgA!="Missing"&busco$Gn!="Missing"),])
nrow(busco[which(busco$PgB!="Missing"&busco$Gn!="Missing"),])


#how many in 3?
nrow(busco[which(busco$OE!="Missing"&busco$PgA!="Missing"&busco$PgB!="Missing"),])
nrow(busco[which(busco$OE!="Missing"&busco$Gn!="Missing"&busco$PgA!="Missing"),])
nrow(busco[which(busco$OE!="Missing"&busco$Gn!="Missing"&busco$PgB!="Missing"),])
nrow(busco[which(busco$Gn!="Missing"&busco$PgA!="Missing"&busco$PgB!="Missing"),])


#in all 4?
nrow(busco[which(busco$OE!="Missing"&busco$Gn!="Missing"&busco$PgA!="Missing"&busco$PgB!="Missing"),])
#297

#OE 1, PgA 2, Pg B 3, Gn 4
plot.new()
draw.quad.venn(366, 356, 355, 397, 320, 315, 348, 343, 339,
               339, 307, 309, 305, 327, 297)


#can we simplify the Porosporas?
nrow(busco[which(busco$PgA!="Missing"&busco$PgB!="Missing"),])
#343
#so at best we have 13 and 12 unique BUSCOs

#so let's just take PgA, the bigger one

#OE 1, PgA 2, Gn 3
#plot.new()
dev.new(width=5, height=5)
draw.triple.venn(366,356,397,320,339,348,309,fill=c("lightsalmon3","firebrick3","lemonchiffon3"),alpha=0.3,cex=1.8)


nrow(busco[which(busco$OE=="Missing"&busco$Gn=="Missing"&busco$PgA=="Missing"),])
#25 missing from all 3
#what's unique to OE?

####now functional classes############
#this is our core conserved set:
trucon<-table(busco[which(busco$OE!="Missing"&busco$Gn!="Missing"&busco$PgA!="Missing"&busco$PgB!="Missing"),3])
#these are missing in OE only:
oelack<-table(busco[which(busco$Gn!="Missing"&busco$PgA!="Missing"&busco$OE=="Missing"),7])


#unique "conserved" orthologs
#OE
oep<-table(busco[which(busco$Gn=="Missing"&busco$PgA=="Missing"&busco$OE!="Missing"),3])

#Pg
pgp<-table(busco[which(busco$Gn=="Missing"&busco$PgA!="Missing"&busco$OE=="Missing"),8])

#Gn
gnp<-table(busco[which(busco$Gn!="Missing"&busco$PgA=="Missing"&busco$OE=="Missing"),7])

nrow(busco[which(busco$OE=="Missing"&busco$PgA!="Missing"),])
nrow(busco[which(busco$OE!="Missing"&busco$Gn!="Missing"&busco$PgA!="Missing"),])
nrow(busco[which(busco$OE!="Missing"&busco$Gn!="Missing"&busco$PgB!="Missing"),])


###Figure 1####
megs<-c(130.2,61.6,61.5,22.9,14,9.2,9.0,8.86,8.81)
gens<-c(7093,8789,7540,5268,6375,3870,4002,2633,5270)
par(mai=c(1,1.2,1,0.1))
plot(megs,gens,xlab="",ylab="",las=1,pch=16,main="Comparative Apicomplexan genomics",cex=1.8,cex.axis=1.4,cex.main=2)
mtext(side=2,line=4.3,cex=1.9,"Number of Annotated Genes")
mtext(side=1,line=3.0,cex=1.9,"Genome Size (Mb)")




#Figure 4########
library("phytools")
library("phangorn")
#atp<-read.newick("C:/Users/andre/Documents/Documents/Documents/OE/Lehane_seqs+Apicomplex_ML_boot_for_plotting.nwk",format="raxml")
atp<-read.newick("C:/Users/andre/Documents/Documents/Documents/OE/Lehane_seqs+Apicomplex_ML_boot_for_plotting.nwk")
#par(xpd = TRUE)


#x.lim=c(-0.05,6),y.lim=c(-0.05,75)
pdf("C:/Users/andre/Documents/Documents/Documents/OE/Figure4.pdf",width=11,height=8,paper="a4r")
plotBS(atp,type="phylogram",method="FBP",cex=0.6)
add.scale.bar(x=-0.1,y=-1.1)
#PMCAs
segments(x0=0.2,x1=6.6, y0=0,y1=0,lwd=2)
segments(x0=0.2,x1=3.5, y0=18,y1=18,lwd=2)
segments(x0=3.5,x1=6.6, y0=19,y1=19,lwd=2)
segments(x0=3.5,x1=3.5, y0=18,y1=19,lwd=2)
segments(x0=0.2,x1=0.2, y0=0,y1=18,lwd=2)
segments(x0=0.2,x1=0.2, y0=18,y1=18,lwd=2)
segments(x0=6.6,x1=6.6, y0=0,y1=19,lwd=2)
segments(x0=6,x1=6, y0=0,y1=19,lwd=2)
rect(xleft=5.1,ybottom=0,xright=6,ytop=19, col = "dodgerblue2")
text(x=5.48,y=10.5,"PMCA ATPases",cex=0.9)
text(x=5.48,y=8,expression(Ca^+2),cex=0.9)
#Type IIs
segments(x0=0.2,x1=6.6, y0=26.45,y1=26.45,lwd=2)
segments(x0=0.2,x1=0.2, y0=18,y1=26.45,lwd=2)
rect(xleft=5.1,ybottom=19,xright=6,ytop=26.45, col = "palevioletred1")
text(x=5.48,y=23.225,"Type II ATPases",cex=0.9)
text(x=5.48,y=21.2,expression(paste(Na^"+"," / K"^"+" , " , ", H^"+"," / ", K^"+")),cex=0.9)
#SERCAs
segments(x0=0.09,x1=6.6, y0=26.45,y1=26.45,lwd=2)
segments(x0=0.09,x1=6.6, y0=52.5,y1=52.5,lwd=2)
segments(x0=0.09,x1=0.09,y0=26.45,y1=52.5,lwd=2)
rect(xleft=5.1,ybottom=26.45,xright=6.6,ytop=52.5, col = "khaki1")
text(x=5.48,y=41.5,"SERCA ATPases",cex=0.9)
text(x=5.45,y=39,expression(Ca^+2),cex=0.9)
#ATP4s
segments(x0=0.09,x1=6.6, y0=67.5,y1=67.5,lwd=2)
segments(x0=0.09,x1=0.09,y0=52.5,y1=67.5,lwd=2)
rect(xleft=5.1,ybottom=52.5,xright=6.6,ytop=67.5, col = "sienna1")
text(x=5.45,y=62.5,"ATP4 ATPases",cex=0.9)
text(x=5.45,y=59,expression(Na^"+"),cex=0.9)
#ENAs
segments(x0=0.09,x1=6.6, y0=73.5,y1=73.5,lwd=2)
segments(x0=0.09,x1=0.09,y0=67.5,y1=73.5,lwd=2)
rect(xleft=5.1,ybottom=67.5,xright=6.6,ytop=73.5, col = "aquamarine")
text(x=5.48,y=71.5,"ENA ATPases",cex=0.9)
text(x=5.48,y=69,expression(paste(Na^"+"," / K"^"+")),cex=0.9)
points(c(3.9,3.9,5.05,5.05),c(15,16,17,18), pch=17, col = "black",cex=1.3)
points(c(1.9,1.9),c(38,39), pch=17, col = "black",cex=1.3)
points(c(2.62,2.62),c(8,9),pch=15, col = "black", cex=1.3)
points(c(1.7,1.7,1.7,2),c(27,28,29,33), pch=15, col = "black",cex=1.3)
points(c(2.3,1.9,1.5,1.5),c(63,64,65,66),pch=15,col="black",cex=1.3)
legend(x=4,y=7, pch=c(15,17),c("Gregarine", "Ophryocystis"),text.font=c(1,3),cex=0.9)
dev.off()



#let's try some summary intron stats..
library(GenomicFeatures) 
#library(rtracklayer)

txdb <- makeTxDbFromGFF(file = "C:/Users/andre/Documents/Documents/Documents/OE/OE_for_sharing/Ophryocystis-elektroscirrha-v1.0.a2.638bc1454374c-genemark.genes.gff3", format = "gff") 
introns <- intronsByTranscript(txdb)

intcount<-elementNROWS(introns)
mean(intcount)

# mean, median all intron 
i <- unlist(introns)
 summary(width(i))


