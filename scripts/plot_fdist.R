#Author: Sarah P. Flanagan
#Date: 18 September 2015
#Revised: 22 April 2016
#Purpose: Analyze data from FDIST2 and numerical analysis and plot figures
#to go along with the fst-heterozygosity paper.

rm(list=ls())

###############################################################################
#DATA FROM LITERATURE
###############################################################################
setwd("E://ubuntushare//fst_outliers//results//data_from_literature")
###################REFORMATTING
smith.dat<-read.delim("smith2015_skewed.txt")
locus.id.index<-c(3,5,7,9,11,13,15,17,19,21)
loci<-colnames(smith.dat)[locus.id.index]

new.dat<-data.frame(matrix(NA,ncol=length(loci),nrow=nrow(smith.dat)))
for(j in 1:length(loci)){
	dat<-smith.dat[,c(locus.id.index[j],(locus.id.index[j]+1))]
	alleles<-levels(as.factor(c(dat[,1],dat[,2])))
	for(i in 1:length(alleles)){
		dat[,1]<- replace(dat[,1],dat[,1]==alleles[i],paste(0,i,sep=""))
		dat[,2]<- replace(dat[,2],dat[,2]==alleles[i],paste(0,i,sep=""))
	}
	new.dat[,j]<-paste(dat[,1],dat[,2],sep="")
}
new.dat<-cbind(smith.dat$Locality,smith.dat[,1],new.dat)
dat.split<-split(new.dat,new.dat[,1])

write.table(loci,"smith.genepop",quote=F,col.names=F,row.names=F)
for(i in 1:length(dat.split)){
	write.table("POP","smith.genepop",quote=F,append=T,
		col.names=F,row.names=F)
	dat<-dat.split[[i]]
	dat$int<-","
	dat<-cbind(dat[,2],dat$int,dat[,3:12])
	write.table(dat,"smith.genepop",quote=F,col.names=F,row.names=F,
		sep='\t',append=T)
}

bourret.dat<-read.delim("Bourret_Genotypes_PowerMarker.txt")
convert.snps<-function(one.col){
	new.dat<-gsub("(\\d)/(\\d)","0\\10\\2",one.col)
	new.dat<-replace(new.dat,new.dat=="0001","0102")
	new.dat<-replace(new.dat,new.dat=="0100","0201")
	new.dat<-replace(new.dat,new.dat=="0101","0202")
	new.dat<-replace(new.dat,new.dat=="0000","0101")
	new.dat<-replace(new.dat,is.na(new.dat),"0000")
	return(new.dat)
}
new.bourret<-apply(bourret.dat[,4:ncol(bourret.dat)],2,convert.snps)
loci.names<-colnames(new.bourret)
new.bourret<-data.frame(cbind(as.character(bourret.dat$Level.2),
	as.character(bourret.dat$Level.1),new.bourret))
bour.split<-split(new.bourret,new.bourret[,1])

write.table(loci.names,"bourret.genepop",quote=F,col.names=F,row.names=F)
for(i in 1:length(bour.split)){
	write.table("POP","bourret.genepop",quote=F,append=T,
		col.names=F,row.names=F)
	dat<-bour.split[[i]]
	dat$int<-","
	dat<-cbind(dat[,2],dat$int,dat[,3:(ncol(dat)-1)])
	write.table(dat,"bourret.genepop",quote=F,col.names=F,row.names=F,
		sep='\t',append=T)
}

hebert.dat<-read.delim("hebert2013_snps.txt")
heb.new<-data.frame(matrix(NA,nrow=nrow(hebert.dat),ncol=ncol(hebert.dat)))
heb.new[,1]<-hebert.dat[,1]
for(i in 2:ncol(hebert.dat)){
	gt<-do.call("rbind",strsplit(as.character(hebert.dat[,i]),"/"))
	alleles<-levels(as.factor(gt))
	alleles<-alleles[!(alleles == "0")]
	for(j in 1:length(alleles)){
		gt<-replace(gt,gt==alleles[j],paste(0,j,sep=""))
	}
	gt<-replace(gt,gt=="0","00")
	heb.new[,i]<-paste(gt[,1],gt[,2],sep="")
}
loci.names<-colnames(hebert.dat)[-1]
heb.new$Pop<-gsub("\\w-(\\w)\\d+","\\1",heb.new[,1])
heb.split<-split(heb.new,heb.new$Pop)

write.table(loci.names,"hebert.genepop",quote=F,col.names=F,row.names=F)
for(i in 1:length(heb.split)){
	write.table("POP","hebert.genepop",quote=F,append=T,
		col.names=F,row.names=F)
	dat<-heb.split[[i]]
	dat$int<-","
	dat<-cbind(dat[,1],dat$int,dat[,2:(ncol(dat)-2)])
	write.table(dat,"hebert.genepop",quote=F,col.names=F,row.names=F,
		sep='\t',append=T)
}

moura.dat<-read.delim("Moura_2014.vcf",comment.char="#")
header.start<-grep("#CHROM",scan("Moura_2014.vcf",what="character"))
header1<-scan("Moura_2014.vcf",what="character")[header.start:
	(header.start+ncol(moura.dat)-1)]
colnames(moura.dat)<-header1
groups<-c(rep("AR",17),rep("SR",13),rep("RU",9),rep("BS",13),rep("AT",21),
	rep("CT",16),rep("OS",7),rep("IC",6),rep("MI",13))

####################PLOTTING
normal.loci<-read.delim("Hess_2013_data_Genepop.genepop.loci")#Hess et al. 2013
normal.ci<-read.delim("Hess_2013_data_Genepop.genepop.ci")
#incline.loci<- read.delim("TrieralpsmaleNEWnosingletons.genepop.loci")#Trier et al. 2014
#incline.loci<-incline.loci[incline.loci$Fst>-1,]
#incline.ci<- read.delim("TrieralpsmaleNEWnosingletons.genepop.ci")
jaggedci.loci<-read.delim("hebert.genepop.loci")#Hebert et al. 2013
jaggedci.ci<-read.delim("hebert.genepop.ci")
#skewed.loci <-read.delim("smith.genepop.loci")#Smith et al. 2015
#skewed.ci <-read.delim("smith.genepop.ci")
skewed<-read.delim("dann.2012.traced.txt")

pdf("../Fig1_literature.pdf",height=7,width=10)
png("../Fig1_literature.png",height=10,width=7,units="in",res=300)
par(mfrow=c(3,1),oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(normal.loci$Het, normal.loci$Fst,xlab="",ylab="",pch=19,las=1)
points(normal.ci$Het,normal.ci[,2],col="red",type="l",lwd=2)
points(normal.ci$Het, normal.ci[,4],col="red",type="l",lwd=2)
text(x=0.115,y=0.55,"A. Well-behaved")

#plot(incline.loci$Het, incline.loci$Fst,xlab="",ylab="",pch=19,las=1)
#points(incline.ci$Het,incline.ci[,2],col="red",type="l",lwd=2)
#points(incline.ci$Het, incline.ci[,4],col="red",type="l",lwd=2)
#text(x=0.05,y=0.21,"B. Incline")

plot(jaggedci.loci$Het, jaggedci.loci$Fst,xlab="",	ylab="",pch=19,las=1)
points(jaggedci.ci$Het,jaggedci.ci[,2],col="red",type="l",lwd=2)
points(jaggedci.ci$Het, jaggedci.ci[,4],col="red",type="l",lwd=2)
text(x=0.25,y=0.85,"B. Incline, Jagged CI")

plot(skewed$points.x, skewed$points.y,xlab="",	ylab="",pch=19,las=1)
points(skewed$lower.x,skewed$lower.y,col="red",type="l",lwd=2)
points(skewed$upper.x, skewed$upper.y,col="red",type="l",lwd=2)
text(x=0.095,y=0.27,"C. Skewed")

mtext(expression(italic(H)[italic(T)]),1,outer=T,cex=0.85)
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.85)

dev.off()

###############################################################################
#FDIST2
###############################################################################
#PLOT ALL THE PARAMETER COMBINATIONS
setwd("E://ubuntushare//fst_outliers//results//fdist2")

file.list<-list.files(pattern="out.dat[2-5]\\d{2}.txt")
for(i in 1:length(file.list)){
	dat<-read.delim(file.list[i], header=F,sep=' ')
	if(dat[1,1]!="NaN"){
	plot_name<-paste(sub('(out.dat\\d+).txt','\\1',file.list[i]),".png",sep="")
	png(plot_name)
	plot(dat$V1,dat$V2,xlab="Exp Het",ylab="Fst",pch=19)
	dev.off()
	}
}

#####MAKE FIG 2
tl<-read.delim("out.dat21.txt",header=F,sep=' ')
bl<-read.delim("out.dat6.txt",header=F,sep=' ')
tr<-read.delim("out.dat171.txt",header=F,sep=' ')
br<-read.delim("out.dat156.txt",header=F,sep=' ')

png("../Figure2_FDISTresults.png",height=7,width=7,units="in",res=300)
pdf("../Figure2_FDISTresults.pdf")
par(mfrow=c(2,2),oma=c(2,6,2,2),mar=c(2,2,2,2))
plot(tl$V1, tl$V2, pch=19, ylab="",xlab="",las=1)
mtext(expression(Estimated~italic(F)[ST]:~0.2),outer=F)
mtext("Number of\nSampled\nDemes:",2,outer=F,line=3,las=1,at=0.6)
mtext("10",2,outer=F,line=4,las=1)
plot(tr$V1, tr$V2, pch=19, ylab="",xlab="",las=1)
mtext(expression(Estimated~italic(F)[ST]:~0.8),outer=F)
plot(bl$V1, bl$V2, pch=19, ylab="",xlab="",las=1)
mtext(expression(italic(F)[ST]),2,outer=T)
mtext("75",2,outer=F,line=4,las=1)
plot(br$V1, br$V2, pch=19, ylab="",xlab="",las=1)
mtext("Expected Heterozygosity",1,outer=T)
dev.off()
###############################################################################
################NUMERICAL ANALYSIS
###############################################################################
#***EXPLORATORY ANALYSIS***#
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis_original_runs")
num<-read.delim( "n10000_s100_m0.01_r1000_0.5freqs.txt")
png("fst.time.png", res=300, width=174, height=87, units="mm")
par(mfrow=c(1,2), oma=c(2,2,1,1),mar=c(2,4,1,1),cex=0.5)
plot(num$Time, num$WrightsFst, xlab="", ylab="Fst",las=1,type="l")
points(num$Time,num$WCFst,type="l",col="blue")
points(num$Time, num$WCFstSimple, type="l", col="red")
legend("bottomright",c("Wright's Fst", "Weir&Cockerham Fst"),
	col=c("black", "blue"), lty=1)
plot(num$Time, num$AvgP, xlab="",ylab="Frequency", 
	type="l",las=1,ylim=c(0,1))
points(num$Time, num$AvgH,type="l",col="purple")
legend("bottomright",c("Avg p", "Ht"),
	col=c("black", "purple"), lty=1)
mtext("Time (Generations)", 1, outer=T,cex=0.5)
dev.off()
#select all rows with pop Hs
pop.het<-num[,grep("Pop\\d+Het",colnames(num))]
hs<-rowMeans(pop.het)
ht<-1-(0.5*0.5+0.5*0.5)
fst<-(ht-hs)/ht
#fst-het
dat<-read.delim("n10000_s100_m0.01_r1000_output.txt")
plot(dat$Ht, dat$WrightsFst)
#check equilibrium for all parameter combos
setwd("E://ubuntushare//fdist2//numerical_analysis//")
out.files<-list.files(pattern="output")
base.names<-unlist(strsplit(out.files,"output.txt"))
colors<-rainbow(18)
for(i in 1:length(base.names)){
	freq.files<-list.files(pattern=base.names[i])
	freq.files<-freq.files[grep("freq",freq.files)]
	png(paste(base.names[i],"eq.png",sep=""), 
		res=300, width=87, height=87, units="mm")
	par(oma=c(2,2,1,1),mar=c(2,4,1,1),cex=0.5)
	for(j in 1:length(freq.files)){
		num<-read.delim(freq.files[j])
		if(j == 1){	
			plot(num$Time, num$WrightsFst, xlab="Generations", 
				ylab="Fst",las=1,type="l",col=colors[j])
		}
		else {
			points(num$Time,num$WrightsFst,type="l",col=colors[j])
		}
	}
	legend("bottomright",c("0.1", "0.15", "0.2", "0.25", "0.3", "0.35", 
		"0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8",
		"0.85", "0.9", 0.95), col=colors, lty=1)
	dev.off()

}
#evaluate output
setwd("E://ubuntushare//fdist2//numerical_analysis")#//New folder//")
out.files<-list.files(pattern="output")
for(i in 1:length(out.files)){
	dat<-read.delim(out.files[i])
	base.name<-strsplit(out.files[i],".txt")
	png(paste(base.name,".png",sep=""), 
		res=300, width=87, height=87, units="mm")
	par(oma=c(2,2,1,1),mar=c(2,4,1,1),cex=0.5)
	plot(dat$Ht, dat$WrightsFst, xlab="Ht",ylab="Wright's Fst", pch=19)
	dev.off()
}
library(ggplot2)
setwd("E://GitHub//fst_outliers//numerical_analysis")
test<-list.files(pattern="n5000_s50_m0.01")
dat<-read.delim(test[19])
png("n5000_s50_m0.01.test.png",width=174, height=174, units="mm",res=300)
qplot(dat$Ht, dat$WrightsFst, colour=as.factor(dat$pbar))
dev.off()
par(oma=c(2,2,1,1),mar=c(2,4,1,1),cex=0.5)
	for(j in 1:(length(test)-1)){
		num<-read.delim(test[j])
		if(j == 1){	
			plot(num$Time, num$AvgP, xlab="Generations", 
				ylab="Fst",las=1,type="l",col=colors[j], ylim=c(0,1))
		}
		else {
			points(num$Time,num$AvgP,type="l",col=colors[j])
		}
	}
	legend("bottomright",c("0.1", "0.15", "0.2", "0.25", "0.3", "0.35", 
		"0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8",
		"0.85", "0.9", 0.95), col=colors, lty=1)
#plot actual vs. sampled het-fst
setwd("E://ubuntushare//fdist2//numerical_analysis//with_sampling//finished")
setwd("E://Docs//fst-het//results")
samp.files<-list.files(pattern="sampledpops.txt")
dat.files<-list.files(pattern="output.txt")
for(i in 1:length(samp.files)){
	samp<-read.delim(samp.files[i])
	dat<-read.delim(dat.files[i])
	base.name<-strsplit(samp.files[i],"sampledpops.txt")
	p<-strsplit(strsplit(samp.files[i], "p")[[1]][2],"[._]")[[1]][1]
	s<-strsplit(strsplit(samp.files[i], "s")[[1]][2],"[._]")[[1]][1]
	png(paste(base.name,"png",sep=""), 
		width=174, height=87, units="mm", res=300)
	par(mfrow=c(1,2), oma=c(2,2,1,1),mar=c(2,4,1,1),cex=0.5)
	plot(dat$Ht, dat$WrightsFst, xlab="", ylab="Wright's Fst", las=1,pch=19)
	legend("topleft", c(paste(p," populations", sep="")), bty="n")
	plot(samp$Ht, samp$WrightsFst, xlab="", ylab="", las=1,pch=19)
	legend("topleft", bty="n",
		c(paste(s, " Samples from each of ",p, " pops",sep="")))
	mtext("Ht", 1, outer=T, cex=0.5)
	dev.off()
}

############################################################################
#PLOT FIG 3-REVISIONS
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis")


nm.list<-list(
	read.delim("Nm0.1.d2.s2.fig.output.txt"),
	read.delim("Nm1.d2.s2.fig.output.txt"),
	read.delim("Nm10.d2.s2.fig.output.txt"),
	read.delim("Nm0.1.d5.s5.fig.output.txt"),
	read.delim("Nm1.d5.s5.fig.output.txt"),
	read.delim("Nm10.d5.s5.fig.output.txt"))
names(nm.list)<-c("Nm0.1.d2.s2","Nm1.d2.s2","Nm10.d2.s2",
	"Nm0.1.d5.s5","Nm1.d5.s5","Nm10.d5.s5")
Nms<-c(0.1,1,10,0.1,1,10)
nm.cis<-list(
	read.delim("Nm0.1.d2.s2.fig.genepop.ci"),
	read.delim("Nm1.d2.s2.fig.genepop.ci"),
	read.delim("Nm10.d2.s2.fig.genepop.ci"),
	read.delim("Nm0.1.d5.s5.fig.genepop.ci"),
	read.delim("Nm1.d5.s5.fig.genepop.ci"),
	read.delim("Nm10.d5.s5.fig.genepop.ci"))

png("../Fig3_NmDemes.png",width=169,height=169,units="mm",res=300)
pdf("../Fig3_NmDemes.pdf")
par(mfrow = c(2, 3),cex = 0.6,mar = c(0, 0, 0, 0), 
	oma = c(4, 4.5,1.5, 0.5), tcl = -0.25,mgp = c(2, 0.6, 0))
for (i in 1:length(nm.list)) {
	#plot(seq(0,0.6,0.06),seq(0,1,0.1), axes = FALSE, type = "n")
	y.min<-min(min(nm.cis[[i]][,2]),min(nm.list[[i]]$WrightsFst))
	y.max<-(max(nm.list[[i]]$WrightsFst)+2*sd(nm.list[[i]]$WrightsFst))
	if(y.max>1)
		y.max<-1
	plot(nm.list[[i]]$Ht,nm.list[[i]]$WrightsFst,las=1,ylab="",xlab="",
		col="black",pch=19,xaxt="n", xlim=c(0,0.6),ylim=c(y.min,y.max))
	points(nm.cis[[i]]$Het,nm.cis[[i]][,2],col="red",type="l")
	points(nm.cis[[i]]$Het, nm.cis[[i]][,4],col="red",type="l")
	Nm<-Nms[i]
	exp.fst<-round(1/((4*as.numeric(Nm))+1),digits=3)
	avg.fst<-round(mean(nm.list[[i]]$WrightsFst),digits=3)
	leg.nms<-c(bquote(Exp.~italic(F)[ST]~"="~.(exp.fst)),
		bquote(Mean~italic(F)[ST]~"="~.(avg.fst)))
	legend("topleft",legend=as.expression(leg.nms),bty="n")
	if (i %in% c(4,5,6))
		axis(1, at = seq(0,0.6,0.1))
	if(i==1){
		mtext("Nm = 0.1", 3,cex=.75)
		mtext("2 Demes",2,line=1.75,cex=.75)}
	if(i==2)
		mtext("Nm = 1",3,cex=.75)
	if(i==3){
		mtext("Nm=10",3,cex=.75)
	}
	if(i==4)
		mtext("5 Demes",2,line=1.75,cex=.75)
}
mtext(expression(italic(H)[italic(T)]),1,outer=T,line=2,cex=.75)
mtext(expression(Wright~"'"~s~italic(F)[ST]),2,outer=T,line=3,cex=.75)
dev.off()

############################################################################
#PLOT FIG 4-REVISIONS (FORMER FIG 5)
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis")

smp.list<-list(
	read.delim("Nm1.d2.s2.fig.sampledpops.txt"),
	read.delim("Nm1.d2.s4.fig.sampledpops.txt"),
	read.delim("Nm1.d2.s8.fig.sampledpops.txt"),
	read.delim("Nm1.d5.s5.fig.sampledpops.txt"),
	read.delim("Nm1.d5.s10.fig.sampledpops.txt"),
	read.delim("Nm1.d5.s20.figs.sampledpops.txt"))
names(smp.list)<-c("Nm1.d2.s2.fig","Nm1.d2.s4.fig","Nm1.d2.s8.fig",
	"Nm1.d5.s5.fig","Nm1.d5.s10.fig","Nm1.d5.s20.fig")
smp.cis<-list(
	read.delim("Nm1.d2.s2.fig.genepop.ci"),
	read.delim("Nm1.d2.s4.fig.genepop.ci"),
	read.delim("Nm1.d2.s8.fig.genepop.ci"),
	read.delim("Nm1.d5.s5.fig.genepop.ci"),
	read.delim("Nm1.d5.s10.fig.genepop.ci"),
	read.delim("Nm1.d5.s20.figs.genepop.ci"))

png("../Fig4_samplesize.png",width=169,height=169,units="mm",res=300)
pdf("../Fig4_samplesize.pdf")
par(mfrow = c(2, 3),cex = 0.6,mar = c(0, 0, 0, 0), 
	oma = c(4, 4.5, 1.5, 0.5), tcl = -0.25,mgp = c(2, 0.6, 0))
for (i in 1:length(smp.list)) {
	#plot(seq(0,0.6,0.06),seq(0,1,0.1), axes = FALSE, type = "n")
	y.min<-min(min(smp.cis[[i]][,2]),min(smp.list[[i]]$WrightsFst))
	y.max<-(max(smp.list[[i]]$WrightsFst)+2*sd(smp.list[[i]]$WrightsFst))
	if(y.max>1)
		y.max<-1
	plot(smp.list[[i]]$Ht,smp.list[[i]]$WrightsFst,las=1,xlab="",ylab="",
		col="black",pch=19,xaxt="n", xlim=c(0,0.6),ylim=c(0,y.max))
	points(smp.cis[[i]]$Het,smp.cis[[i]][,2],col="red",type="l")
	points(smp.cis[[i]]$Het, smp.cis[[i]][,4],col="red",type="l")
	Nm<-strsplit(strsplit(names(smp.list)[i], "Nm")[[1]][2],"[._]")[[1]][1]
	exp.fst<-round(1/((4*as.numeric(Nm))+1),digits=4)
	avg.fst<-round(mean(smp.list[[i]]$WrightsFst),digits=4)
	leg.nms<-c(bquote(Exp.~italic(F)[ST]~"="~.(exp.fst)),
		bquote(Mean~italic(F)[ST]~"="~.(avg.fst)))
	legend("topleft",legend=as.expression(leg.nms),bty="n")

	if (i %in% c(4,5,6))
		axis(1, at = seq(0,0.6,0.1))
	if(i==1){
		mtext("1 Sample Per Deme", 3,cex=.75)
		mtext("2 Demes",2,line=1.75,cex=.75)}
	if(i==2)
		mtext("2 Samples Per Deme",3,cex=.75)
	if(i==3){
		mtext("4 Samples Per Deme",3,cex=.75)
	}
	if(i==4)
		mtext("5 Demes",2,line=1.75,cex=.75)
}
mtext(expression(italic(H)[italic(T)]),1,outer=T,line=2,cex=.75)
mtext(expression(Wright~"'"~s~italic(F)[ST]),2,outer=T,line=3,cex=.75)

dev.off()

############################################################################
#PLOT FIG 5-REVISIONS (FORMERLY PART OF S1-4)
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis_genepop")

find.los.sig<-function(ci.dat,loci.dat){
	bal<-NULL
	pos<-NULL
	loci.dat<-loci.dat[loci.dat$Het>0,]
	for(ii in 1:(nrow(ci.dat)-1)){
		fst.het<-loci.dat[loci.dat$Het >= ci.dat[ii,"Het"] & 
			loci.dat$Het <= ci.dat[(ii+1),"Het"],]
		max.bound<-min(ci.dat[ii:(ii+1),4])
		min.bound<-max(ci.dat[ii:(ii+1),2])
		bal<-rbind(bal, fst.het[fst.het$Fst <= min.bound,]) 
		pos<-rbind(pos, fst.het[fst.het$Fst >= max.bound,])
	}
	props<-c((nrow(bal)/nrow(loci.dat)), (nrow(pos)/nrow(loci.dat)),
		((nrow(bal)+nrow(pos))/nrow(loci.dat)))	
	return(props)
}

ci.files<-list.files(pattern="genepop.ci")
loci.files<-list.files(pattern="genepop.loci")
plot.ci.files<-ci.files[c(grep(".d2.",ci.files,fixed=T),
	grep(".d5.",ci.files,fixed=T), grep(".d10.",ci.files,fixed=T),
	grep(".d50.",ci.files,fixed=T))]
plot.loci.files<-loci.files[c(grep(".d2.",loci.files,fixed=T),
	grep(".d5.",loci.files,fixed=T), grep(".d10.",loci.files,fixed=T),
	grep(".d50.",loci.files,fixed=T))]


s2.ci<-plot.ci.files[grep(".s2.",plot.ci.files,fixed=T)]
s5.ci<-plot.ci.files[grep(".s5.",plot.ci.files,fixed=T)]
s10.ci<-plot.ci.files[grep(".s10.",plot.ci.files,fixed=T)]
s20.ci<-plot.ci.files[grep(".s20.",plot.ci.files,fixed=T)]
s2.loci<-plot.loci.files[grep(".s2.",plot.loci.files,fixed=T)]
s5.loci<-plot.loci.files[grep(".s5.",plot.loci.files,fixed=T)]
s10.loci<-plot.loci.files[grep(".s10.",plot.loci.files,fixed=T)]
s20.loci<-plot.loci.files[grep(".s20.",plot.loci.files,fixed=T)]

ci.list<-c(s2.ci[grep("Nm0.1",s2.ci,fixed=T)],
	s5.ci[grep("Nm0.1",s5.ci,fixed=T)],
	s10.ci[grep("Nm0.1",s10.ci,fixed=T)],
	s20.ci[grep("Nm0.1",s20.ci,fixed=T)])
loci.list<-c(s2.loci[grep("Nm0.1",s2.loci,fixed=T)],
	s5.loci[grep("Nm0.1",s5.loci,fixed=T)],
	s10.loci[grep("Nm0.1",s10.loci,fixed=T)],
	s20.loci[grep("Nm0.1",s20.loci,fixed=T)])

ss<-c(".s2.",".s5.",".s10.",".s20")
ds<-c(".d2.",".d5.",".d10.",".d50.")


png("../Fig5_jaggedPattern.png",height=225,width=225,units="mm",res=300)
pdf("../Fig5_jaggedPattern.pdf")
par(mfrow=c(4,4), oma=c(2,3,1,1),mar=c(1,1,1,0),mgp=c(3,0.65,0),cex=0.5)
for(i in 1:length(ds)){
	for(j in 1:length(ss)){
		ci.dat<-read.delim(ci.list[grep(ds[i], ci.list,fixed=T)[
			grep(ds[i], ci.list,fixed=T) %in% 
			grep(ss[j],ci.list,fixed=T)]])
		sig.dat<-read.delim(loci.list[grep(ds[i],loci.list,fixed=T)[
			grep(ds[i], loci.list,fixed=T) %in% 
			grep(ss[j],loci.list,fixed=T)]])
		plot(sig.dat$Het, sig.dat$Fst,las=1,ylim=c(0,1),
			xlab="",ylab="",pch=19, xaxt="n",yaxt="n")
		points(ci.dat$Het,ci.dat[,2],col="red",type="l")
		points(ci.dat$Het, ci.dat[,4],col="red",type="l")

		axis(1)
		props<-find.los.sig(ci.dat,sig.dat)
		legend("topleft",bty="n",
			c(paste(round(props[1],3)," balancing"),
			paste(round(props[2],3)," positive"),
			paste(round(props[3],3)," total")))
		if(i == 1 & j == 1){
			mtext("2 Demes",2,cex=0.5,line=2)
			mtext("2 Samples",3,cex=0.5)
		}
		if(i == 1 & j == 2){
			mtext("5 Samples",3,cex=0.5)
		}
		if(i == 1 & j == 3){
			mtext("10 Samples",3,cex=0.5)
		}
		if(i == 1 & j == 4){
			mtext("20 Samples",3,cex=0.5)
		}
		if(i == 2 & j == 1){
			mtext("5 Demes",2,cex=0.5,line=2)
		}
		if(i == 3 & j == 1){
			mtext("10 Demes",2,cex=0.5,line=2)
		}
		if(i == 4 & j == 1){
			mtext ("50 Demes",2,cex=0.5,line=2)
		}
		#if(i == 4){
		#	axis(1)}
		if(j == 1)
			axis(2,las=1)
	}
}
mtext(expression(italic(H)[italic(T)]),1,line=0.5,outer=T,cex=0.5)
mtext(expression(italic(F)[ST]),2,outer=T,line=2,cex=0.5)
dev.off()

############################################################################
#PLOT FIG 6-REVISIONS
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis_selection")
source("../../fhetboot/R/fhetboot.R")
sel.list<-c("Nm1.d2.s20.ds0.output.txt",
	"Nm1.d2.s20.ds0.01.output.txt",
	"Nm1.d2.s20.ds0.1.output.txt",
	"Nm1.d2.s20.ds0.5.output.txt",
	"Nm1.d5.s20.ds0.output.txt",
	"Nm1.d5.s20.ds0.01.output.txt",
	"Nm1.d5.s20.ds0.1.output.txt",
	"Nm1.d5.s20.ds0.5.output.txt")
sel.gpop.cis<-c("Nm1.d2.s20.ds0.genepop.ci",
	"Nm1.d2.s20.ds0.01.genepop.ci",
	"Nm1.d2.s20.ds0.1.genepop.ci",
	"Nm1.d2.s20.ds0.5.genepop.ci",
	"Nm1.d5.s20.ds0.genepop.ci",
	"Nm1.d5.s20.ds0.01.genepop.ci",
	"Nm1.d5.s20.ds0.1.genepop.ci",
	"Nm1.d5.s20.ds0.5.genepop.ci")
sel.gpop<-list(
	my.read.genepop("Nm1.d2.s20.ds0.genepop"),
	my.read.genepop("Nm1.d2.s20.ds0.01.genepop"),
	my.read.genepop("Nm1.d2.s20.ds0.1.genepop"),
	my.read.genepop("Nm1.d2.s20.ds0.5.genepop"),
	my.read.genepop("Nm1.d5.s20.ds0.genepop"),
	my.read.genepop("Nm1.d5.s20.ds0.01.genepop"),
	my.read.genepop("Nm1.d5.s20.ds0.1.genepop"),
	my.read.genepop("Nm1.d5.s20.ds0.5.genepop"))
fsts<-lapply(sel.gpop,calc.actual.fst)
names(fsts)<-c("Nm1.d2.s20.ds0.genepop","Nm1.d2.s20.ds0.01.genepop",
	"Nm1.d2.s20.ds0.1.genepop","Nm1.d2.s20.ds0.5.genepop" ,
	"Nm1.d5.s20.ds0.genepop","Nm1.d5.s20.ds0.01.genepop",
	"Nm1.d5.s20.ds0.1.genepop","Nm1.d5.s20.ds0.5.genepop")
sig.loci<-c("","Nm1.d2.s20.ds0.01.sigloci.txt",
	"Nm1.d2.s20.ds0.1.sigloci.txt",
	"Nm1.d2.s20.ds0.5.sigloci.txt","",
	"Nm1.d5.s20.ds0.01.sigloci.txt",
	"Nm1.d5.s20.ds0.1.sigloci.txt",
	"Nm1.d5.s20.ds0.5.sigloci.txt")



sels<-c(".ds0.o",".ds0.01.",".ds0.1.",".ds0.5")
selsg<-c(".ds0.g",".ds0.01.",".ds0.1.",".ds0.5")
ds<-c(".d2.",".d5.")

#calculate cis from fhetboot if they haven't already
sel.ci<-as.list(rep("",length(sel.gpop)))
for(i in 1:length(sel.gpop)){
	gpop<-my.read.genepop(sel.gpop[i])
	fsts<-calc.actual.fst(gpop)
	boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop))))
	avg.ci95<-ci.means(boot.out[[2]])
	outdat<-data.frame(Het=rownames(sel.ci[[i]][[1]]),
		Low95=sel.ci[[i]][[1]],High95=sel.ci[[i]][[2]])
	write.csv(outdat,paste(sel.gpop[i],"fhetboot.95ci.csv",sep="."))
	sel.ci[[i]]<-avg.ci95
}
names(sel.ci)<-sel.gpop


ci.list<-c("Nm1.d2.s20.ds0.genepop.fhetboot.95ci.csv",
	"Nm1.d2.s20.ds0.01.genepop.fhetboot.95ci.csv",
	"Nm1.d2.s20.ds0.1.genepop.fhetboot.95ci.csv",
	"Nm1.d2.s20.ds0.5.genepop.fhetboot.95ci.csv",
	"Nm1.d5.s20.ds0.genepop.fhetboot.95ci.csv",
	"Nm1.d5.s20.ds0.01.genepop.fhetboot.95ci.csv",
	"Nm1.d5.s20.ds0.1.genepop.fhetboot.95ci.csv",
	"Nm1.d5.s20.ds0.5.genepop.fhetboot.95ci.csv")

png("../Fig6_fhetboot.png",height=169,width=225,units="mm",res=300)
pdf("../Fig6_fhetboot.pdf")
par(mfrow=c(2,4), oma=c(2,3,3,1),mar=c(1,1.5,1,1),mgp=c(3,0.65,0),cex=1)
for(i in 1:length(ds)){
	for(j in 1:length(sels)){
		ci.dat<-read.csv(ci.list[grep(ds[i], ci.list,fixed=T)[
			grep(ds[i], ci.list,fixed=T) %in% 
			grep(selsg[j],ci.list,fixed=T)]])
		sig.dat<-fsts[[grep(ds[i],names(fsts),fixed=T)[
			grep(ds[i], names(fsts),fixed=T) %in% 
			grep(selsg[j],names(fsts),fixed=T)]]]
		gci.dat<-read.delim(sel.gpop.cis[grep(ds[i], sel.gpop.cis,fixed=T)[
			grep(ds[i], sel.gpop.cis,fixed=T) %in% 
			grep(selsg[j],sel.gpop.cis,fixed=T)]])
		plot(sig.dat$Ht,sig.dat$Fst,las=1,xlab="",ylab="",
			pch=19,xaxt='n',yaxt='n')
		points(ci.dat$Het,ci.dat$Low95,col="dodgerblue",type="l",lty=1,lwd=2)
		points(ci.dat$Het,ci.dat$High95,col="dodgerblue",type="l",
			lty=1,lwd=2)
		points(gci.dat$Het,gci.dat[,2],col="red",type="l",
			lty=2,lwd=2)
		points(gci.dat$Het,gci.dat[,4],col="red",type="l",
			lty=2,lwd=2)
		if(j>1){
			sig<-read.delim(sig.loci[grep(ds[i],sig.loci,fixed=T)[
				grep(ds[i], sig.loci,fixed=T) %in% 
				grep(sels[j],sig.loci,fixed=T)]],header=F)
			points(sig.dat[sig.dat$Locus %in% sig$V1,c("Ht","Fst")],
				col="red",pch=8)
		}
		axis(1)
		axis(2,las=1)
		#props<-find.los.sig(ci.dat,sig.dat)
		#legend("topleft",bty="n",
		#	c(paste(round(props[1],3)," balancing"),
		#	paste(round(props[2],3)," positive"),
		#	paste(round(props[3],3)," total")))
		if(i == 1 & j == 1){
			mtext("2 Demes",2,cex=1,line=2)
			mtext(expression(italic(s)~"="~0),3,cex=1)
		}
		if(i == 1 & j == 2){
			mtext(expression(italic(s)~"="~0.01),3,cex=1)
		}
		if(i == 1 & j == 3){
			mtext(expression(italic(s)~"="~0.1),3,cex=1)
		}
		if(i == 1 & j == 4){
			mtext(expression(italic(s)~"="~0.5),3,cex=1)
		}
		if(i == 2 & j == 1){
			mtext("5 Demes",2,cex=1,line=2)
		}
		#if(j == 1)
		#	axis(2,las=1)
	}
}
mtext(expression(italic(H)[italic(T)]),1,line=0.5,outer=T,cex=1)
mtext(expression(italic(F)[ST]),2,outer=T,line=1.5,cex=1)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("dodgerblue","red","red"),pch=c("","","*"),lty=c(1,2,0),
	c("fhetboot","LOSITAN","Selected Loci"),cex=1,
	bg="white",ncol=3,box.lty=0)

dev.off()

############################################################################
#PLOT FIG 7-REVISIONS (HIGHER NM)
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis_selection")
source("../../fhetboot/R/fhetboot.R")
sel.list<-c("Nm10.d2.s20.ds0.output.txt",
	"Nm10.d2.s20.ds0.01.output.txt",
	"Nm10.d2.s20.ds0.1.output.txt",
	"Nm10.d2.s20.ds0.5.output.txt",
	"Nm10.d5.s20.ds0.output.txt",
	"Nm10.d5.s20.ds0.01.output.txt",
	"Nm10.d5.s20.ds0.1.output.txt",
	"Nm10.d5.s20.ds0.5.output.txt")
sel.gpop.cis<-c("Nm10.d2.s20.ds0.genepop.ci",
	"Nm10.d2.s20.ds0.01.genepop.ci",
	"Nm10.d2.s20.ds0.1.genepop.ci",
	"Nm10.d2.s20.ds0.5.genepop.ci",
	"Nm10.d5.s20.ds0.genepop.ci",
	"Nm10.d5.s20.ds0.01.genepop.ci",
	"Nm10.d5.s20.ds0.1.genepop.ci",
	"Nm10.d5.s20.ds0.5.genepop.ci")
sel.gpop<-list(
	my.read.genepop("Nm10.d2.s20.ds0.genepop"),
	my.read.genepop("Nm10.d2.s20.ds0.01.genepop"),
	my.read.genepop("Nm10.d2.s20.ds0.1.genepop"),
	my.read.genepop("Nm10.d2.s20.ds0.5.genepop"),
	my.read.genepop("Nm10.d5.s20.ds0.genepop"),
	my.read.genepop("Nm10.d5.s20.ds0.01.genepop"),
	my.read.genepop("Nm10.d5.s20.ds0.1.genepop"),
	my.read.genepop("Nm10.d5.s20.ds0.5.genepop"))
fsts<-lapply(sel.gpop,calc.actual.fst)
names(fsts)<-c("Nm10.d2.s20.ds0.genepop","Nm10.d2.s20.ds0.01.genepop",
	"Nm10.d2.s20.ds0.1.genepop","Nm10.d2.s20.ds0.5.genepop" ,
	"Nm10.d5.s20.ds0.genepop","Nm10.d5.s20.ds0.01.genepop",
	"Nm10.d5.s20.ds0.1.genepop","Nm10.d5.s20.ds0.5.genepop")
sig.loci<-c("","Nm10.d2.s20.ds0.01.sigloci.txt",
	"Nm10.d2.s20.ds0.1.sigloci.txt",
	"Nm10.d2.s20.ds0.5.sigloci.txt","",
	"Nm10.d5.s20.ds0.01.sigloci.txt",
	"Nm10.d5.s20.ds0.1.sigloci.txt",
	"Nm10.d5.s20.ds0.5.sigloci.txt")



sels<-c(".ds0.o",".ds0.01.",".ds0.1.",".ds0.5")
selsg<-c(".ds0.g",".ds0.01.",".ds0.1.",".ds0.5")
ds<-c(".d2.",".d5.")

names<-c("Nm10.d2.s20.ds0.genepop","Nm10.d2.s20.ds0.01.genepop",
	"Nm10.d2.s20.ds0.1.genepop","Nm10.d2.s20.ds0.5.genepop" ,
	"Nm10.d5.s20.ds0.genepop","Nm10.d5.s20.ds0.01.genepop",
	"Nm10.d5.s20.ds0.1.genepop","Nm10.d5.s20.ds0.5.genepop")

#calculate cis from fhetboot if they haven't already
sel.ci<-as.list(rep("",length(sel.gpop)))
for(i in 1:length(sel.gpop)){
	gpop<-sel.gpop[[i]]
	#fsts<-fsts[[i]]
	boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop))))
	avg.ci95<-ci.means(boot.out[[2]])
	outdat<-data.frame(Het=rownames(avg.ci95[[1]]),
		Low95=avg.ci95[[1]],High95=avg.ci95[[2]])
	write.csv(outdat,paste(names[i],"fhetboot.95ci.csv",sep="."))
	sel.ci[[i]]<-avg.ci95
}
names(sel.ci)<-names


ci.list<-c("Nm10.d2.s20.ds0.genepop.fhetboot.95ci.csv",
	"Nm10.d2.s20.ds0.01.genepop.fhetboot.95ci.csv",
	"Nm10.d2.s20.ds0.1.genepop.fhetboot.95ci.csv",
	"Nm10.d2.s20.ds0.5.genepop.fhetboot.95ci.csv",
	"Nm10.d5.s20.ds0.genepop.fhetboot.95ci.csv",
	"Nm10.d5.s20.ds0.01.genepop.fhetboot.95ci.csv",
	"Nm10.d5.s20.ds0.1.genepop.fhetboot.95ci.csv",
	"Nm10.d5.s20.ds0.5.genepop.fhetboot.95ci.csv")

png("../Fig7_fhetboot_Nm10.png",height=169,width=225,units="mm",res=300)
pdf("../Fig7_fhetboot_Nm10.pdf")
par(mfrow=c(2,4), oma=c(2,3,3,1),mar=c(1,1.5,1,1),mgp=c(3,0.65,0),cex=1)
for(i in 1:length(ds)){
	for(j in 1:length(sels)){
		ci.dat<-read.csv(ci.list[grep(ds[i], ci.list,fixed=T)[
			grep(ds[i], ci.list,fixed=T) %in% 
			grep(selsg[j],ci.list,fixed=T)]])
		sig.dat<-fsts[[grep(ds[i],names(fsts),fixed=T)[
			grep(ds[i], names(fsts),fixed=T) %in% 
			grep(selsg[j],names(fsts),fixed=T)]]]
		gci.dat<-read.delim(sel.gpop.cis[grep(ds[i], sel.gpop.cis,fixed=T)[
			grep(ds[i], sel.gpop.cis,fixed=T) %in% 
			grep(selsg[j],sel.gpop.cis,fixed=T)]])
		plot(sig.dat$Ht,sig.dat$Fst,las=1,xlab="",ylab="",
			pch=19,xaxt='n',yaxt='n')
		axis(2,las=1)
		points(ci.dat$Het,ci.dat$Low95,col="dodgerblue",type="l",lty=1,lwd=2)
		points(ci.dat$Het,ci.dat$High95,col="dodgerblue",type="l",
			lty=1,lwd=2)
		points(gci.dat$Het,gci.dat[,2],col="red",type="l",
			lty=2,lwd=2)
		points(gci.dat$Het,gci.dat[,4],col="red",type="l",
			lty=2,lwd=2)
		if(j>1){
			sig<-read.delim(sig.loci[grep(ds[i],sig.loci,fixed=T)[
				grep(ds[i], sig.loci,fixed=T) %in% 
				grep(sels[j],sig.loci,fixed=T)]],header=F)
			points(sig.dat[sig.dat$Locus %in% sig$V1,c("Ht","Fst")],
				col="red",pch=8)
		}
		axis(1)
		#props<-find.los.sig(ci.dat,sig.dat)
		#legend("topleft",bty="n",
		#	c(paste(round(props[1],3)," balancing"),
		#	paste(round(props[2],3)," positive"),
		#	paste(round(props[3],3)," total")))
		if(i == 1 & j == 1){
			mtext("2 Demes",2,cex=1,line=2)
			mtext(expression(italic(s)~"="~0),3,cex=1)
		}
		if(i == 1 & j == 2){
			mtext(expression(italic(s)~"="~0.01),3,cex=1)
		}
		if(i == 1 & j == 3){
			mtext(expression(italic(s)~"="~0.1),3,cex=1)
		}
		if(i == 1 & j == 4){
			mtext(expression(italic(s)~"="~0.5),3,cex=1)
		}
		if(i == 2 & j == 1){
			mtext("5 Demes",2,cex=1,line=2)
		}
		#if(j == 1)
		#	axis(2,las=1)
	}
}
mtext(expression(italic(H)[italic(T)]),1,line=0.5,outer=T,cex=1)
mtext(expression(italic(F)[ST]),2,outer=T,line=1.5,cex=1)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("dodgerblue","red","red"),pch=c("","","*"),lty=c(1,2,0),
	c("fhetboot","LOSITAN","Selected Loci"),cex=1,
	bg="white",ncol=3,box.lty=0)

dev.off()


############################################################################
#PLOT FIG S2-REVISIONS (FORMER FIG 4)
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis")

pop.list<-list(
	read.delim("Nm1.d2.n100.output.txt"),
	read.delim("Nm1.d2.n500.fig.output.txt"),
	read.delim("Nm1.d2.n1000.fig.output.txt"),
	read.delim("Nm1.d5.n100.fig.output.txt"),
	read.delim("Nm1.d5.n500.fig.output.txt"),
	read.delim("Nm1.d5.n1000.fig.output.txt"))
names(pop.list)<-c("Nm1.d2.n100","Nm1.d2.n500","Nm1.d2.n1000",
	"Nm1.d5.n100","Nm1.d5.n500","Nm1.d5.n1000")
pop.cis<-list(
	read.delim("Nm1.d2.n100.genepop.ci"),
	read.delim("Nm1.d2.n500.fig.genepop.ci"),
	read.delim("Nm1.d2.n1000.fig.genepop.ci"),
	read.delim("Nm1.d5.n100.fig.genepop.ci"),
	read.delim("Nm1.d5.n500.fig.genepop.ci"),
	read.delim("Nm1.d5.n1000.fig.genepop.ci"))

png("FigS2_popsize.png",height=225,width=169,units="mm",res=300)
pdf("FigS2_popsize.pdf")
par(mfrow = c(2, 3),cex = 0.6,mar = c(0, 0, 0, 0), 
	oma = c(4, 4.5, 1.5, 0.5), tcl = -0.25,mgp = c(2, 0.6, 0))
for (i in 1:length(pop.list)) {
	#plot(seq(0,0.6,0.06),seq(0,1,0.1), axes = FALSE, type = "n")
	y.min<-min(min(pop.cis[[i]][,2]),min(pop.list[[i]]$WrightsFst))
	y.max<-(max(pop.list[[i]]$WrightsFst)+2*sd(pop.list[[i]]$WrightsFst))
	if(y.max>1)
		y.max<-1
	plot(pop.list[[i]]$Ht,pop.list[[i]]$WrightsFst,las=1,,xlab="",ylab="",
		col="black",pch=19,xaxt="n", xlim=c(0,0.6),ylim=c(0,y.max))
	points(pop.cis[[i]]$Het,pop.cis[[i]][,2],col="red",type="l")
	points(pop.cis[[i]]$Het, pop.cis[[i]][,4],col="red",type="l")

	Nm<-strsplit(strsplit(names(pop.list)[i], "Nm")[[1]][2],"[._]")[[1]][1]
	exp.fst<-round(1/((4*as.numeric(Nm))+1),digits=4)
	avg.fst<-round(mean(pop.list[[i]]$WrightsFst),digits=4)
	leg.nms<-c(bquote(Exp.~italic(F)[ST]~"="~.(exp.fst)),
		bquote(Mean~italic(F)[ST]~"="~.(avg.fst)))
	legend("topleft",legend=as.expression(leg.nms),bty="n")

	if (i %in% c(4,5,6))
		axis(1, at = seq(0,0.6,0.1))
	if(i==1){
		mtext("N = 100", 3,cex=.75)
		mtext("2 Demes",2,line=1.75,cex=.75)}
	if(i==2)
		mtext("N = 500",3,cex=.75)
	if(i==3){
		mtext("N=1000",3,cex=.75)
	}
	if(i==4)
		mtext("5 Demes",2,line=1.75,cex=.75)
}
mtext(expression(italic(H)[italic(T)]),1,outer=T,line=2,cex=.75)
mtext(expression(Wright~"'"~s~italic(F)[ST]),2,outer=T,line=3,cex=.75)
dev.off()

############################################################################
#PLOT FIGS S3-S6-REVISIONS (FORMER S1-S4)
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis_genepop")
ci.files<-list.files(pattern="genepop.ci")
loci.files<-list.files(pattern="genepop.loci")
plot.ci.files<-ci.files[c(grep(".d2.",ci.files,fixed=T),
	grep(".d5.",ci.files,fixed=T), grep(".d10.",ci.files,fixed=T),
	grep(".d50.",ci.files,fixed=T))]
plot.loci.files<-loci.files[c(grep(".d2.",loci.files,fixed=T),
	grep(".d5.",loci.files,fixed=T), grep(".d10.",loci.files,fixed=T),
	grep(".d50.",loci.files,fixed=T))]


s2.ci<-plot.ci.files[grep(".s2.",plot.ci.files,fixed=T)]
s5.ci<-plot.ci.files[grep(".s5.",plot.ci.files,fixed=T)]
s10.ci<-plot.ci.files[grep(".s10.",plot.ci.files,fixed=T)]
s20.ci<-plot.ci.files[grep(".s20.",plot.ci.files,fixed=T)]
s2.loci<-plot.loci.files[grep(".s2.",plot.loci.files,fixed=T)]
s5.loci<-plot.loci.files[grep(".s5.",plot.loci.files,fixed=T)]
s10.loci<-plot.loci.files[grep(".s10.",plot.loci.files,fixed=T)]
s20.loci<-plot.loci.files[grep(".s20.",plot.loci.files,fixed=T)]

nms<-c("Nm1.","Nm10.")
ds<-c(".d2.",".d5.",".d10.",".d50.")

find.los.sig<-function(ci.dat,loci.dat){
	bal<-NULL
	pos<-NULL
	loci.dat<-loci.dat[loci.dat$Het>0,]
	for(ii in 1:(nrow(ci.dat)-1)){
		fst.het<-loci.dat[loci.dat$Het >= ci.dat[ii,"Het"] & 
			loci.dat$Het <= ci.dat[(ii+1),"Het"],]
		max.bound<-min(ci.dat[ii:(ii+1),4])
		min.bound<-max(ci.dat[ii:(ii+1),2])
		bal<-rbind(bal, fst.het[fst.het$Fst <= min.bound,]) 
		pos<-rbind(pos, fst.het[fst.het$Fst >= max.bound,])
	}
	props<-c((nrow(bal)/nrow(loci.dat)), (nrow(pos)/nrow(loci.dat)),
		((nrow(bal)+nrow(pos))/nrow(loci.dat)))	
	return(props)
}

plot.lositan<-function(ci.list,loci.list,out.name,colgroups,rowgroups,pdf=F){
	if(pdf==FALSE){
		png(paste(out.name,"png",sep="."),height=225,width=169,
			units="mm",res=300)
	} else{
		pdf(paste(out.name,"pdf",sep="."))
	}
	numcol<-length(colgroups)
	numrow<-length(rowgroups)
	par(mfrow=c(numrow,numcol), oma=c(2,3,1,1),mar=c(1,1,1,0),
		mgp=c(3,0.65,0),cex=0.5)
	for(i in 1:length(rowgroups)){
		for(j in 1:length(colgroups)){
			ci.dat<-read.delim(ci.list[grep(rowgroups[i], ci.list,fixed=T)[
				grep(rowgroups[i], ci.list,fixed=T) %in% 
				grep(colgroups[j],ci.list,fixed=T)]])
			sig.dat<-read.delim(loci.list[grep(rowgroups[i],loci.list,fixed=T)[
				grep(rowgroups[i], loci.list,fixed=T) %in% 
				grep(colgroups[j],loci.list,fixed=T)]])
			plot(sig.dat$Het, sig.dat$Fst,las=1,ylim=c(0,1),
				xlab="",ylab="",pch=19, xaxt="n",yaxt="n")
			points(ci.dat$Het,ci.dat[,2],col="red",type="l")
			points(ci.dat$Het, ci.dat[,4],col="red",type="l")
			axis(1)
			props<-find.los.sig(ci.dat,sig.dat)
			legend("topleft",bty="n",
				c(paste(round(props[1],3)," balancing"),
				paste(round(props[2],3)," positive"),
				paste(round(props[3],3)," total")))
			cnames<-gsub("Nm(\\d+).","\\1",nms)
			rnames<-gsub(".d(\\d+).","\\1",ds)
			if(i == 1 & j == 1){
				mtext(bquote(.(rnames[i])~Demes),2,cex=0.5,line=2)
				mtext(bquote(Nm~"="~.(cnames[j])),3,cex=0.5)
			}
			if(i == 1 & j == 2){
				mtext(bquote(Nm~"="~.(cnames[j])),3,cex=0.5)
			}
			if(i == 1 & j == 3){
				mtext(bquote(Nm~"="~.(cnames[j])),3,cex=0.5)
			}
			if(i == 2 & j == 1){
				mtext(bquote(.(rnames[i])~Demes),2,cex=0.5,line=2)
			}
			if(i == 3 & j == 1){
				mtext(bquote(.(rnames[i])~Demes),2,cex=0.5,line=2)
			}
			if(i == 4 & j == 1){
				mtext(bquote(.(rnames[i])~Demes),2,cex=0.5,line=2)
			}
			#if(i == 4){
			#	axis(1)}
			if(j == 1)
				axis(2,las=1)
		}
	}
	mtext(expression(italic(H)[italic(T)]),1,line=0.5,outer=T,cex=0.5)
	mtext(expression(italic(F)[ST]),2,outer=T,line=2,cex=0.5)
	dev.off()
}

plot.lositan(s2.ci,s2.loci,"../S3.s2",nms,ds,pdf=TRUE)
plot.lositan(s5.ci,s5.loci,"../S5.s5",nms,ds,pdf=TRUE)
plot.lositan(s10.ci,s10.loci,"../S10.s10",nms,ds,pdf=TRUE)
plot.lositan(s20.ci,s20.loci,"../S20.s20",nms,ds,pdf=TRUE)

plot.lositan(s2.ci,s2.loci,"../S3.s2",nms,ds,pdf=FALSE)
plot.lositan(s5.ci,s5.loci,"../S5.s5",nms,ds,pdf=FALSE)
plot.lositan(s10.ci,s10.loci,"../S10.s10",nms,ds,pdf=FALSE)
plot.lositan(s20.ci,s20.loci,"../S20.s20",nms,ds,pdf=FALSE)


############################################################################
#EXPLORING LOSITAN PARAMETER SPACE
############################################################################
setwd("E://ubuntushare//fst_outliers//results//numerical_analysis_genepop")

ci.files<-list.files(pattern="genepop.ci")
loci.files<-list.files(pattern="genepop.loci")

ci.files<-sort(ci.files)
loci.files<-sort(loci.files)
sig.loci<-list()
num.sig<-NULL

for(i in 1:length(ci.files)){
	ci.dat<-read.delim(ci.files[i])
	ci.dat<-ci.dat[order(ci.dat[,1]),]
	loci.dat<-read.delim(loci.files[i])
	bal<-NULL
	pos<-NULL
	loci.dat<-loci.dat[loci.dat$Fst >= 0,]
	for(ii in 1:(nrow(ci.dat)-1)){
		fst.het<-loci.dat[loci.dat$Het >= ci.dat[ii,"Het"] & 
			loci.dat$Het <= ci.dat[(ii+1),"Het"],]
		max.bound<-min(ci.dat[ii:(ii+1),4])
		min.bound<-max(ci.dat[ii:(ii+1),2])
		bal<-rbind(bal, fst.het[fst.het$Fst <= min.bound,]) 
		pos<-rbind(pos, fst.het[fst.het$Fst >= max.bound,])
	}
	sig.loci[[i]]<-list(as.data.frame(bal),as.data.frame(pos))
	names(sig.loci[[i]])<-c("bal","pos")
	num.sig<-as.data.frame(rbind(num.sig,rbind(lapply(sig.loci[[i]],nrow))))
}
names(sig.loci)<-ci.files

num.sig<-as.data.frame(cbind(as.numeric(num.sig$bal),as.numeric(num.sig$pos)))
rownames(num.sig)<-ci.files
colnames(num.sig)<-c("bal","pos")
num.sig$total<-rowSums(num.sig)
write.table(num.sig,"sig.loci.sim.LOSITAN.txt",quote=F,col.names=T,
	row.names=T)

#total number of loci = 2000
#If I expect 0.05 to be outliers due to threshold, should have ~100 loci.
#Ideally they'd be evenly distributed but I'm not sure that's really true.

ok.sig<-num.sig[num.sig$total >50 & num.sig$total < 150,]
los.sig<-sig.loci[names(sig.loci) %in% rownames(ok.sig)]
plot.titles<-c("Nm=10,d=10,s=10","Nm=10,d=10,s=2","Nm=10,d=10,s=5","Nm=10,d=100,s=20",
	"Nm=10,d=100,s=10","Nm=10,d=100,s=5","Nm=10,d=2,s=2","Nm=10,d=2,s=5",
	"Nm=10,d=25,s=10","Nm=10,d=25,s=20","Nm=10,d=25,s=5","Nm=10,d=5,s=2",
	"Nm=10,d=5,s=5","Nm=10,d=50,s=10","Nm=10,d=50,s=20","Nm=10,d=50,s=5",
	"Nm=10,d=75,s=10","Nm=10,d=75,s=20","Nm=10,d=75,s=5")
#plot the significant ones
png("fig6.png",height=282,width=225,units="mm",res=300)
par(mfrow=c(5,4), oma=c(2,3,1,1),mar=c(1,2,1,0),mgp=c(3,0.65,0),cex=0.5)
for(i in 1:nrow(ok.sig)){
	sig.name<-paste(strsplit(rownames(ok.sig)[i],"genepop.ci"),
		"genepop.loci",sep="")
	ci.dat<-read.delim(rownames(ok.sig)[i])
	sig.dat<-read.delim(sig.name)
	#sig.file<-read.delim(sig.name)
	#plot(sig.file$Ht,sig.file$WrightsFst,ylim=c(0,1),xlim=c(0,0.5))
	plot(sig.dat$Het, sig.dat$Fst,ylim=c(0,1),las=1,xlab="",ylab="",
		pch=19,yaxt="n")
	points(ci.dat$Het,ci.dat[,2],col="yellow",type="l")
	points(ci.dat$Het, ci.dat[,4],col="red",type="l")
	props<-find.los.sig(ci.dat,sig.dat)
	legend("topleft",bty="n",c(paste(round(props[1],3)," balancing"),
		paste(round(props[2],3)," positive"),
		paste(round(props[3],3)," total")))
	legend("topright",plot.titles[i],bty="n")
	if(i %in% c(1,5,9,13,17))
		axis(2,las=1)
	axis(1)
}
mtext("Heterozygosity",1,line=0.5,outer=T,cex=0.5)
mtext("Fst",2,outer=T,line=1,cex=0.5)
dev.off()

#testing genepop
s2d2.gene<-read.delim("E://ubuntushare//numerical_analysis_genepop//samp2deme2.txt")
s2d2.loci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp2deme2genepop.loci")
s2d2.ci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp2deme2genepop.ci")
s2d2.mydat<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10.samp2.deme2.sampledpops.txt")
s2d100.gene<-read.delim("E://ubuntushare//numerical_analysis_genepop//samp2deme100.txt")
s2d100.loci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp2deme100genepop.loci")
s2d100.ci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp2deme100genepop.ci")
s2d100.mydat<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10.samp2.deme100.sampledpops.txt")
s5d2.gene<-read.delim("E://ubuntushare//numerical_analysis_genepop//samp5deme2.txt")
s5d2.loci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp5deme2genepop.loci")
s5d2.ci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp5deme2genepop.ci")
s5d2.mydat<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10.samp5.deme2.sampledpops.txt")
s5d100.gene<-read.delim("E://ubuntushare//numerical_analysis_genepop//samp5deme100.txt")
s5d100.loci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp5deme100genepop.loci")
s5d100.ci<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10samp5deme100genepop.ci")
s5d100.mydat<-read.delim("E://ubuntushare//numerical_analysis_genepop//loc10.samp5.deme100.sampledpops.txt")

png("genepop.comparison.png",height=225,width=169,units="mm",res=300)
par(mfrow=c(4,3), oma=c(2,3,1,1),mar=c(1,1,1,0),mgp=c(3,0.65,0),cex=0.5)
plot(s2d2.gene$Het, s2d2.gene$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
mtext("Genepop",3)
mtext("2 Demes, 2 Samples",2)
plot(s2d2.loci$Het, s2d2.loci$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
points(s2d2.ci$Het,s2d2.ci[,2],col="yellow",type="l")
points(s2d2.ci$Het, s2d2.ci[,4],col="red",type="l")
mtext("Lositan", 3)
plot(s2d2.mydat$Ht, s2d2.mydat$WrightsFst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
mtext("My Calcs",3)
			
plot(s2d100.gene$Het, s2d100.gene$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
mtext("100 Demes, 2 Samples",2)
plot(s2d100.loci$Het, s2d100.loci$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
points(s2d100.ci$Het,s2d100.ci[,2],col="yellow",type="l")
points(s2d100.ci$Het, s2d100.ci[,4],col="red",type="l")
plot(s2d100.mydat$Ht, s2d100.mydat$WrightsFst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)

plot(s5d2.gene$Het, s5d2.gene$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
mtext("2 Demes, 5 Samples",2)
plot(s5d2.loci$Het, s5d2.loci$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
points(s5d2.ci$Het,s5d2.ci[,2],col="yellow",type="l")
points(s5d2.ci$Het,s5d2.ci[,4],col="red",type="l")
plot(s5d2.mydat$Ht, s5d2.mydat$WrightsFst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)

plot(s5d100.gene$Het,s5d100.gene$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
mtext("100 Demes, 5 Samples",2)
plot(s5d100.loci$Het, s5d100.loci$Fst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
points(s5d100.ci$Het,s5d100.ci[,2],col="yellow",type="l")
points(s5d100.ci$Het,s5d100.ci[,4],col="red",type="l")
plot(s5d100.mydat$Ht, s5d100.mydat$WrightsFst,las=1,ylim=c(0,1),
	xlab="",ylab="",pch=19)
dev.off()

#my dat

plot.mydat.lositan<-function(ci.list,loci.list,out.name){
	png(out.name,height=225,width=169,units="mm",res=300)
	par(mfrow=c(4,3), oma=c(2,3,1,1),mar=c(1,1,1,0),mgp=c(3,0.65,0),cex=0.5)
	for(i in 1:length(ds)){
		for(j in 1:length(nms)){
			ci.dat<-read.delim(ci.list[grep(ds[i], ci.list,fixed=T)[
				grep(ds[i], ci.list,fixed=T) %in% 
				grep(nms[j],ci.list,fixed=T)]])
			sig.dat<-read.delim(loci.list[grep(ds[i],loci.list,fixed=T)[
				grep(ds[i], loci.list,fixed=T) %in% 
				grep(nms[j],loci.list,fixed=T)]])
			colnames(sig.dat)[2]<-"Het"
			colnames(sig.dat)[4]<-"Fst"
			plot(sig.dat$Het, sig.dat$Fst,las=1,ylim=c(0,1),
				xlab="",ylab="",pch=19, xaxt="n",yaxt="n")
			points(ci.dat$Het,ci.dat[,2],col="yellow",type="l")
			points(ci.dat$Het, ci.dat[,4],col="red",type="l")
			props<-find.los.sig(ci.dat,sig.dat)
			legend("topleft",bty="n",c(paste(props[1]," balancing"),
				paste(props[2]," positive"),
				paste(props[3]," total")))
			if(i == 1 & j == 1){
				mtext("2 Demes",2,cex=0.5,line=2)
				mtext("Nm = 0.1",3,cex=0.5)
			}
			if(i == 1 & j == 2){
				mtext("Nm = 1",3,cex=0.5)
			}
			if(i == 1 & j == 3){
				mtext("Nm = 10",3,cex=0.5)
			}
			if(i == 2 & j == 1){
				mtext("5 Demes",2,cex=0.5,line=2)
			}
			if(i == 3 & j == 1){
				mtext("10 Demes",2,cex=0.5,line=2)
			}
			if(i == 4 & j == 1){
				mtext ("50 Demes",2,cex=0.5,line=2)
			}
			axis(1)
			if(j == 1)
				axis(2,las=1)
		}
	}
	mtext("Heterozygosity",1,line=0.5,outer=T,cex=0.5)
	mtext("Fst",2,outer=T,line=2,cex=0.5)
	dev.off()
}
ci.files<-list.files(pattern="genepop.ci")
loci.files<-list.files(pattern="sampledpops.txt")
plot.ci.files<-ci.files[c(grep(".d2.",ci.files,fixed=T),
	grep(".d5.",ci.files,fixed=T), grep(".d10.",ci.files,fixed=T),
	grep(".d50.",ci.files,fixed=T))]
plot.loci.files<-loci.files[c(grep(".d2.",loci.files,fixed=T),
	grep(".d5.",loci.files,fixed=T), grep(".d10.",loci.files,fixed=T),
	grep(".d50.",loci.files,fixed=T))]


s2.ci<-plot.ci.files[grep(".s2.",plot.ci.files,fixed=T)]
s5.ci<-plot.ci.files[grep(".s5.",plot.ci.files,fixed=T)]
s10.ci<-plot.ci.files[grep(".s10.",plot.ci.files,fixed=T)]
s20.ci<-plot.ci.files[grep(".s20.",plot.ci.files,fixed=T)]
s2.loci<-plot.loci.files[grep(".s2.",plot.loci.files,fixed=T)]
s5.loci<-plot.loci.files[grep(".s5.",plot.loci.files,fixed=T)]
s10.loci<-plot.loci.files[grep(".s10.",plot.loci.files,fixed=T)]
s20.loci<-plot.loci.files[grep(".s20.",plot.loci.files,fixed=T)]

nms<-c("Nm0.1.","Nm1.","Nm10.")
ds<-c(".d2.",".d5.",".d10.",".d50.")

plot.mydat.lositan(s2.ci,s2.loci,"fig6.s2.mydat.png")
plot.mydat.lositan(s5.ci,s5.loci,"fig6.s5.mydat.png")
plot.mydat.lositan(s10.ci,s10.loci,"fig6.s10.mydat.png")
plot.mydat.lositan(s20.ci,s20.loci,"fig6.s20.mydat.png")

############################################################################
#ADAM'S CODE FIGS
############################################################################
setwd("E://ubuntushare//fdist2//AdamsFdist")
plot.adam.lositan<-function(ci.list,loci.list,out.name){
	png(out.name,height=225,width=169,units="mm",res=300)
	par(mfrow=c(4,3), oma=c(2,3,1,1),mar=c(1,1,1,0),mgp=c(3,0.65,0),cex=0.5)
	for(i in 1:length(ds)){
		for(j in 1:length(nms)){
			ci.dat<-read.delim(ci.list[grep(ds[i], ci.list,fixed=T)[
				grep(ds[i], ci.list,fixed=T) %in% 
				grep(nms[j],ci.list,fixed=T)]])
			sig.dat<-read.delim(loci.list[grep(ds[i],loci.list,fixed=T)[
				grep(ds[i], loci.list,fixed=T) %in% 
				grep(nms[j],loci.list,fixed=T)]])
			colnames(sig.dat)[colnames(sig.dat)=="Ht"]<-"Het"
			plot(sig.dat$Het, sig.dat$Fst,las=1,ylim=c(-0.01,1),
				xlab="",ylab="",pch=19, xaxt="n",yaxt="n")
			points(ci.dat$Het,ci.dat[,2],col="yellow",type="l")
			points(ci.dat$Het, ci.dat[,4],col="red",type="l")
			props<-find.los.sig(ci.dat,sig.dat)
			legend("topleft",bty="n",c(paste(props[1]," balancing"),
				paste(props[2]," positive"),
				paste(props[3]," total")))
			if(i == 1 & j == 1){
				mtext("2 Demes",2,cex=0.5,line=2)
				mtext("Nm = 0.1",3,cex=0.5)
			}
			if(i == 1 & j == 2){
				mtext("Nm = 1",3,cex=0.5)
			}
			if(i == 1 & j == 3){
				mtext("Nm = 10",3,cex=0.5)
			}
			if(i == 2 & j == 1){
				mtext("5 Demes",2,cex=0.5,line=2)
			}
			if(i == 3 & j == 1){
				mtext("10 Demes",2,cex=0.5,line=2)
			}
			if(i == 4 & j == 1){
				mtext ("50 Demes",2,cex=0.5,line=2)
			}
			if(i == 4){
				axis(1)}
			if(j == 1)
				axis(2,las=1)
		}
	}
	mtext("Heterozygosity",1,line=0.5,outer=T,cex=0.5)
	mtext("Fst",2,outer=T,line=2,cex=0.5)
	dev.off()
}

ci.files<-list.files(pattern=".ci")
loci.files<-list.files(pattern=".loci")
ci.files<-ci.files[!(ci.files %in% loci.files)]
adam.files<-list.files(pattern="Fdistcheckerresults")
plot.ci.files<-ci.files[c(grep(".d2.",ci.files,fixed=T),
	grep(".d5.",ci.files,fixed=T), grep(".d10.",ci.files,fixed=T),
	grep(".d50.",ci.files,fixed=T))]
plot.loci.files<-loci.files[c(grep(".d2.",loci.files,fixed=T),
	grep(".d5.",loci.files,fixed=T), grep(".d10.",loci.files,fixed=T),
	grep(".d50.",loci.files,fixed=T))]
plot.adam.files<-adam.files[c(grep(".d2.",adam.files,fixed=T),
	grep(".d5.",adam.files,fixed=T), grep(".d10.",adam.files,fixed=T),
	grep(".d50.",adam.files,fixed=T))]

s2.ci<-plot.ci.files[grep(".s2.",plot.ci.files,fixed=T)]
s5.ci<-plot.ci.files[grep(".s5.",plot.ci.files,fixed=T)]
s10.ci<-plot.ci.files[grep(".s10.",plot.ci.files,fixed=T)]
s20.ci<-plot.ci.files[grep(".s20.",plot.ci.files,fixed=T)]
s2.loci<-plot.loci.files[grep(".s2.",plot.loci.files,fixed=T)]
s5.loci<-plot.loci.files[grep(".s5.",plot.loci.files,fixed=T)]
s10.loci<-plot.loci.files[grep(".s10.",plot.loci.files,fixed=T)]
s20.loci<-plot.loci.files[grep(".s20.",plot.loci.files,fixed=T)]
s2.adam<-plot.adam.files[grep(".s2.",plot.adam.files,fixed=T)]
s5.adam<-plot.adam.files[grep(".s5.",plot.adam.files,fixed=T)]
s10.adam<-plot.adam.files[grep(".s10.",plot.adam.files,fixed=T)]
s20.adam<-plot.adam.files[grep(".s20.",plot.adam.files,fixed=T)]

nms<-c("Nm0.1.","Nm1.","Nm10.")
ds<-c(".d2.",".d5.",".d10.",".d50.")

plot.adam.lositan(s2.ci,s2.adam,"fig6.s2.adam.png")
plot.adam.lositan(s5.ci,s5.adam,"fig6.s5.adam.png")
plot.adam.lositan(s10.ci,s10.adam,"fig6.s10.adam.png")
plot.adam.lositan(s20.ci,s20.adam,"fig6.s20.adam.png")


plot.lositan(s2.ci,s2.loci,"fig6.s2.adam.los.png")
plot.lositan(s5.ci,s5.loci,"fig6.s5.adam.los.png")
plot.lositan(s10.ci,s10.loci,"fig6.s10.adam.los.png")
plot.lositan(s20.ci,s20.loci,"fig6.s20.adam.los.png")
