#Author: Sarah P. Flanagan
#Date: 18 September 2015
#Purpose: Create a bunch of plots of the fdist2 results
#I ran fdist2 with a bunch of different parameter settings.

rm(list=ls())

###############################################################################
#FDIST2
###############################################################################

setwd("E://ubuntushare//fst-het")

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

###############################################################################
################NUMERICAL ANALYSIS
###############################################################################
#***EXPLORATORY ANALYSIS***#
setwd("E://ubuntushare//fst-het//numerical_analysis//original_runs")
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
#PLOT FIG 3
############################################################################
setwd("E://Docs//fst-het//results")


dat.list<-list(
	read.delim("Nm0.1.2pops.1sample.output.txt"),
	read.delim("Nm1.2pops.1sample.output.txt"),
	read.delim("Nm10.2pops.1sample.output.txt"),
	read.delim("Nm0.1.5pops.1sample.output.txt"),
	read.delim("Nm1.5pops.1sample.output.txt"),
	read.delim("Nm10.5pops.1sample.output.txt"))
names(dat.list)<-c("Nm0.1.2pops.1sample",
	"Nm1.2pops.1sample","Nm10.2pops.1sample",
	"Nm0.1.5pops.1sample","Nm1.5pops.1sample",
	"Nm10.5pops.1sample")
Nms<-c(0.1,1,10,0.1,1,10)


png("fig3.png",width=169,height=169,units="mm",res=300)
par(mfrow = c(2, 3),cex = 0.6,mar = c(0, 0, 0, 0), 
	oma = c(4, 4, 1.5, 0.5), tcl = -0.25,mgp = c(2, 0.6, 0))
for (i in 1:length(dat.list)) {
	#plot(seq(0,0.6,0.06),seq(0,1,0.1), axes = FALSE, type = "n")
	
	y.max<-(max(dat.list[[i]]$WrightsFst)+2*sd(dat.list[[i]]$WrightsFst))
	if(y.max>1)
		y.max<-1
	plot(dat.list[[i]]$Ht,dat.list[[i]]$WrightsFst,las=1,
		col="black",pch=19,xaxt="n", xlim=c(0,0.6),ylim=c(0,y.max))
	Nm<-Nms[i]
	exp.fst<-round(1/((4*as.numeric(Nm))+1),digits=3)
	avg.fst<-round(mean(dat.list[[i]]$WrightsFst),digits=3)
	legend("topleft",c(paste("Exp. Fst = ",exp.fst),
		paste("Mean Fst = ",avg.fst)),bty="n")
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
mtext("Ht",1,outer=T,line=2,cex=.75)
mtext("Wright's Fst",2,outer=T,line=3,cex=.75)
dev.off()

############################################################################
#PLOT FIG 4
############################################################################
setwd("E://Docs//fst-het//results")

dat.list<-list(
	read.delim("Nm1.2pops.1sample.n100.output.txt"),
	read.delim("Nm1.2pops.1sample.n500.output.txt"),
	read.delim("Nm1.2pops.1sample.n1000.output.txt"),
	read.delim("Nm1.5pops.1sample.n100.output.txt"),
	read.delim("Nm1.5pops.1sample.n500.output.txt"),
	read.delim("Nm1.5pops.1sample.n1000.output.txt"))
names(dat.list)<-c("Nm1.2pops.1sample.n100",
	"Nm1.2pops.1sample.n500","Nm1.2pops.1sample.n1000",
	"Nm1.5pops.1sample.n100","Nm1.5pops.1sample.n500",
	"Nm1.5pops.1sample.n1000")
png("fig4.png",width=169,height=169,units="mm",res=300)
par(mfrow = c(2, 3),cex = 0.6,mar = c(0, 0, 0, 0), 
	oma = c(4, 4, 1.5, 0.5), tcl = -0.25,mgp = c(2, 0.6, 0))
for (i in 1:length(dat.list)) {
	#plot(seq(0,0.6,0.06),seq(0,1,0.1), axes = FALSE, type = "n")
	
	y.max<-(max(dat.list[[i]]$WrightsFst)+2*sd(dat.list[[i]]$WrightsFst))
	if(y.max>1)
		y.max<-1
	plot(dat.list[[i]]$Ht,dat.list[[i]]$WrightsFst,las=1,
		col="black",pch=19,xaxt="n", xlim=c(0,0.6),ylim=c(0,y.max))
	Nm<-strsplit(strsplit(names(dat.list)[i], "Nm")[[1]][2],"[._]")[[1]][1]
	exp.fst<-round(1/((4*as.numeric(Nm))+1),digits=4)
	avg.fst<-round(mean(dat.list[[i]]$WrightsFst),digits=4)
	legend("topleft",c(paste("Exp. Fst = ",exp.fst),
		paste("Mean Fst = ",avg.fst)),bty="n")
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
mtext("Ht",1,outer=T,line=2,cex=.75)
mtext("Wright's Fst",2,outer=T,line=3,cex=.75)
dev.off()

############################################################################
#PLOT FIG 5
############################################################################
setwd("E://Docs//fst-het//results")

dat.list<-list(
	read.delim("Nm1.2pops.1sample.n1000.sampledpops.txt"),
	read.delim("Nm1.2pops.2sample.n1000.sampledpops.txt"),
	read.delim("Nm1.2pops.4sample.n1000.sampledpops.txt"),
	read.delim("Nm1.5pops.1sample.n1000.sampledpops.txt"),
	read.delim("Nm1.5pops.2sample.n1000.sampledpops.txt"),
	read.delim("Nm1.5pops.4sample.n1000.sampledpops.txt"))
names(dat.list)<-c("Nm1.2pops.1sample.n1000",
	"Nm1.2pops.2sample.n1000","Nm1.2pops.4sample.n1000",
	"Nm1.5pops.1sample.n1000","Nm1.5pops.2sample.n1000",
	"Nm1.5pops.4sample.n1000")

png("fig5.png",width=169,height=169,units="mm",res=300)
par(mfrow = c(2, 3),cex = 0.6,mar = c(0, 0, 0, 0), 
	oma = c(4, 4, 1.5, 0.5), tcl = -0.25,mgp = c(2, 0.6, 0))
for (i in 1:length(dat.list)) {
	#plot(seq(0,0.6,0.06),seq(0,1,0.1), axes = FALSE, type = "n")
	
	y.max<-(max(dat.list[[i]]$WrightsFst)+2*sd(dat.list[[i]]$WrightsFst))
	if(y.max>1)
		y.max<-1
	plot(dat.list[[i]]$Ht,dat.list[[i]]$WrightsFst,las=1,
		col="black",pch=19,xaxt="n", xlim=c(0,0.6),ylim=c(0,y.max))
	Nm<-strsplit(strsplit(names(dat.list)[i], "Nm")[[1]][2],"[._]")[[1]][1]
	exp.fst<-round(1/((4*as.numeric(Nm))+1),digits=4)
	avg.fst<-round(mean(dat.list[[i]]$WrightsFst),digits=4)
	legend("topleft",c(paste("Exp. Fst = ",exp.fst),
		paste("Mean Fst = ",avg.fst)),bty="n")
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
mtext("Ht",1,outer=T,line=2,cex=.75)
mtext("Wright's Fst",2,outer=T,line=3,cex=.75)
dev.off()


############################################################################
#EXPLORING PARAMETER SPACE
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

############################################################################
#FIGS 6
############################################################################
#samples.ci<-lapply(strsplit(ci.files,".s"),function(x){
#	strsplit(x[2],".",fixed=TRUE)[[1]][1]})
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

nms<-c("Nm0.1.","Nm1.","Nm10.")
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

plot.lositan<-function(ci.list,loci.list,out.name){
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
			plot(sig.dat$Het, sig.dat$Fst,las=1,ylim=c(0,1),
				xlab="",ylab="",pch=19, xaxt="n",yaxt="n")
			points(ci.dat$Het,ci.dat[,2],col="yellow",type="l")
			points(ci.dat$Het, ci.dat[,4],col="red",type="l")
			axis(1)
			props<-find.los.sig(ci.dat,sig.dat)
			legend("topleft",bty="n",
				c(paste(round(props[1],3)," balancing"),
				paste(round(props[2],3)," positive"),
				paste(round(props[3],3)," total")))
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
			#if(i == 4){
			#	axis(1)}
			if(j == 1)
				axis(2,las=1)
		}
	}
	mtext("Heterozygosity",1,line=0.5,outer=T,cex=0.5)
	mtext("Fst",2,outer=T,line=2,cex=0.5)
	dev.off()
}

plot.lositan(s2.ci,s2.loci,"fig6.s2.png")
plot.lositan(s5.ci,s5.loci,"fig6.s5.png")
plot.lositan(s10.ci,s10.loci,"fig6.s10.png")
plot.lositan(s20.ci,s20.loci,"fig6.s20.png")

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
