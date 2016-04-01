source("E:/ubuntushare/fst_outliers/fhetboot/R/fhetboot.R")
setwd("E:/ubuntushare/fst_outliers/results/numerical_analysis_selection")

all.files<-list.files(pattern=".genepop$")

proportions<-do.call(rbind,lapply(all.files, function(x) {
	gpop<-my.read.genepop(x)
	sig<-read.table(paste(gsub("(.*).genepop","\\1",all.files[1]),
		"sigloci.txt",sep="."))
	fsts<-calc.actual.fst(gpop)
	boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop))))
	png(paste(x,".png",sep=""),height=7,width=7,units="in",res=300)
	plot.cis(fsts,boot.out,make.file=F)
	points(fsts[fsts$Locus %in% sig$V1,c("Ht","Fst")],col="red",pch=8)
	dev.off()
	outliers<-find.outliers(fsts,boot.out=boot.out, 
		file.name=x)
	prop95<-nrow(outliers[[1]])/(ncol(gpop)-2)
	prop99<-nrow(outliers[[2]])/(ncol(gpop)-2)
	propSig95<-length(sig$V1[sig$V1 %in% outliers[[1]]$Locus])/nrow(sig)
	propSig99<-length(sig$V1[sig$V1 %in% outliers[[2]]$Locus])/nrow(sig)
	return(cbind(prop95,prop99,propSig95,propSig99))
}))
rownames(proportions)<-all.files
write.table(proportions,"ProportionOutliers.txt",sep="\t",quote=F,row.names=T,
	col.names=T)

####If you already made the files:
all.files<-list.files(pattern=".genepop95")

proportions<-do.call(rbind,lapply(all.files, function(x) {
	prop95.dat<-read.csv(x,row.names=1)
	prop99.dat<-read.csv(gsub("95","99",x),row.names=1)
	prop95<-nrow(prop95.dat)/2000
	prop99<-nrow(prop99.dat)/2000
	sig<-read.table(gsub("genepop95.csv","sigloci.txt",x))
	propSig95<-length(sig$V1[sig$V1 %in% prop95.dat$Locus])/nrow(sig)
	propSig99<-length(sig$V1[sig$V1 %in% prop99.dat$Locus])/nrow(sig)
	return(cbind(prop95,prop99,propSig95,propSig99))
}))
rownames(proportions)<-all.files
write.table(proportions,"ProportionOutliers.txt",sep="\t",quote=F,row.names=T,
	col.names=T)
proportions$Nm<-as.vector(gsub("Nm(.*).d\\d+.*","\\1",rownames(proportions)))
proportions$d<-as.vector(gsub(".*.d(\\d+).*","\\1",rownames(proportions)))
proportions$s<-as.vector(gsub("Nm.*d.*.s(\\d+)[:punct:]ds.*","\\1",rownames(proportions)))
proportions$ds<-as.vector(gsub(".*.ds(.*\\d).genepop95.csv","\\1",rownames(proportions)))
proportions$ds[props$ds=="Nm20_d20_s20_ds01_test.genepop95.csv"]<-"0.01"
props$s[props$d=="Nm20_d20_s20_ds01_test.genepop95.csv"]<-"0.01"
proportions$Nm<-factor(proportions$Nm)
proportions$d<-factor(proportions$d)
proportions$s<-factor(proportions$d)
proportions$ds<-factor(proportions$ds)

gsub("Nm.*s(\\d+)[ds.*","\\1",rownames(proportions))
###########TESTING
gpop<-my.read.genepop("Nm20_d20_s20_ds01_test.genepop")
sig<-read.table("Nm20_d20_s20_ds01_test.sigloci.txt",sep='\t')
fsts<-calc.actual.fst(gpop)
plot(fsts$Ht,fsts$Fst)
points(fsts[fsts$Locus %in% sig$V1,"Ht"],fsts[fsts$Locus %in% sig$V1,"Fst"],pch=8,col="red")
deltaq<-read.table("Nm20_d20_s20_ds01_test.delatq.txt",header=T)
head(deltaq) 
dqcalc<-function(x){
	m<-0.02
	s<-0.01
	h<-0.5
	dq<-((s*x$q)*(1-x$q)*(1-x$q+(h*((2*x$q)-1))))-(m*(x$q-x$Qbar))
}
rdq<-dqcalc(deltaq)
  colnames(deltaq)<-c("locus","pop","q","qhat","p","qbar","m","s")

deltas<-data.frame(dq=c(deltaq$Deltaq,deltaq$ObsDeltaq),
	qtype=c(rep("Exp",nrow(deltaq)),rep("Obs",nrow(deltaq))),
	locus=c(as.factor(as.character(deltaq$Locus)),
		as.factor(as.character(deltaq$Locus))))

qhat<-function(s,qbar,m){
	a<-s*-1*0.5
	b<-(s/2)-m
	c<-m*qbar
	q1<-((-1*b)+sqrt((b^2)-(4*a*c)))/(2*a)
	q2<-((-1*b)-sqrt((b^2)-(4*a*c)))/(2*a)
	return(cbind(q1,q2))
}

require(ggplot2)
png("exp_obs_dq.png",height=7,width=8,units="in",res=300)
par(lwd=1.3,cex=1.3)
p<-ggplot(data = deltas, aes(x=as.factor(locus), y=dq)) + 
	geom_boxplot(aes(fill=qtype))
p <- p + xlab("Locus") + ylab("Change in q")
p <- p + guides(fill=guide_legend(title="Type of delta-q"))
p
dev.off()

dq.split<-split(deltaq,deltaq$Locus)
png("checking_q.png",height=11,width=8,units="in",res=300)
par(mfrow=c(5,2),oma=c(2,2,2,2),mar=c(2,2,2,2),lwd=1.3,cex=0.5)
for(i in 1:length(dq.split)){
	qhat<-sig[sig$V1 %in% factor(dq.split[[i]]$Locus),"V3"]
	meanobs<-rowMeans(sig[sig$V1 %in% factor(dq.split[[i]]$Locus),4:23],
		na.rm=T)
	plot(dq.split[[i]]$Gen,dq.split[[i]]$q,xlab="",ylab="",type="l",las=1,
		ylim=c(0,1))
	mtext(dq.split[[i]]$Locus[1],3)
	if(meanobs>0.5){
		legend("bottom",c(paste("qbar=",dq.split[[i]]$Qbar[1]),
			paste("qhat=",qhat),paste("Mean q at gen 5000=",meanobs)),
			bty="n")
	} else {
		legend("top",c(paste("qbar=",dq.split[[i]]$Qbar[1]),
			paste("qhat=",qhat),paste("Mean q at gen 5000=",meanobs)),
			bty="n")
	} 
}
mtext("Generation",1,outer=T)
mtext("q",2,outer=T)
dev.off()