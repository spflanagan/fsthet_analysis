source("E:/ubuntushare/fst_outliers/scripts/bootstrap_fst_het.R")
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