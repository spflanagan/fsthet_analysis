source("E:/ubuntushare/fst_outliers/scripts/bootstrap_fst_het.R")
setwd("B:/ubuntushare/fst_outliers/results/numerical_analysis_genepop")

all.files<-list.files(pattern=".genepop$")
proportions<-do.call(rbind,lapply(all.files, function(x) {
	gpop<-my.read.genepop(x)
	fsts.wcc<-calc.actual.fst(gpop,"WCC")
	fsts<-calc.actual.fst(gpop)
	boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop))))
	wcc.boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop,"WCC"))))
	plotting.cis(fsts,boot.out,make.file=T,file.name=paste(x,"wcc.png",sep=""))
	outliers<-find.outliers(fsts,boot.out=boot.out, 
		file.name=x)
	wcc.outliers<-find.outliers(fsts,boot.out=wcc.boot.out, 
		file.name=x)
	wcc.prop<-nrow(wcc.outliers)/(ncol(gpop)-2)
	prop<-nrow(outliers)/(ncol(gpop)-2)
	return(cbind(prop,wcc.prop))
}))
rownames(proportions)<-all.files
write.table(proportions,"ProportionOutliers_WCC.txt",sep="\t",quote=F,row.names=T,
	col.names=T)

#compare to lositan analysis?
los.sig<-read.table("B:/ubuntushare/fst_outliers/results/numerical_analysis_genepop/sig.loci.sim.LOSITAN.txt",skip=1)
colnames(los.sig)<-c("file","bal","pos","total","balp","posp","totalp")
#these are all at the 95% level
#do a paired t-test.
los.sig<-los.sig[order(los.sig$file),]
proportions<-proportions[order(rownames(proportions)),]
t.test(los.sig$totalp, proportions[,1], paired = TRUE)
t.test(los.sig$totalp,proportions[,1],paired=T,alternative="greater")

max.boot<-read.csv("Nm10.d50.s2.genepop95.csv")