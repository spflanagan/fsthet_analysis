source("E:/ubuntushare/fst_outliers/scripts/bootstrap_fst_het.R")
setwd("E:/ubuntushare/fst_outliers/results/numerical_analysis_selection")

sel<-read.table("Nm0.1.d20.s20.ds0.01.sampledpops.txt",header=T)
sig<-read.table('Nm0.1.d20.s20.ds0.005.sigloci.txt')

gpop<-my.read.genepop("Nm0.1.d20.s20.ds0.005.genepop")
	fsts<-calc.actual.fst(gpop)
	boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop))))
	plot.cis(fsts,boot.out,make.file=F)
points(fsts[fsts$Locus %in% sig$V1,c("Ht","Fst")],col="red",pch=8)
outliers<-find.outliers(fsts,boot.out=boot.out, 
		file.name="Nm0.1.d20.s20.ds0.01.genepop")
prop95<-nrow(outliers[[1]])/(ncol(gpop)-2)
	prop99<-nrow(outliers[[2]])/(ncol(gpop)-2)