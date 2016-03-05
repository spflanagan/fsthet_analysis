
source("E:/ubuntushare/fst_outliers/scripts/bootstrap_fst_het.R")
setwd("E:/ubuntushare/fst_outliers/results/numerical_analysis_selection")
gpop<-my.read.genepop("Nm0.1.d2.s20.ds.genepop")

#gpop<-my.read.genepop("E:/ubuntushare/fst_outliers/results/numerical_analysis_genepop/Nm0.1.d2.s20.genepop")
#gpop<-my.read.genepop("Hess_2013_data_Genepop.gen") #this worked.
fsts<-calc.actual.fst(gpop)

boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
plot.cis(fsts,boot.out,make.file=F)

#example with many replicates
boot.out<-as.data.frame(t(replicate(100, fst.boot(gpop))))
plot.cis(fsts,boot.out)
boot1000.ci95<-mean.cis(boot.out[[2]])
boot1000.ci99<-mean.cis(boot.out[[3]])
cis<-as.data.frame(t(do.call(rbind,c(boot1000.ci95,boot1000.ci99))))
colnames(cis)<-c("low95","upp95","low99","upp99")
cis$Ht<-as.numeric(rownames(cis))
write.csv(cis,"Bootstrapping1000CIs.csv")
#recalculate bins
outliers1000<-find.outliers(fsts,boot.out=boot.out)
outliers100<-find.outliers(fsts,boot.out=bootres)
plot.cis(fsts,bootres)
#example with just one bootstrap rep
test<-fst.boot(gpop)
plot.cis(fsts,
	ci.list=list(test$CI95[,1],test$CI95[,2],test$CI99[,1],test$CI99[,2]))
ci.df<-data.frame(low95=test$CI95[,1],upp95=test$CI95[,2],
	low99=test$CI99[,1],upp99=test$CI99[,2])
ci.df$Ht<-as.numeric(rownames(ci.df))
write.csv(ci.df,"Bootstrap1_CIs.csv")
out1<-find.outliers(fsts,ci.df=ci.df)


 bootres<-as.data.frame(t(replicate(100,fst.boot(gpop))))