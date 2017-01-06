#Example use of fhetboot
#Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
#Last updated: 5 December 2016


library("fhetboot",lib.loc = "B:/ubuntushare/fst_outliers/fhetboot")
gfile<-system.file("extdata", "example.genepop.txt",package = 'fhetboot')
gpop<-my.read.genepop(gfile)

#basic analysis
fsts<-calc.actual.fst(gpop)
fsts.wc<-calc.actual.fst(gpop,"wc")
fsts.wcc<-calc.actual.fst(gpop,"wcc")

#compare
plot(fsts$Ht,fsts$Fst,col="black",pch=19,ylim=c(-0.05,1))
points(fsts.wc$Ht,fsts.wc$Fst,col="blue",cex=1.3)
points(fsts.wcc$Ht,fsts.wcc$Fst,pch=5,col="red")

boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
boot.pvals<-p.boot(fsts,boot.out=boot.out)
boot.cor.pvals<-p.adjust(boot.pvals,method="BH")
boot.sig<-boot.cor.pvals[boot.cor.pvals <= 0.05]
plotting.cis(fsts,boot.out,make.file=F,sig.list=names(boot.sig),smooth.ci = TRUE,smoothing.rate = 0.02)

#example with just one bootstrap rep
test<-fst.boot(gpop)
plotting.cis(fsts, ci.list=test[[3]])
out1<-find.outliers(fsts,ci.df=ci.df)

#running it without bootstrapping
non.boot.out<-as.data.frame(t(replicate(1, fst.boot(gpop,bootstrap = FALSE))))
non.boot.pvals<-p.boot(fsts.wcc,boot.out=non.boot.out)
non.boot.cor.pvals<-p.adjust(non.boot.pvals,method="BH")
non.boot.sig<-non.boot.cor.pvals[non.boot.cor.pvals <= 0.05]
non.boot.outliers<-find.outliers(fsts.wcc,non.boot.out)
plotting.cis(fsts.wcc,non.boot.out,make.file=F,sig.list=names(non.boot.sig))
