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

ci.boot<-ci.means(boot.out[[3]])
boot.means<-fst.boot.means(boot.out)

#example with just one bootstrap rep
test<-fst.boot(gpop)
plotting.cis(fsts, ci.list=test[[3]])
out1<-find.outliers(fsts,ci.df=ci.df)

#example with 100 replicates
boot100<-as.data.frame(t(replicate(100, fst.boot(gpop))))
plotting.cis(fsts,boot100)
boot100.ci95<-ci.means(boot100[[3]])
cis100<-as.data.frame(t(do.call(rbind,boot100.ci95)))
colnames(cis100)<-c("low95","upp95")
cis100$Ht<-as.numeric(rownames(cis100))
write.csv(cis100,"Bootstrapping100CIs.csv")

#example with 100 replicates
boot1000<-as.data.frame(t(replicate(1000, fst.boot(gpop))))
plotting.cis(fsts,boot1000)
boot1000.ci95<-ci.means(boot1000[[3]])
cis1000<-as.data.frame(t(do.call(rbind,boot1000.ci95)))
colnames(cis1000)<-c("low95","upp95")
cis1000$Ht<-as.numeric(rownames(cis1000))
write.csv(cis1000,"Bootstrapping1000CIs.csv")

#recalculate bins
outliers1000<-find.outliers(fsts,boot.out=boot1000)
outliers100<-find.outliers(fsts,boot.out=boot100)
plot.cis(fsts,boot100)



