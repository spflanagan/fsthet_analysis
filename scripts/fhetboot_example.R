#Example use of fhetboot
#Author: Sarah P. Flanagan
#Last updated: 14 July 2016


library("fhetboot",lib.loc = "E:/ubuntushare/fst_outliers/fhetboot")
gfile<-system.file("extdata", "example.genepop.txt",package = 'fhetboot')
gpop<-my.read.genepop(gfile)

#basic analysis
fsts<-calc.actual.fst(gpop)
boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
plotting.cis(fsts,boot.out,make.file=F)

#example with just one bootstrap rep
test<-fst.boot(gpop)
plotting.cis(fsts,
             ci.list=list(test$CI95[,1],test$CI95[,2],test$CI99[,1],test$CI99[,2]))
ci.df<-data.frame(low95=test$CI95[,1],upp95=test$CI95[,2],
                  low99=test$CI99[,1],upp99=test$CI99[,2])
ci.df$Ht<-as.numeric(rownames(ci.df))
write.csv(ci.df,"Bootstrap1_CIs.csv")
out1<-find.outliers(fsts,ci.df=ci.df)

#example with 100 replicates
boot100<-as.data.frame(t(replicate(100, fst.boot(gpop))))
plotting.cis(fsts,boot100)
boot100.ci95<-ci.means(boot100[[2]])
boot100.ci99<-ci.means(boot100[[3]])
cis100<-as.data.frame(t(do.call(rbind,c(boot100.ci95,boot100.ci99))))
colnames(cis100)<-c("low95","upp95","low99","upp99")
cis100$Ht<-as.numeric(rownames(cis100))
write.csv(cis100,"Bootstrapping100CIs.csv")

#example with 100 replicates
boot1000<-as.data.frame(t(replicate(1000, fst.boot(gpop))))
plotting.cis(fsts,boot1000)
boot1000.ci95<-ci.means(boot1000[[2]])
boot1000.ci99<-ci.means(boot1000[[3]])
cis1000<-as.data.frame(t(do.call(rbind,c(boot1000.ci95,boot1000.ci99))))
colnames(cis1000)<-c("low95","upp95","low99","upp99")
cis1000$Ht<-as.numeric(rownames(cis1000))
write.csv(cis1000,"Bootstrapping1000CIs.csv")

#recalculate bins
outliers1000<-find.outliers(fsts,boot.out=boot1000)
outliers100<-find.outliers(fsts,boot.out=boot100)
plot.cis(fsts,boot100)



