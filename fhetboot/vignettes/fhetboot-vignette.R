## ------------------------------------------------------------------------
library("fhetboot",lib.loc = "B:/ubuntushare/fst_outliers/fhetboot")
gfile<-system.file("extdata", "example.genepop.txt",package = 'fhetboot')
gpop<-my.read.genepop(gfile)

## ------------------------------------------------------------------------
fsts<-calc.actual.fst(gpop)
head(fsts)

#Plot the actual values to see what your distribution looks like
par(mar=c(4,4,1,1))
plot(fsts$Ht, fsts$Fst,xlab="Ht",ylab="Fst",pch=19)

## ------------------------------------------------------------------------
boot.test<-fst.boot(gpop)
str(boot.test)
head(boot.test[[3]][[1]])

