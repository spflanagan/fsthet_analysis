library(fsthet)

####TESTING####
setwd("~/Projects/fst_outliers/results/numerical_analysis_selection/")
gpop<-my.read.genepop("Nm1.d2.s20.ds0.genepop")#"Nm0.1.d2.s2.genepop.step.loci"
f1<-calc.actual.fst(gpop,"fst")
f2<-calc.actual.fst(gpop,"var")
f3<-calc.actual.fst(gpop,"theta")
f4<-calc.actual.fst(gpop,"betahat")
f<-calc.actual.fst(gpop,"wcc") #run with corrected weir & cockerham
f.n<-calc.actual.fst(gpop,"wright") #run with nei's formulation of wright's fst
f.w<-calc.actual.fst(gpop,"wc") #run with weir and cockerham uncorrected
f.fsthet<-as.data.frame(t(replicate(1,fst.boot(gpop,"wcc",bootstrap = F)))) #not bootstrapping
f.fhetboot<-as.data.frame(t(replicate(10,fst.boot(gpop,"wcc",bootstrap=T)))) #bootstrapping
los<-read.delim("Nm1.d2.s20.ds0.genepop.step.loci") #the lositan data
plot(los$Het,los$Fst,ylim=c(0,1),pch=19)
points(f1$Ht,f1$Fst,pch=1,col="orchid")
points(f2$Ht,f2$Fst,pch=2,col="cornflowerblue")
points(f3$Ht,f3$Fst,pch=3,col="forestgreen")
points(f4$Ht,f4$Fst,pch=4,col="goldenrod")
gpop<-my.read.genepop("~/Desktop/genepopdata.txt")
adams<-read.delim("~/Desktop/Fdistcheckerresults.xls")
points(adams$Ht,adams$Fst,col="blue")
points(adams$Ht,adams$LosFst,col="purple")
points(adams$Ht,adams$theta,col="grey")
points(adams$popHt,adams$popFst,col="goldenrod")
points(adams$HB,adams$b.hat,col="green")
points(f$Ht,f$Fst,col="cornflowerblue",pch=5)
points(f.fsthet$Fsts[[1]]$Ht,f.fsthet$Fsts[[1]]$Fst,col="orchid",pch=3)

####PLOT THE CORRELATION####
setwd("../numerical_analysis_genepop/")
nm1.gpop<-list.files(pattern="Nm1\\..*.genepop$")
png("supplemental_fig2.png",height=8,width=8,units="in",res=300)
pdf("supplemental_fig2.pdf")
par(mfrow=c(7,4),mar=c(1,1,1,1),oma=c(2.5,2.5,1,1),cex=0.5)
lapply(nm1.gpop,function(gpop.file){
  gpop<-my.read.genepop(gpop.file)
  f1<-calc.actual.fst(gpop,"fst")
  f4<-calc.actual.fst(gpop,"betahat")
  plot(f1$Ht,f4$Ht,pch=19,bty="L",cex.axis=1.25)
  r<-cor.test(f1$Ht,f4$Ht)$estimate
  p<-cor.test(f1$Ht,f4$Ht)$p.value
  legend("topleft",bty='n',cex=1.25,
         legend=c(paste("p = ", round(p,3),sep=""),paste("r = ",round(r,3),sep="")))
  legend(legend=sub("(.*).genepop","\\1",gpop.file),bty='n',"bottomright",cex=1.25)
})
  mtext(expression(italic(H)[T]),1,cex=0.75,outer=T,line=1.5)
  mtext(expression(italic(H)[B]),2,cex=0.75,outer=T,line=1)
dev.off()

