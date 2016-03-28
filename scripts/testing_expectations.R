source("E:/ubuntushare/fst_outliers/scripts/bootstrap_fst_het.R")
setwd("E:/ubuntushare/fst_outliers/results/numerical_analysis_testingexp")

files<-list.files(pattern=".genepop$")

gpops<-lapply(files, my.read.genepop)
fsts<-lapply(gpops, calc.actual.fst)

#Exp Fst = 1/(4Nm+1)
exp<-c((1/((4*0.1)+1)),(1/(4+1)),(1/((4*10)+1)))

obs<-c(mean(as.numeric(fsts[[1]]$Fst)),mean(as.numeric(fsts[[2]]$Fst)),
	mean(as.numeric(fsts[[3]]$Fst)))

low.ci<-c(mean(as.numeric(fsts[[1]]$Fst))-(sd(as.numeric(fsts[[1]]$Fst))*1.96),
	mean(as.numeric(fsts[[2]]$Fst))-(sd(as.numeric(fsts[[2]]$Fst))*1.96),
	mean(as.numeric(fsts[[3]]$Fst))-(sd(as.numeric(fsts[[3]]$Fst))*1.96))

upp.ci<-c(mean(as.numeric(fsts[[1]]$Fst))+(sd(as.numeric(fsts[[1]]$Fst))*1.96),
	mean(as.numeric(fsts[[2]]$Fst))+(sd(as.numeric(fsts[[2]]$Fst))*1.96),
	mean(as.numeric(fsts[[3]]$Fst))+(sd(as.numeric(fsts[[3]]$Fst))*1.96))

nm<-c(0.1,1,10)

png("TestingExpectations.png",height=7,width=7,units="in",res=300)
par(lwd=1.3,cex=1.5)
plot(obs, xaxt='n',xlab="",ylab="",las=1,ylim=c(0,1),pch=19)
arrows(x1=1,y1=low.ci[1],x0=1,y0=upp.ci[1],code=3,angle=90)
arrows(x1=2,y1=low.ci[2],x0=2,y0=upp.ci[2],code=3,angle=90)
arrows(x1=3,y1=low.ci[3],x0=3,y0=upp.ci[3],code=3,angle=90)
points(exp,pch=8,col="blue")
axis(1,labels=c(0.1,1,10),at=c(1,2,3))
legend("top",ncol=2,c("Observed","Expected"),pch=c(19,8),
	col=c("black","blue"),bty='n')
mtext(expression(italic("Nm")),1,line=2,cex=1.5)
mtext(expression(italic("Fst")),2,line=2,cex=1.5)
dev.off()


