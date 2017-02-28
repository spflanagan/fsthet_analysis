source("B:/ubuntushare/fst_outliers/fhetboot/R/fhetboot.R")
setwd("B:/ubuntushare/fst_outliers/results/numerical_analysis_selection")

rm(list=ls())
sel.all.files<-list.files(pattern="d*.s20.ds.*.genepop$")
ds0.files<-sel.all.files[grep("ds0.genepop",sel.all.files)]

sel.proportions<-data.frame(wcc.prop.out=numeric(),wcc.prop.sig=numeric())

#do.call(rbind,lapply(sel.all.files, function(x) {
for(i in 1:length(sel.all.files)){
	gpop<-my.read.genepop(sel.all.files[i])
	fsts.wcc<-calc.actual.fst(gpop,"WCC")
#	fsts<-calc.actual.fst(gpop)
#	boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop))))
	wcc.boot.out<-as.data.frame(t(replicate(1,fst.boot(gpop,"WCC", bootstrap = FALSE))))
	plotting.cis(fsts.wcc,wcc.boot.out,make.file=T,file.name=paste(sel.all.files[i],"wcc.noboot.png",sep=""))
#	outliers<-find.outliers(fsts,boot.out=boot.out, 
#		file.name=sel.all.files[i])
	wcc.outliers<-find.outliers(fsts.wcc,boot.out=wcc.boot.out, 
		file.name=sel.all.files[i])
	wcc.outliers<-wcc.outliers[wcc.outliers$Ht != 0,]
	wcc.prop.out<-nrow(wcc.outliers)/(ncol(gpop)-2)
	if(sel.all.files[i] %in% ds0.files){
	  wcc.prop.Sig<-0
	} else {
  	sig<-read.table(gsub("genepop","sigloci.txt",sel.all.files[i]))
  	wcc.prop.Sig<-length(sig$V1[sig$V1 %in% wcc.outliers$Locus])/nrow(sig)
	}
	sel.proportions[i,]<-cbind(wcc.prop.out,wcc.prop.Sig)
	#return(cbind(wcc.prop.out,wcc.prop.sig))
}
rownames(sel.proportions)<-sel.all.files
sel.proportions$Demes<-as.numeric(gsub("Nm(\\d+.*).d(\\d+).s(\\d+).ds(\\d.*).genepop","\\2",rownames(sel.proportions)))
sel.proportions$Nm<-as.numeric(gsub("Nm(\\d+.*).d(\\d+).s(\\d+).ds(\\d.*).genepop","\\1",rownames(sel.proportions)))
sel.proportions$Selection<-as.numeric(gsub("Nm(\\d+.*).d(\\d+).s(\\d+).ds(\\d+.*).genepop","\\4",rownames(sel.proportions)))
props<-data.frame(#Selection0=sel.proportions[sel.proportions$Selection == 0,],
                  Selection0.01=sel.proportions[sel.proportions$Selection == 0.01,],
                  Selection0.1=sel.proportions[sel.proportions$Selection == 0.1,],
                  Selection0.5=sel.proportions[sel.proportions$Selection == 0.5,],
                  Selection1=sel.proportions[sel.proportions$Selection == 1,])
props.out<-data.frame(Demes=props$Selection0.01.Demes,Nm=props$Selection0.01.Nm,
                      props[,grep("wcc.prop",colnames(props))])

props.out<-props.out[order(props.out$Demes,props.out$Nm),]

write.table(props.out,"SelectedProportionOutliers.24.02.2017.txt",sep="\t",quote=F,row.names=F,
	col.names=T)

####If you already made the files:
out.files<-list.files(pattern=".genepop.csv")
ds0.files<-out.files[grep("ds0.genepop",out.files)]
out.files<-out.files[!(out.files %in% ds0.files)]
out.files<-out.files[grep("ds",out.files)]

proportions<-do.call(rbind,lapply(out.files, function(x) {
	prop.dat<-read.csv(x)
	prop<-nrow(prop.dat)/2000
	sig<-read.table(gsub("genepop.csv","sigloci.txt",x))
	propSig<-length(sig$V1[sig$V1 %in% prop.dat$Locus])/nrow(sig)
	return(cbind(prop,propSig))
}))
rownames(proportions)<-out.files
write.table(proportions,"ProportionOutliers_SelLoci.txt",sep="\t",quote=F,row.names=T,
	col.names=T)
proportions$Nm<-as.vector(gsub("Nm(.*).d\\d+.*","\\1",rownames(proportions)))
proportions$d<-as.vector(gsub(".*.d(\\d+).*","\\1",rownames(proportions)))
proportions$s<-as.vector(gsub("Nm.*d.*.s(\\d+)[:punct:]ds.*","\\1",rownames(proportions)))
proportions$ds<-as.vector(gsub(".*.ds(.*\\d).genepop95.csv","\\1",rownames(proportions)))
proportions$ds[props$ds=="Nm20_d20_s20_ds01_test.genepop95.csv"]<-"0.01"
props$s[props$d=="Nm20_d20_s20_ds01_test.genepop95.csv"]<-"0.01"
proportions$Nm<-factor(proportions$Nm)
proportions$d<-factor(proportions$d)
proportions$s<-factor(proportions$d)
proportions$ds<-factor(proportions$ds)

gsub("Nm.*s(\\d+)[ds.*","\\1",rownames(proportions))
###########TESTING
gpop<-my.read.genepop("Nm20_d20_s20_ds01_test.genepop")
sig<-read.table("Nm20_d20_s20_ds01_test.sigloci.txt",sep='\t')
fsts<-calc.actual.fst(gpop)
plot(fsts$Ht,fsts$Fst)
points(fsts[fsts$Locus %in% sig$V1,"Ht"],fsts[fsts$Locus %in% sig$V1,"Fst"],pch=8,col="red")
deltaq<-read.table("Nm20_d20_s20_ds01_test.delatq.txt",header=T)
head(deltaq) 
dqcalc<-function(x){
	m<-0.02
	s<-0.01
	h<-0.5
	dq<-((s*x$q)*(1-x$q)*(1-x$q+(h*((2*x$q)-1))))-(m*(x$q-x$Qbar))
}
rdq<-dqcalc(deltaq)
  colnames(deltaq)<-c("locus","pop","q","qhat","p","qbar","m","s")

deltas<-data.frame(dq=c(deltaq$Deltaq,deltaq$ObsDeltaq),
	qtype=c(rep("Exp",nrow(deltaq)),rep("Obs",nrow(deltaq))),
	locus=c(as.factor(as.character(deltaq$Locus)),
		as.factor(as.character(deltaq$Locus))))

qhat<-function(s,qbar,m){
	a<-s*-1*0.5
	b<-(s/2)-m
	c<-m*qbar
	q1<-((-1*b)+sqrt((b^2)-(4*a*c)))/(2*a)
	q2<-((-1*b)-sqrt((b^2)-(4*a*c)))/(2*a)
	return(cbind(q1,q2))
}

require(ggplot2)
png("exp_obs_dq.png",height=7,width=8,units="in",res=300)
par(lwd=1.3,cex=1.3)
p<-ggplot(data = deltas, aes(x=as.factor(locus), y=dq)) + 
	geom_boxplot(aes(fill=qtype))
p <- p + xlab("Locus") + ylab("Change in q")
p <- p + guides(fill=guide_legend(title="Type of delta-q"))
p
dev.off()

dq.split<-split(deltaq,deltaq$Locus)
png("checking_q.png",height=11,width=8,units="in",res=300)
par(mfrow=c(5,2),oma=c(2,2,2,2),mar=c(2,2,2,2),lwd=1.3,cex=0.5)
for(i in 1:length(dq.split)){
	qhat<-sig[sig$V1 %in% factor(dq.split[[i]]$Locus),"V3"]
	meanobs<-rowMeans(sig[sig$V1 %in% factor(dq.split[[i]]$Locus),4:23],
		na.rm=T)
	plot(dq.split[[i]]$Gen,dq.split[[i]]$q,xlab="",ylab="",type="l",las=1,
		ylim=c(0,1))
	mtext(dq.split[[i]]$Locus[1],3)
	if(meanobs>0.5){
		legend("bottom",c(paste("qbar=",dq.split[[i]]$Qbar[1]),
			paste("qhat=",qhat),paste("Mean q at gen 5000=",meanobs)),
			bty="n")
	} else {
		legend("top",c(paste("qbar=",dq.split[[i]]$Qbar[1]),
			paste("qhat=",qhat),paste("Mean q at gen 5000=",meanobs)),
			bty="n")
	} 
}
mtext("Generation",1,outer=T)
mtext("q",2,outer=T)
dev.off()

#LOSITAN FILES

iam.ci<-list.files(pattern="genepop.ci")
iam.loci<-list.files(pattern="genepop.loci")
stp.ci<-list.files(pattern="step.ci")
stp.loci<-list.files(pattern="step.loci")
ss.ci<-data.frame(filename=character(),ParamSet=character(),
       Nm=numeric(),demes=numeric(),sampled=numeric(), 
	PropBal=numeric(),PropPos=numeric(),PropOutliers=numeric(),
	stringsAsFactors=F)

for(i in 1:length(stp.ci)){
  nm<-gsub("Nm(\\d.*).d.*","\\1",stp.ci[i])
  d<-gsub("Nm\\d.*.d(\\d+).*","\\1",stp.ci[i])
  s<-gsub("Nm\\d.*.d\\d+.*.s(\\d+).*","\\1",stp.ci[i])
  params<-gsub("(.*).genepop.*","\\1",stp.ci[i])
  dat<-read.delim(stp.ci[i])
  low.ci<-dat[,c(1,2)]
  upp.ci<-dat[,c(1,4)]
  locus.name<-paste(gsub("(Nm\\d.*.genepop.step).ci","\\1",stp.ci[i]),"loci",sep=".")
  loc.name<-stp.loci[stp.loci %in% locus.name]
  loc.dat<-read.delim(loc.name)
  loc.dat<-loc.dat[loc.dat$Fst>0,]
  out<-data.frame(t(apply(loc.dat,1,function(x){
	pos<-0
	bal<-0
    #get the two confidence interval values closest to x
	low.fst<-low.ci[which.min(abs(as.numeric(low.ci$Het)-as.numeric(x["Het"]))),]
	upp.fst<-upp.ci[which.min(abs(as.numeric(upp.ci$Het)-as.numeric(x["Het"]))),]
   if(as.numeric(x["Fst"]) > as.numeric(upp.fst[2])){
	pos<-1
	}
	if(as.numeric(x["Fst"]) < as.numeric(low.fst[2])){
      	bal<-1
    }
    return(cbind(bal,pos))
  })))
  prop.bal<-sum(out$X1)/nrow(loc.dat)
  prop.pos<-sum(out$X2)/nrow(loc.dat)
	prop.out<-sum(out$X1,out$X2)/nrow(loc.dat)
  ss.ci[i,]<-cbind(locus.name,params,nm,d,s,prop.bal,prop.pos,prop.out)
}

write.csv(ss.ci,"StepwiseLositanOutliersSelection.csv")

si.ci<-data.frame(filename=character(),ParamSet=character(),
                 Nm=numeric(),demes=numeric(),s=numeric(), PropOutliers=numeric(),stringsAsFactors=F)
for(i in 1:length(iam.ci)){
  nm<-gsub("Nm(\\d.*).d.*","\\1",iam.ci[i])
  d<-gsub("Nm\\d.*.d(\\d+).*","\\1",iam.ci[i])
  s<-gsub("Nm\\d.*.d\\d+.*.s(\\d+).*","\\1",iam.ci[i])
  params<-gsub("(.*).genepop.*","\\1",iam.ci[i])
  dat<-read.delim(iam.ci[i])
  low.ci<-dat[,c(1,2)]
  upp.ci<-dat[,c(1,4)]
  locus.name<-paste(gsub("(Nm\\d.*.genepop).ci","\\1",iam.ci[1]),"loci",sep=".")
  loc.name<-iam.loci[iam.loci %in% locus.name]
  loc.dat<-read.delim(loc.name)
  loc.dat<-loc.dat[loc.dat$Fst>-100,]
  out<-apply(loc.dat,1,function(x){
    outliers<-0
    #get the two confidence interval values closest to x
    lowl<-low.ci[as.numeric(low.ci$Het) <= as.numeric(x["Het"]),]
    lowl<-lowl[nrow(lowl),]
    uppl<-low.ci[as.numeric(low.ci$Het) >= as.numeric(x["Het"]),][1,]

    if(nrow(lowl) > 0 & nrow(uppl) > 0) { low.fst<-mean(lowl$Fst,uppl$Fst) 
	} else{
		if(nrow(lowl) > 0){ low.fst<-lowl$Fst }
		if(nrow(uppl) > 0){ low.fst<-uppl$Fst } 
	}
    lowu<-upp.ci[as.numeric(upp.ci$Het) <= as.numeric(x["Het"]),]
    lowu<-lowu[nrow(lowu),]
    uppu<-upp.ci[as.numeric(upp.ci$Het) >= as.numeric(x["Het"]),][1,]
    if(nrow(lowu) > 0 & nrow(uppu) > 0) {  upp.fst<-mean(lowu[,2],uppu[,2]) 
	} else{
		if(nrow(lowu) > 0){ upp.fst<-lowu[,2] }
		if(nrow(uppu) > 0){ upp.fst<-uppu[,2] } 
	}

    if(as.numeric(x["Fst"]) > as.numeric(upp.fst) | as.numeric(x["Fst"]) < as.numeric(low.fst)){
      outliers<-1
    }
    return(outliers)
  })
  outliers<-sum(out)/nrow(loc.dat)
  si.ci[i,]<-cbind(locus.name,params,nm,d,s,outliers)
}

sis.ci<-merge(s.ci,si.ci,by="ParamSet")
t.test(as.numeric(sis.ci$PropOutliers.x),as.numeric(sis.ci$PropOutliers.y),paired=T)

#t = 1.8257, df = 83, p-value = 0.07149


write.csv(si.ci,"InfiniteAllelesModel_Selection.csv")

#compare to lositan analysis
ss.ci<-read.csv("StepwiseLositanOutliers_Selection.csv")
ss.ci$filename<-gsub("(Nm\\d+.*.genepop).step.loci","\\1",ss.ci$filename)
si.ci<-read.csv("InfiniteAllelesModel_Selection.csv")
si.ci$filename<-gsub("(Nm\\d+.*.genepop).loci","\\1",si.ci$filename)
sel.proportions<-read.table("SelectedProportionOutliers.24.02.2017.txt",header=T)


props<-data.frame(Demes=rep(sel.proportions$Demes,4),Nm=rep(sel.proportions$Nm,4),
                  Selection=c(rep(0.01,nrow(sel.proportions)),rep(0.1,nrow(sel.proportions)),
                            rep(0.5,nrow(sel.proportions)),rep(1,nrow(sel.proportions))),
                  wcc.out=c(sel.proportions$Selection0.01.wcc.prop.out,sel.proportions$Selection0.1.wcc.prop.out,
                            sel.proportions$Selection0.5.wcc.prop.out,sel.proportions$Selection1.wcc.prop.out))
props$filename<-paste("Nm",props$Nm,".d",props$Demes,".s20.ds",props$Selection,".genepop",sep="")


#step
step.prop<-merge(ss.ci,props,by="filename")
t.test(step.prop$PropOutliers,step.prop$wcc.out,"greater")
#t = 76.727, df = 18, p-value < 2.2e-16
