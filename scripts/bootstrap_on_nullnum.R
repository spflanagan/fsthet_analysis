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

##Analyze lositan output
setwd("~/Desktop/numerical_analysis_genepop/")
iam.ci<-list.files(pattern="genepop.ci")
iam.loci<-list.files(pattern="genepop.loci")
stp.ci<-list.files(pattern="step.ci")
stp.loci<-list.files(pattern="step.loci")
s.ci<-data.frame(filename=character(),ParamSet=character(),
                 Nm=numeric(),demes=numeric(),selection=numeric(), PropOutliers=numeric(),stringsAsFactors=F)
for(i in 1:length(stp.ci)){
  nm<-gsub("Nm(\\d.*).d.*","\\1",stp.ci[i])
  d<-gsub("Nm\\d.*.d(\\d+).*","\\1",stp.ci[i])
  s<-gsub("Nm\\d.*.d\\d+.*.s(\\d+).*","\\1",stp.ci[i])
  params<-gsub("(.*).genepop.*","\\1",stp.ci[i])
  dat<-read.delim(stp.ci[i])
  low.ci<-dat[,c(1,2)]
  upp.ci<-dat[,c(1,4)]
  locus.name<-paste(gsub("(Nm\\d.*.genepop.step).ci","\\1",stp.ci[1]),"loci",sep=".")
  loc.name<-stp.loci[stp.loci %in% locus.name]
  loc.dat<-read.delim(loc.name)
  loc.dat<-loc.dat[loc.dat$Fst>-100,]
  out<-apply(loc.dat,1,function(x){
    outliers<-0
    #get the two confidence interval values closest to x
    lowl<-low.ci[as.numeric(low.ci$Het) <= as.numeric(x["Het"]),]
    lowl<-lowl[nrow(lowl),]
    uppl<-low.ci[as.numeric(low.ci$Het) >= as.numeric(x["Het"]),][1,]
    low.fst<-mean(lowl$Fst,uppl$Fst)
    lowu<-upp.ci[as.numeric(upp.ci$Het) <= as.numeric(x["Het"]),]
    lowu<-lowu[nrow(lowu),]
    uppu<-upp.ci[as.numeric(upp.ci$Het) >= as.numeric(x["Het"]),][1,]
    upp.fst<-mean(lowu[,2],uppu[,2])
    if(as.numeric(x["Fst"]) > as.numeric(upp.fst) | as.numeric(x["Fst"]) < as.numeric(low.fst)){
      outliers<-1
    }
    return(outliers)
  })
  outliers<-sum(out)/nrow(loc.dat)
  s.ci[i,]<-cbind(locus.name,params,nm,d,s,outliers)
}


i.ci<-data.frame(filename=character(),ParamSet=character(),
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
    low.fst<-mean(lowl$Fst,uppl$Fst)
    lowu<-upp.ci[as.numeric(upp.ci$Het) <= as.numeric(x["Het"]),]
    lowu<-lowu[nrow(lowu),]
    uppu<-upp.ci[as.numeric(upp.ci$Het) >= as.numeric(x["Het"]),][1,]
    upp.fst<-mean(lowu[,2],uppu[,2])
    if(as.numeric(x["Fst"]) > as.numeric(upp.fst) | as.numeric(x["Fst"]) < as.numeric(low.fst)){
      outliers<-1
    }
    return(outliers)
  })
  outliers<-sum(out)/nrow(loc.dat)
  i.ci[i,]<-cbind(locus.name,params,nm,d,s,outliers)
}

is.ci<-merge(s.ci,i.ci,by="ParamSet")
t.test(as.numeric(is.ci$PropOutliers.x),as.numeric(is.ci$PropOutliers.y),paired=T)

#Paired t-test
#
#data:  as.numeric(is.ci$PropOutliers.x) and as.numeric(is.ci$PropOutliers.y)
#t = 2.6616, df = 83, p-value = 0.009335
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.001338050 0.009250578
#sample estimates:
#  mean of the differences 
#0.005294314 

#selection
setwd("~/Desktop/numerical_analysis_selection/")
iam.ci<-list.files(pattern="genepop.ci")
iam.loci<-list.files(pattern="genepop.loci")
stp.ci<-list.files(pattern="step.ci")
stp.loci<-list.files(pattern="step.loci")
ss.ci<-data.frame(filename=character(),ParamSet=character(),
                 Nm=numeric(),demes=numeric(),selection=numeric(), ds=numeric(),PropOutliers=numeric(),stringsAsFactors=F)
for(i in 1:length(stp.ci)){
  nm<-gsub("Nm(\\d.*).d.*","\\1",stp.ci[i])
  d<-gsub("Nm\\d.*.d(\\d+).*","\\1",stp.ci[i])
  s<-gsub("Nm\\d.*.d\\d+.*.s(\\d+).*","\\1",stp.ci[i])
  ds<-gsub("Nm\\d.*.d\\d+.*.s\\d+.ds(\\d.*).genepop.*","\\1",stp.ci[i])
  params<-gsub("(.*).genepop.*","\\1",stp.ci[i])
  dat<-read.delim(stp.ci[i])
  low.ci<-dat[,c(1,2)]
  upp.ci<-dat[,c(1,4)]
  locus.name<-paste(gsub("(Nm\\d.*.genepop.step).ci","\\1",stp.ci[1]),"loci",sep=".")
  loc.name<-stp.loci[stp.loci %in% locus.name]
  loc.dat<-read.delim(loc.name)
  loc.dat<-loc.dat[loc.dat$Fst>-100,]
  out<-apply(loc.dat,1,function(x){
    outliers<-0
    #get the two confidence interval values closest to x
    lowl<-low.ci[as.numeric(low.ci$Het) <= as.numeric(x["Het"]),]
    lowl<-lowl[nrow(lowl),]
    uppl<-low.ci[as.numeric(low.ci$Het) >= as.numeric(x["Het"]),][1,]
    low.fst<-mean(lowl$Fst,uppl$Fst)
    lowu<-upp.ci[as.numeric(upp.ci$Het) <= as.numeric(x["Het"]),]
    lowu<-lowu[nrow(lowu),]
    uppu<-upp.ci[as.numeric(upp.ci$Het) >= as.numeric(x["Het"]),][1,]
    upp.fst<-mean(lowu[,2],uppu[,2])
    if(as.numeric(x["Fst"]) > as.numeric(upp.fst) | as.numeric(x["Fst"]) < as.numeric(low.fst)){
      outliers<-1
    }
    return(outliers)
  })
  outliers<-sum(out)/nrow(loc.dat)
  ss.ci[i,]<-cbind(locus.name,params,nm,d,s,ds,outliers)
}


si.ci<-data.frame(filename=character(),ParamSet=character(),
                 Nm=numeric(),demes=numeric(),s=numeric(), ds=numeric(),PropOutliers=numeric(),stringsAsFactors=F)
for(i in 1:length(iam.ci)){
  nm<-gsub("Nm(\\d.*).d.*","\\1",iam.ci[i])
  d<-gsub("Nm\\d.*.d(\\d+).*","\\1",iam.ci[i])
  s<-gsub("Nm\\d.*.d\\d+.*.s(\\d+).*","\\1",iam.ci[i])
  ds<-gsub("Nm\\d.*.d\\d+.*.s\\d+.ds(\\d.*).genepop.*","\\1",iam.ci[i])
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
    low.fst<-mean(lowl$Fst,uppl$Fst)
    lowu<-upp.ci[as.numeric(upp.ci$Het) <= as.numeric(x["Het"]),]
    lowu<-lowu[nrow(lowu),]
    uppu<-upp.ci[as.numeric(upp.ci$Het) >= as.numeric(x["Het"]),][1,]
    upp.fst<-mean(lowu[,2],uppu[,2])
    if(as.numeric(x["Fst"]) > as.numeric(upp.fst) | as.numeric(x["Fst"]) < as.numeric(low.fst)){
      outliers<-1
    }
    return(outliers)
  })
  outliers<-sum(out)/nrow(loc.dat)
  si.ci[i,]<-cbind(locus.name,params,nm,d,s,ds,outliers)
}

sis.ci<-merge(ss.ci,si.ci,by="ParamSet")
t.test(as.numeric(sis.ci$PropOutliers.x),as.numeric(sis.ci$PropOutliers.y),paired=T)

#Paired t-test
#
#data:  as.numeric(sis.ci$PropOutliers.x) and as.numeric(sis.ci$PropOutliers.y)
#t = 83.581, df = 7, p-value = 9.24e-12
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.7782770 0.8235961
#sample estimates:
#  mean of the differences 
#0.8009366 