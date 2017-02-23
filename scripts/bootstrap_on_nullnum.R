library("fhetboot",lib.loc = "~/Projects/fst_outliers/fhetboot")
setwd("~/Projects/fst_outliers/results/numerical_analysis_genepop")

find.los.sig<-function(ci.dat,loci.dat){
  bal<-NULL
  pos<-NULL
  loci.dat<-loci.dat[loci.dat$Het>0,]
  for(ii in 1:(nrow(ci.dat)-1)){
    fst.het<-loci.dat[loci.dat$Het >= ci.dat[ii,"Het"] & 
                        loci.dat$Het <= ci.dat[(ii+1),"Het"],]
    max.bound<-min(ci.dat[ii:(ii+1),4])
    min.bound<-max(ci.dat[ii:(ii+1),2])
    bal<-rbind(bal, fst.het[fst.het$Fst <= min.bound,]) 
    pos<-rbind(pos, fst.het[fst.het$Fst >= max.bound,])
  }
  props<-c((nrow(bal)/nrow(loci.dat)), (nrow(pos)/nrow(loci.dat)),
           ((nrow(bal)+nrow(pos))/nrow(loci.dat)))	
  return(props)
}


all.files<-list.files(pattern=".genepop$")
proportions<-data.frame(wcc.prop=numeric())

#proportions<-do.call(rbind,lapply(all.files, function(x) {
for(i in 1:length(all.files)){
	print(all.files[i])
	gpop<-my.read.genepop(all.files[i])
	fsts.wcc<-calc.actual.fst(gpop,"WCC")
	#fsts<-calc.actual.fst(gpop)
	#boot.out<-as.data.frame(t(replicate(10,fst.boot(gpop))))
	wcc.boot.out<-as.data.frame(t(replicate(1,fst.boot(gpop,"WCC", bootstrap = FALSE))))
	print("plotting.cis")
	plotting.cis(fsts.wcc,wcc.boot.out,make.file=T,file.name=paste(all.files[i],"wcc.noboot.png",sep=""))
#	outliers<-find.outliers(fsts,boot.out=boot.out, 
#		file.name=all.files[i])
	print("finding outliers")
	wcc.outliers<-find.outliers(fsts.wcc,boot.out=wcc.boot.out, 
		file.name=all.files[i])
	if(nrow(wcc.outliers)>0){
		wcc.outliers<-wcc.outliers[wcc.outliers$Ht != 0,]
	wcc.prop<-nrow(wcc.outliers)/(ncol(gpop)-2)
	} else {
		wcc.prop<-0
	}
#	print("p values")
#	wcc.sig<-p.adjust(p.boot(fsts.wcc,wcc.boot.out),"BH")
#	wcc.prop.sig<-length(wcc.sig[wcc.sig<=0.05])/(ncol(gpop)-2)
#	prop<-nrow(outliers)/(ncol(gpop)-2)
	print("return")
	proportions[i,]<-wcc.prop
}

rownames(proportions)<-all.files
proportions$Demes<-as.numeric(gsub("Nm(\\d+.*).d(\\d+).s(\\d+).genepop","\\2",rownames(proportions)))
proportions$Nm<-as.numeric(gsub("Nm(\\d+.*).d(\\d+).s(\\d+).genepop","\\1",rownames(proportions)))
proportions$Samples<-as.numeric(gsub("Nm(\\d+.*).d(\\d+).s(\\d+).genepop","\\3",rownames(proportions)))
props<-data.frame(Samples2=proportions[proportions$Samples == 2,],
                  Samples5=proportions[proportions$Samples == 5,],
                  Samples10=proportions[proportions$Samples == 10,],
                  Samples20=proportions[proportions$Samples == 20,])
props.out<-data.frame(Demes=props$Samples2.Demes,Nm=props$Samples2.Nm,
                        Samples2=props$Samples2.wcc.prop,Samples5=props$Samples5.wcc.prop,
                        Samples10=props$Samples10.wcc.prop,Samples20=props$Samples20.wcc.prop)

props.out<-props.out[order(props.out$Demes,props.out$Nm),]

write.table(props.out,"ProportionOutliers_WCC.23.02.2017.txt",sep="\t",quote=F,row.names=T,
	col.names=T)


##Analyze lositan output
setwd("~/Desktop/numerical_analysis_genepop/")

iam.ci<-list.files(pattern="genepop.ci")
iam.loci<-list.files(pattern="genepop.loci")
stp.ci<-list.files(pattern="step.ci")
stp.loci<-list.files(pattern="step.loci")
s.ci<-data.frame(filename=character(),ParamSet=character(),
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
  props<-find.los.sig(dat,loc.dat)
  s.ci[i,]<-cbind(locus.name,params,nm,d,s,props[1],props[2],props[3])
}

write.csv(s.ci,"StepwiseLositanOutliers.csv")

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
  i.ci[i,]<-cbind(locus.name,params,nm,d,s,outliers)
}

is.ci<-merge(s.ci,i.ci,by="ParamSet")
t.test(as.numeric(is.ci$PropOutliers.x),as.numeric(is.ci$PropOutliers.y),paired=T)

#t = 1.8257, df = 83, p-value = 0.07149


write.csv(i.ci,"InfiniteAllelesModel.csv")

#compare to lositan analysis
s.ci<-read.csv("StepwiseLositanOutliers.csv")
s.ci$filename<-gsub("(Nm\\d+.*.genepop).step.loci","\\1",s.ci$filename)
i.ci<-read.csv("InfiniteAllelesModel.csv")
i.ci$filename<-gsub("(Nm\\d+.*.genepop).loci","\\1",i.ci$filename)
proportions<-read.table("ProportionOutliers_WCC.txt")
proportions$filename<-rownames(proportions)

#step
step.prop<-merge(s.ci,proportions,by="filename")
t.test(x=step.prop$PropOutliers,y=step.prop$wcc.prop,paired=T,alternative="less")

#        Paired t-test
#
#data:  step.prop$PropOutliers and step.prop$wcc.prop
#t = 5.4663, df = 83, p-value = 2.363e-07
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
# 0.06049898        Inf
#sample estimates:
#mean of the differences 
#             0.08696166 

#Correlation between fst and beta
library(Hmisc)
fst.cor<-list()
ht.cor<-list()
h.plots<-list()
f.plots<-list()
  #data.frame(file=character(),Ht=numeric(),Fst=numeric(),Hb=numeric(),beta=numeric(),
          #          stringsAsFactors = FALSE)
nm1.files<-all.files[grep("Nm1\\..*",all.files)]
png("Nm1_correlations_fst.png",height=10,width=7.5,units="in",res=300)
par(mfrow=c(7,4),oma=c(1.5,1.5,1,1),mar=c(1,1,1,1))
for(i in 1:length(nm1.files))
{
  gpop<-my.read.genepop(nm1.files[i])
  fsts.wcc<-calc.actual.fst(gpop,"WCC")
  fsts<-calc.actual.fst(gpop)
  #hp<-plot(fsts$Ht,fsts.wcc$Ht,plot=FALSE)
  plot(fsts$Fst,fsts.wcc$Fst,pch=19,axes=F,xlab="",ylab="")
  axis(1)
  axis(2)
  text(0.35,0.8,nm1.files[i],col="red")
  text(0.35,0.7,paste("r=",round(fst.cor[[i]][[1]][2],3),", p=",fst.cor[[i]][[3]][2]),col="red")
  #h.plots[[i]]<-hp
#  f.plots[[i]]<-fp
#  ht.cor[[i]]<- rcorr(fsts$Ht,fsts.wcc$Ht)
 # fst.cor[[i]]<-rcorr(as.numeric(fsts$Fst[fsts$Fst != NaN]),fsts.wcc$Fst[fsts.wcc$Fst != NaN])
}
mtext(expression(italic(F)[HT]),1,outer=T,line=0.5)
mtext(expression(italic(F)[HB]),2,outer=T,line=0.5)

png("Nm1_correlations_ht.png",height=10,width=7.5,units="in",res=300)
par(mfrow=c(7,4),oma=c(1.5,1.5,1,1),mar=c(1,1,1,1))
for(i in 1:length(nm1.files))
{
  gpop<-my.read.genepop(nm1.files[i])
  fsts.wcc<-calc.actual.fst(gpop,"WCC")
  fsts<-calc.actual.fst(gpop)
  #hp<-plot(fsts$Ht,fsts.wcc$Ht,plot=FALSE)
  plot(fsts$Ht,fsts.wcc$Ht,pch=19,axes=F,xlab="",ylab="")
  axis(1)
  axis(2)
  text(0.35,0.8,nm1.files[i],col="red")
  text(0.35,0.7,paste("r=",round(ht.cor[[i]][[1]][2],3),", p=",ht.cor[[i]][[3]][2]),col="red")
  #h.plots[[i]]<-hp
  #  f.plots[[i]]<-fp
  #  ht.cor[[i]]<- rcorr(fsts$Ht,fsts.wcc$Ht)
  # fst.cor[[i]]<-rcorr(as.numeric(fsts$Fst[fsts$Fst != NaN]),fsts.wcc$Fst[fsts.wcc$Fst != NaN])
}
mtext(expression(italic(H)[T]),1,outer=T,line=0.5)
mtext(expression(italic(H)[B]),2,outer=T,line=0.5)
#rcorr(fst.cor$Ht,fst.cor$Hb)
#los.sig<-read.table("B:/ubuntushare/fst_outliers/results/numerical_analysis_genepop/sig.loci.sim.LOSITAN.txt",skip=1)
#colnames(los.sig)<-c("file","bal","pos","total","balp","posp","totalp")
#these are all at the 95% level
#do a paired t-test.
#los.sig<-los.sig[order(los.sig$file),]
#proportions<-proportions[order(rownames(proportions)),]
#t.test(los.sig$totalp, proportions[,2], paired = TRUE)
#t.test(los.sig$totalp,proportions[,2],paired=T,alternative="greater")
#max.boot<-read.csv("Nm10.d50.s2.genepop95.csv")

## Does bootstrapping reduce the influence of outlier clusters?
par(mfrow=c(1,2),oma=c(2,2,2,2))

gpop<-my.read.genepop("Nm0.1.d5.s20.genepop")
fsts.wcc<-calc.actual.fst(gpop,"wcc")
boot.out<-as.data.frame(t(replicate(10, fst.boot(gpop))))
boot.pvals<-p.boot(fsts.wcc,boot.out=boot.out)
boot.cor.pvals<-p.adjust(boot.pvals,method="BH")
boot.sig<-boot.cor.pvals[boot.cor.pvals <= 0.05]
boot.outliers<-find.outliers(fsts.wcc,boot.out)
non.boot.out<-as.data.frame(t(replicate(1, fst.boot(gpop,bootstrap = FALSE))))
non.boot.pvals<-p.boot(fsts.wcc,boot.out=non.boot.out)
non.boot.cor.pvals<-p.adjust(non.boot.pvals,method="BH")
non.boot.sig<-non.boot.cor.pvals[non.boot.cor.pvals <= 0.05]
non.boot.outliers<-find.outliers(fsts.wcc,non.boot.out)

plotting.cis(fsts.wcc,boot.out,make.file=F,sig.list=names(boot.sig))
plotting.cis(fsts.wcc,non.boot.out,make.file=F,sig.list=names(non.boot.sig))
dim(boot.outliers)
dim(non.boot.outliers)

