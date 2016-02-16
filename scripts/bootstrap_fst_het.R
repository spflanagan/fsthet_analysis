#Author: Sarah P. Flanagan
#Date: 15 February 2016
#Purpose: bootstrap over a genepop file to generate a null distribution
#and then use Johnson distribution to generate confidence intervals.

#install.packages("SuppDists")
library(SuppDists)
library(boot)

remove.spaces<-function (charvec) 
{#adapted from adegenet
    charvec <- gsub("^([[:blank:]]*)([[:space:]]*)", "", charvec)
    charvec <- gsub("([[:blank:]]*)([[:space:]]*)$", "", charvec)
    return(charvec)
}

my.read.genepop<-function (file, ncode = 2L, quiet = FALSE) 
{#modified from the adegenet function read.genepop
	if (!quiet) 
		cat("\nParsing Genepop file...\n\n")
	prevcall <- match.call()
	txt <- scan(file, sep = "\n", what = "character", quiet = TRUE)
	if (!quiet) 
		cat("\nFile description: ", txt[1], "\n")
	txt <- txt[-1]
	txt <- gsub("\t", " ", txt)
	#extract the locus info
	locinfo.idx <- 1:(min(grep("POP", toupper(txt))) - 1)
	locinfo <- txt[locinfo.idx]
	locinfo <- paste(locinfo, collapse = ",")
	loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
	loc.names <- remove.spaces(loc.names)
	nloc <- length(loc.names)
	#Now discard the locus info.
	txt <- txt[-locinfo.idx]
	pop.idx <- grep("POP", toupper(txt))
	npop <- length(pop.idx)
	nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
	splited <- nocomma[which(!nocomma %in% pop.idx)]
	if (length(splited) > 0) {
		for (i in sort(splited, decreasing = TRUE)) {
			txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
		}
		txt <- txt[-splited]
	}
	#extract the population information
	pop.idx <- grep("POP", toupper(txt))
	txt[length(txt) + 1] <- "POP"
	pops<-txt[pop.idx]
	nind<-diff(c(pop.idx,length(txt)))-1
	popinfo<-unlist(apply(cbind(pops,nind),1,function(x){ rep(x[1],x[2])}))	
	#Now keep just the genotype info
	txt <- txt[-c(pop.idx, length(txt))]
	temp <- sapply(1:length(txt), function(i) strsplit(txt[i], ","))
	ind.names <- sapply(temp, function(e) e[1])
	ind.names <- remove.spaces(ind.names)
	vec.genot <- sapply(temp, function(e) e[2])
	vec.genot <- remove.spaces(vec.genot)	
	X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), 
		ncol = nloc, byrow = TRUE)
	rownames(X) <- 1:nrow(X)
	colnames(X) <- loc.names
	#combine pop info, ind names, and genotypes
	res<-data.frame(popinfo,ind.names,X)
      if (!quiet) 
		cat("\n...done.\n\n")
	return(res)
}

calc.allele.freq<-function(genotypes){
	obs.gen<-summary(as.factor(genotypes))
	alleles<-cbind(substr(genotypes,1,2),substr(genotypes,3,4))
	obs.af<-summary(as.factor(alleles))/sum(summary(as.factor(alleles)))
	return(obs.af)
}

calc.exp.het<-function(af){ 
	sqaf<-af*af
	ht<-1-sum(sqaf)
	return(ht)
}

calc.fst<-function(df,i){
	df.split<-split(df[,i],df[,1])
	af<-lapply(df.split,calc.allele.freq)
	hexp<-unlist(lapply(af,calc.exp.het))
	n<-unlist(lapply(df.split,length))
	hs<-sum(hexp*n)/sum(n)
	ht<-calc.exp.het(calc.allele.freq(df[,i]))
	fst<-(ht-hs)/ht
	return(c(ht,fst))
}

bin.cis<-function(ht.fst){
	ht.fst<-as.data.frame(ht.fst[order(ht.fst[,1]),])
	
}

fst.onerow<-function(df){
	row<-sample(3:ncol(df),1)
	ht.fst<-calc.fst(df,row)
	return(ht.fst)
}

fst.boot<-function(df){	
	nloci<-(ncol(df)-2)
	boot.out<-as.data.frame(t(replicate(nloci, fst.boot(gpop))))
	colnames(boot.out)<-c("Ht","Fst")
	#order by het
	boot.out<-as.data.frame(boot.out[order(boot.out$Ht),])
	#create overlapping bins 
	breaks<-hist(boot.out$Ht,breaks=20)$breaks
	br.rate<-breaks[2]-breaks[1]
	newbreaks<-breaks-(br.rate/2)
	bins<-rbind(cbind(breaks[seq(1,length(breaks),2)],
		breaks[seq(2,length(breaks),2)]),
		cbind(newbreaks[seq(1,length(newbreaks),2)],
		newbreaks[seq(2,length(newbreaks),2)]))
	bins<-bins[order(bins[,1]),]
	mids<-apply(bins,1,mean)
	#bin fsts
	bin.fst<-apply(bins, 1, function(x){ #this returns a list of Fst vectors
		out<-boot.out[boot.out$Ht > x[1] &	boot.out$Ht < x[2],"Fst"] })
	names(bin.fst)<-mids
	for(i in 1:(length(bin.fst)-1)){
		if(length(bin.fst[[i]]) < 20){	
			bin.fst[[(i+1)]]<-c(bin.fst[[i]],bin.fst[[(i+1)]])
			bin.fst<-bin.fst[-i]
		}
		if((i+1)==length(bin.fst)){
			if(length(bin.fst[[i+1]]) < 20){	
				bin.fst[[i]]<-c(bin.fst[[i]],bin.fst[[(i+1)]])
				bin.fst<-bin.fst[-(i+1)]
			}
		}
	}
	bin.fst<-lapply(bin.fst,sort)	
	#find 95% quantile
	fst95<-lapply(bin.fst, function(x){
		keep.fst<-x[round(length(x)*0.025):round(length(x)*0.975)]
		fst.thresh<-c(min(keep.fst),max(keep.fst))
		names(fst.thresh)<-c("Low95","Upp95")
		return(fst.thresh)
	})
	fst.95<-do.call(rbind,fst95)
	#find 99% quantile
	fst.99<-lapply(bin.fst, function(x){
		keep.fst<-x[round(length(x)*0.005):round(length(x)*0.995)]
		fst.thresh<-c(min(keep.fst),max(keep.fst))
		names(fst.thresh)<-c("Low99","Upp99")
		return(fst.thresh)
	})
	fst.99<-do.call(rbind,fst99)
	return(list(Fsts=boot.out,CI95=fst.95,CI99=fst.99))
}
setwd("E:/ubuntushare/fst-het/data_from_literature")
gpop<-my.read.genepop("Hess_2013_data_Genepop.gen") #this worked.

fsts<-data.frame(Locus=character(),Ht=numeric(),Fst=numeric(),stringsAsFactors=F)
for(i in 3:ncol(gpop)){
	fsts[i-2,]<-c(colnames(gpop)[i],calc.fst(gpop,i))
}

boot.out<-as.data.frame(t(replicate(10000, fst.boot(gpop))))
colnames(boot.out)<-c("Ht","Fst")
boot.out<-boot.out[order(boot.out$Ht),]

boot.sub<-boot.out[1:400,]
johns.dist<-FitJohnsonDistribution(mean(boot.sub$Fst),var(boot.sub$Fst),
	skewness(boot.sub$Fst),kurtosis(boot.sub$Fst))

