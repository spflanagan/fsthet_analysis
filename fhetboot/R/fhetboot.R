#Author: Sarah P. Flanagan
#Last updated 24 June 2016
#Latest update allows haploid genepop format and omits missing data
#Date: 15 February 2016
#Purpose: bootstrap over a genepop file to generate a null distribution
#and then use Johnson distribution to generate confidence intervals.


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
	txt <- txt[-1] #remove the file desription
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
	if(pops[1] == "POP"){
		for(i in 1:length(pops)){ #if there aren't pop names they're given
			pops[i]<-paste("POP ",i,sep="") } #numbers
	}
	popinfo<-as.vector(unlist(apply(cbind(
		as.character(pops),as.numeric(nind)),
		1,function(x){ rep(x[1],x[2])})))	
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

allele.counts<-function(genotypes){
  obs.gen<-summary(as.factor(genotypes))
  obs.gen<-obs.gen[names(obs.gen) != "000000" & names(obs.gen) 
                   != "0000" & names(obs.gen) != "0"]
  if(nchar(names(obs.gen[1])) %% 2 == 0){
    splitnum<-nchar(names(obs.gen[1]))/2
    allele.names<-levels(as.factor(c(substr(names(obs.gen),1,splitnum),
                                     substr(names(obs.gen),splitnum+1,nchar(names(obs.gen[1]))))))
    alleles<-cbind(substr(genotypes,1,splitnum),
                   substr(genotypes,splitnum+1,nchar(names(obs.gen[1]))))
    alleles<-alleles[alleles != "0" & alleles != "00" &
                       alleles != "000" & alleles!= "0000" &  
                       alleles != "00000" & alleles != "000000"]
  } else { #assume it's haploid
    allele.names<-levels(as.factor(names(obs.gen)))
    alleles<-cbind(genotypes,genotypes)
    alleles<-alleles[alleles != "0" & alleles != "00" &
                       alleles != "000" & alleles!= "0000" &  
                       alleles != "00000" & alleles != "000000"]
  }
  obs<-summary(as.factor(alleles))
  if(length(obs) < length(allele.names)){
    num.missing<-length(allele.names[!(allele.names%in% names(obs))])
    AlleleCounts<-c(obs,rep(0,num.missing))
    names(AlleleCounts)<-c(names(obs)[names(obs) %in% allele.names],
                           allele.names[!(allele.names%in% names(obs))])   
  }else{
    AlleleCounts<-obs
  }
  return(AlleleCounts)
}

calc.allele.freq<-function(genotypes){
  counts<-allele.counts(genotypes)
  obs.af<-counts/sum(counts)
  return(obs.af)
}

calc.exp.het<-function(af){ 
	sqaf<-af*af
	ht<-1-sum(sqaf)
	return(ht)
}

calc.fst<-function(df,i){
	df.split<-split(df[,i],df[,1])
	af<-do.call("rbind",lapply(df.split,calc.allele.freq))
	hexp<-apply(af,1,calc.exp.het)
	hs<-mean(hexp)
	ht<-calc.exp.het(calc.allele.freq(df[,i]))
	fst<-(ht-hs)/ht
	return(c(ht,fst))
}

wc.fst<-function(df,i){
  ###ADAPTED FROM LOTTERHOS & WHITLOCK (2014) 
  #Data from: Evaluation of demographic history and neutral parameterization on the performance of Fst outlier tests. 
  #Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.v8d05
  #FSTcalcs.R code
  
  #This is the calculation of beta from Weir and Cockerham (1993)
  df.split<-split(df[,i],df[,1])
  AllCounts<-do.call("rbind",lapply(df.split,allele.counts))
  pops <- rownames(AllCounts)
  numpops <- dim(AllCounts)[1]
  numalleles <- dim(AllCounts)[2]	
  sample.sizes <- rowSums(AllCounts)
  num.inds<-unlist(lapply(df.split,length))
  af<-do.call("rbind",lapply(df.split,calc.allele.freq))
  X <- sum(af^2)
  Y <- sum(colSums(af)^2)
  n<-unlist(lapply(df.split,length))
  M <- num.inds[1]	#uncorrected for sample size
  F0 <- (M*X-numpops)/((M-1)*numpops) 
  F1 <- (Y-X)/(numpops*(numpops-1))
  He <- 1-F1
  fst <- (F0-F1)/(1-F1)
  #FST <- 1 - (1-F0)/(1-F1); an alternative way of writing the previous line
  return(c(He, fst))
}

wc.corr.fst<-function(df,i){
  ###ADAPTED FROM LOTTERHOS & WHITLOCK (2014) 
  #Data from: Evaluation of demographic history and neutral parameterization on the performance of Fst outlier tests. 
  #Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.v8d05
  #FSTcalcs.R code
  
  #This is the code implemented by FDIST2 in my.thetacal
  #Beaumont & Nichols (1996) changed the calculation of beta from
  #Cockerham & Weir (1993) by adding a sample size correction.
  df.split<-split(df[,i],df[,1])
  AllCounts<-do.call("rbind",lapply(df.split,allele.counts))
  af<-do.call("rbind",lapply(df.split,calc.allele.freq))
  sample.sizes <- rowSums(AllCounts) #these sample sizes are the number of alleles sampled.
  pops <- rownames(AllCounts)
  numpops <- dim(AllCounts)[1]
  numalleles <- dim(AllCounts)[2]
  X <- sum(af^2)
  Y <- sum(colSums(af)^2)
  n<-unlist(lapply(df.split,length))
  x0 <- sum((rowSums(AllCounts^2)-sample.sizes)/(sample.sizes*(sample.sizes-1)))
  yy=0
  for (j in 1:(numpops-1)){
    for (k in (j+1):numpops){
      y1=0;
      for (i in 1:numalleles){
        y1 <- y1 + AllCounts[j,i]*AllCounts[k,i]
        #print(y1)
      }
      yy <- yy + y1/(sample.sizes[j]*sample.sizes[k])
    }
  }
  
  q2<- x0/numpops
  q3 <- 2*yy/(numpops*(numpops-1))
  
  het0 <- 1-q2
  het1 <- 1-q3 #this is what is output as heterozygosity by fdist2
  FST <- (q2-q3)/(1-q3)
  #FST <- 1-het0/het1 #alternative way of writing    
  return(c(het1, FST))
}

fst.boot.onecol<-function(df){
	col<-sample(3:ncol(df),1)
	ht.fst<-calc.fst(df,col)
	return(ht.fst)
}

wc.boot.onecol<-function(df){
  col<-sample(3:ncol(df),1)
  ht.fst<-wc.fst(df,col)
  return(ht.fst)
}

wcc.boot.onecol<-function(df){
  col<-sample(3:ncol(df),1)
  ht.fst<-wc.corr.fst(df,col)
  return(ht.fst)
}

fst.options.print<-function()
{
  print("For Nei's Fst formulation: nei, Nei, NEI, N",quote=F)
  print("For Weir and Cockerham (1993)'s Fst formulation: WeirCockerham, WC, weircockerham, wc",quote=F)
  print("For Beaumont's sample-size-corrected version of Weir and Cockerhams (1993)'s Fst formulation:
        WeirCockerhamCorrected, WCC, weircockerhamcorrected, wcc, corrected",quote=F)
}

fst.boot<-function(df, fst.choice="nei"){	
  #Fst options are Nei, WeirCockerham, or WeirCockerhamCorrected
  fst.options<-c("nei", "Nei","NEI","N","WeirCockerham","WC", "weircockerham","wc",
                 "WeirCockerhamCorrected","WCC", "weircockerhamcorrected","wcc", "corrected")
  if(!(fst.choice %in% fst.options)) { stop("Fst choice not an option. Use fst.options.print() to see options.")}
	nloci<-(ncol(df)-2)
	if(fst.choice == "nei" | fst.choice == "Nei" | fst.choice == "NEI" | fst.choice == "N"){
	  boot.out<-as.data.frame(t(replicate(nloci, fst.boot.onecol(df))))}
	if(fst.choice == "WeirCockerham" | fst.choice == "WC" | fst.choice == "weircockerham" | fst.choice == "wc"){
	  boot.out<-as.data.frame(t(replicate(nloci, wc.boot.onecol(df))))}
	if(fst.choice == "WeirCockerhamCorrected" | fst.choice == "WCC" | 
	   fst.choice == "weircockerhamcorrected" | fst.choice == "wcc" | fst.choice == "corrected"){
	  boot.out<-as.data.frame(t(replicate(nloci, wcc.boot.onecol(df))))}
	colnames(boot.out)<-c("Ht","Fst")
	print("Bootstrapping done. Now Calculating CIs")
	#order by het
	boot.out<-as.data.frame(boot.out[order(boot.out$Ht),])
	#create overlapping bins 
	breaks<-hist(boot.out$Ht,breaks=25,plot=F)$breaks
	br.rate<-0.1
	newbreaks<-breaks-(br.rate/2)
#	if((length(breaks) %% 2)==0){ #if it's even
#		bins<-rbind(cbind(breaks[seq(1,length(breaks),2)],
#			breaks[seq(2,length(breaks),2)]),
#			cbind(newbreaks[seq(1,length(newbreaks),2)],
#			newbreaks[seq(2,length(newbreaks),2)]))
#	} else {
#		bins<-rbind(cbind(breaks[seq(1,length(breaks)+1,2)],
#			breaks[seq(2,length(breaks)+1,2)]),
#			cbind(newbreaks[seq(1,length(newbreaks)+1,2)],
#			newbreaks[seq(2,length(newbreaks)+1,2)]))
#	}
#	bins<-bins[order(bins[,1]),]
	bins<-as.data.frame(cbind(newbreaks,breaks))
	mids<-apply(bins,1,mean)
	#bin fsts
	bin.fst<-apply(bins, 1, function(x){ #this returns a list of Fst vectors
		out<-boot.out[boot.out$Ht > x[1] &	boot.out$Ht < x[2],"Fst"] })
	names(bin.fst)<-bins$breaks
	#merge those with too few with the next bin up
	rmvec<-NULL
	for(i in 1:(length(bin.fst)-1)){
		if(length(bin.fst[[i]]) < 20){	
			bin.fst[[(i+1)]]<-c(bin.fst[[i]],bin.fst[[(i+1)]])
			rmvec<-c(rmvec,i)
		}
		if((i+1)==length(bin.fst)){
			if(length(bin.fst[[i+1]]) < 20){	
				bin.fst[[i]]<-c(bin.fst[[i]],bin.fst[[(i+1)]])
				rmvec<-c(rmvec,(i+1))
			}
		}
	}
	#remove bins with too few
	bin.fst<-bin.fst[-rmvec]
	bin.fst<-lapply(bin.fst,sort)	
	#find 95% quantile
	fst95<-lapply(bin.fst, function(x){
		keep.fst<-x[round(length(x)*0.025):round(length(x)*0.975)]
		if(length(keep.fst)>0){
			fst.thresh<-c(min(keep.fst),max(keep.fst))
			names(fst.thresh)<-c("Low95","Upp95")
		} else{
			fst.thresh<-c("","")
			names(fst.thresh)<-c("Low95","Upp95")
		}
		return(fst.thresh)
	})
	fst.95<-do.call(rbind,fst95)
	#find 99% quantile
	fst.99<-lapply(bin.fst, function(x){
		keep.fst<-x[round(length(x)*0.005):round(length(x)*0.995)]
		if(length(keep.fst)>0){
			fst.thresh<-c(min(keep.fst),max(keep.fst))
		} else {
			fst.thresh<-c("","")
		}
		names(fst.thresh)<-c("Low99","Upp99")
		return(fst.thresh)
	})
	fst.99<-do.call(rbind,fst.99)
	return(list(Fsts=boot.out,CI95=fst.95,CI99=fst.99))
}


ci.means<-function(boot.out.list){ #should be boot.out[[2]] or boot.out[[3]]
	if(class(boot.out.list)=="list") {
		boot.ci<-as.data.frame(do.call(rbind,boot.out.list))
	} else {
		boot.ci<-as.data.frame(boot.out.list)
	}
	boot.ci$Ht<-rownames(boot.ci)
	avg.cil<-tapply(boot.ci[,1],boot.ci$Ht,mean)
	avg.ciu<-tapply(boot.ci[,2],boot.ci$Ht,mean)
	return(list(avg.cil,avg.ciu))
}

plotting.cis<-function(df,boot.out=NULL,ci.list=NULL,Ht.name="Ht",Fst.name="Fst",
	ci.col=c("red","gold"), pt.pch=1,file.name=NULL,
	make.file=TRUE) {
#This function takes a dataframe with empirical Fst and Ht measurements
#It must have at least two columns, one named "Ht" and one named "Fst"
#Or you must pass the column names to the function
#You must give it bootstrap output or a list of confidence interval values.
	if(is.null(boot.out) & is.null(ci.list)){
		 stop("Must provide bootstrap output or a list of CI values") 
	} else if(is.null(ci.list)){
		avg.ci95<-ci.means(boot.out[[2]])
		avg.ci99<-ci.means(boot.out[[3]])
	} else {
		avg.ci95<-list(as.numeric(ci.list$low95),as.numeric(ci.list$upp95))
		names(avg.ci95[[1]])<-rownames(ci.list)	
		names(avg.ci95[[2]])<-rownames(ci.list)
		avg.ci99<-list(as.numeric(ci.list$low99),as.numeric(ci.list$upp99))
		names(avg.ci99[[1]])<-rownames(ci.list)	
		names(avg.ci99[[2]])<-rownames(ci.list)
	}
	if(names(avg.ci95[[1]])[1] != "0"){
		avg.ci95[[1]]<-c(0,avg.ci95[[1]])
		avg.ci95[[2]]<-c(0,avg.ci95[[2]])
		avg.ci99[[1]]<-c(0,avg.ci99[[1]])
		avg.ci99[[2]]<-c(0,avg.ci99[[2]])
		names(avg.ci95[[1]])[1]<-0
		names(avg.ci95[[2]])[1]<-0
		names(avg.ci99[[1]])[1]<-0
		names(avg.ci99[[2]])[1]<-0
	}
	if(make.file==TRUE){
		if(!is.null(file.name)) {
			png(file.name,height=8,width=9,units="in",res=300) }
		else {
			png("OutlierLoci.png",height=8,width=9,units="in",res=300) }
	}
	x.max<-round(as.numeric(max(df[,Ht.name]))+0.1,1)
	plot(df[,Ht.name],df[,Fst.name],xlab="",ylab="",las=1,pch=pt.pch,axes=F,
		xlim=c(0,x.max))
	axis(1,pos=0,at=seq(0,x.max,0.1))
	axis(2,pos=0,las=1)
	mtext(expression("F"["ST"]),2,line=2.5)
	mtext(expression("H"["T"]),1,line=2.5)
	points(names(avg.ci95[[1]]),avg.ci95[[1]],type="l",col=ci.col[1])
	points(names(avg.ci95[[2]]),avg.ci95[[2]],type="l",col=ci.col[1])
	points(names(avg.ci99[[1]]),avg.ci99[[1]],type="l",col=ci.col[2])
	points(names(avg.ci99[[2]]),avg.ci99[[2]],type="l",col=ci.col[2])
	legend(x=0.01,y=max(df[,Fst.name]),c("95% CI","99% CI"),
		col=ci.col,lty=1,bty='n')
	if(make.file==TRUE) dev.off()
	
}

find.outliers<-function(df,ci.df=NULL,boot.out=NULL,file.name=NULL){
#Need to either give this function bootstrap output or a list of CIs
	if(is.null(boot.out) & is.null(ci.df)){
		stop("Must provide bootstrap output or a list of CI values") 
	} else if(is.null(ci.df)){
		avg.ci95<-ci.means(boot.out[[2]])
		avg.ci99<-ci.means(boot.out[[3]])
		ci.df<-as.data.frame(t(do.call(rbind,c(avg.ci95,avg.ci99))))
		colnames(ci.df)<-c("low95","upp95","low99","upp99")
		ci.df$Ht<-as.numeric(rownames(ci.df))
	}
	diff<-0
	for(i in 2:nrow(ci.df)){
		diff<-c(diff,ci.df$Ht[i]-ci.df$Ht[(i-1)])
	}
	bin<-cbind(ci.df$Ht-(diff/2),ci.df$Ht+(diff/2))
	actual.bin<-apply(bin, 1, function(x){ #this returns a list of Fst vectors
		out<-df[df$Ht > x[1] &	df$Ht < x[2],] })
	out95<-NULL
	out99<-NULL
	for(i in 1:nrow(ci.df)){
		out95<-rbind(out95,actual.bin[[i]][
			actual.bin[[i]]$Fst< ci.df[i,"low95"] | 
			actual.bin[[i]]$Fst> ci.df[i,"upp95"],])
		out99<-rbind(out99,actual.bin[[i]][
			actual.bin[[i]]$Fst< ci.df[i,"low99"] | 
			actual.bin[[i]]$Fst> ci.df[i,"upp99"],])
	}
	if(!is.null(file.name)){
		write.csv(out95,paste(file.name,"95.csv",sep=""))
		write.csv(out99,paste(file.name,"99.csv",sep=""))
	}
	return(list(out95,out99))
}

calc.actual.fst<-function(df, fst.choice="N"){
	fst.options<-c("nei", "Nei","NEI","N","WeirCockerham","WC", "weircockerham","wc",
	               "WeirCockerhamCorrected","WCC", "weircockerhamcorrected","wcc", "corrected")
	if(!(fst.choice %in% fst.options)) { stop("Fst choice not an option. Use fst.options.print() to see options.")}
	fsts<-data.frame(Locus=character(),Ht=numeric(),Fst=numeric(),
	                 stringsAsFactors=F)
	
	if(fst.choice == "nei" | fst.choice == "Nei" | fst.choice == "NEI" | fst.choice == "N"){
	  for(i in 3:ncol(df)){
	    fsts[i-2,]<-c(colnames(df)[i],calc.fst(df,i))
	  }
	}
	if(fst.choice == "WeirCockerham" | fst.choice == "WC" | fst.choice == "weircockerham" | fst.choice == "wc"){
	  for(i in 3:ncol(df)){
	  fsts[i-2,]<-c(colnames(df)[i],wc.fst(df,i))
	  }
	}
	if(fst.choice == "WeirCockerhamCorrected" | fst.choice == "WCC" | 
	   fst.choice == "weircockerhamcorrected" | fst.choice == "wcc" | fst.choice == "corrected"){
	  for(i in 3:ncol(df)){
	    fsts[i-2,]<-c(colnames(df)[i],wc.corr.fst(df,i))
	  }
	}
	
	return(fsts)
}