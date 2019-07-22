#Author: Sarah P. Flanagan
#27 Feb 2017: change to focus on quantiles instead of bootstrapping
#3 Dec 2016: add fhetboot function to conduct the basic analysis
#2 Dec 2016: allow user to specify confidence intervals
#2 Dec 2016: Calculate p-values from bootstrapped distribution
#27 June 2016: allows haploid genepop format and omits missing data
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
  txt <- txt[-1]
  txt <- gsub("\t", " ", txt)
  locinfo.idx <- 1:(min(grep("POP", toupper(txt))) - 1)
  locinfo <- txt[locinfo.idx]
  locinfo <- paste(locinfo, collapse = ",")
  loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
  loc.names <- remove.spaces(loc.names)
  nloc <- length(loc.names)
  txt <- txt[-locinfo.idx]
  pop.idx <- grep("POP$", toupper(txt))
  npop <- length(pop.idx)
  nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
  splited <- nocomma[which(!nocomma %in% pop.idx)]
  if (length(splited) > 0) {
    for (i in sort(splited, decreasing = TRUE)) {
      txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
    }
    txt <- txt[-splited]
  }
  txt[length(txt) + 1] <- "POP"
  pops <- txt[pop.idx]
  nind <- diff(c(pop.idx, length(txt))) - 1
  if (pops[1] %in% c("POP", "pop", "Pop")) {
    for (i in 1:length(pops)) {
      pops[i] <- paste("POP ", i, sep = "")
    }
  }
  popinfo <- as.vector(unlist(apply(cbind(as.character(pops), 
                                          as.numeric(nind)), 1, function(x) {
                                            rep(x[1], x[2])
                                          })))
  txt <- txt[-c(pop.idx, length(txt))]
  temp <- sapply(1:length(txt), function(i) strsplit(txt[i], 
                                                     ","))
  ind.names <- sapply(temp, function(e) e[1])
  ind.names <- remove.spaces(ind.names)
  vec.genot <- sapply(temp, function(e) e[2])
  vec.genot <- remove.spaces(vec.genot)
  X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), 
              ncol = nloc, byrow = TRUE)
  rownames(X) <- 1:nrow(X)
  colnames(X) <- loc.names
  res <- data.frame(popinfo, ind.names, X)
  if (!quiet) 
    cat("\n...done.\n\n")
  return(res)
}

allele.counts<-function(genotypes){
  obs.gen<-table(as.factor(genotypes))
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
  AlleleCounts<-AlleleCounts[order(names(AlleleCounts))]
  return(AlleleCounts)
}

calc.allele.freq<-function(genotypes){
  counts<-allele.counts(genotypes)
  obs.af<-counts/sum(counts)
  return(obs.af)
}

calc.exp.het<-function(af){ 
	#sqaf<-af*af
	#ht<-1-sum(sqaf)
  ht<-2*af[1]*(1-af[1])
	return(ht)
}

calc.fst<-function(df,i){
	df.split<-split(df[,i],df[,1])
	af<-do.call("rbind",lapply(df.split,calc.allele.freq))
	hexp<-apply(af,1,calc.exp.het)
	if(length(hexp[is.na(hexp)]) > 0){
	  print(paste("WARNING:", length(hexp[is.na(hexp)]),"populations are missing locus",colnames(df)[i]))
	}
	hs<-mean(hexp,na.rm = TRUE)
	pbar<-mean(af[,1],na.rm=TRUE)
	ht<-2*pbar*(1-pbar)
	fst<-1-(hs/ht)
	return(c(ht,fst))
}

var.fst<-function(df,i){
  df.split<-split(df[,i],df[,1])
  N<-length(df.split)
  af<-do.call("rbind",lapply(df.split,calc.allele.freq))
  if(length(af[is.na(af)]) > 0){
    print(paste("WARNING:", nrow(af[is.na(af[,1]),]),"populations are missing locus",colnames(df)[i]))
    N<-nrow(af[!is.na(af[,1]),])
  }
  pbar<-mean(af[,1],na.rm=TRUE)
  vp<-(af[!is.na(af[,1]),1]-pbar)^2
  varp<-(1/N)*sum(vp)
  fst<-(varp/(pbar*(1-pbar)))-(1/(2*N))
  ht<-2*pbar*(1-pbar)
  return(c(ht,fst))
}

calc.theta<-function(df,i){
  df.split<-split(df[,i],df[,1])
  N<-length(df.split)
  af<-do.call("rbind",lapply(df.split,calc.allele.freq))
  if(length(af[is.na(af)]) > 0){
    print(paste("WARNING:", nrow(af[is.na(af[,1]),]),"populations are missing locus",colnames(df)[i]))
    N<-nrow(af[!is.na(af[,1]),])
  }
  pbar<-mean(af[!is.na(af[,1]),1])
  vp<-(af[!is.na(af[,1]),1]-pbar)^2
  ns<-unlist(lapply(df.split,length))
  nbar<-mean(unlist(lapply(df.split,length)))
  s2a<-(1/((N-1)*nbar))*sum(vp)
  nc<-(1/(N-1))*(sum(ns)-(sum(ns^2)*(1/sum(ns))))
  T1<-s2a-((1/(nbar-1))*((pbar*(1-pbar))-((N-1)/N)*s2a))
  T2<-((nc-1)/(nbar-1))*(pbar*(1-pbar))+(s2a/N)*(1+(((N-1)*(nbar-nc))/(nbar-1)))
  theta<-T1/T2
  return(c(T2,theta))
}

calc.betahat<-function(df,i){
  df.split<-split(df[,i],df[,1])
  r<-length(df.split)
  af<-do.call("rbind",lapply(df.split,calc.allele.freq))
  if(length(af[is.na(af)]) > 0){
    print(paste("WARNING:", nrow(af[is.na(af[,1]),]),"populations are missing locus",colnames(df)[i]))
    N<-nrow(af[!is.na(af[,1]),])
  }
  M<-mean(unlist(lapply(df.split[which(!is.na(af[,1]))],length)))
  
  Y<-sum(af[!is.na(af[,1]),1])*sum(af[!is.na(af[,1]),1])+sum(1-af[!is.na(af[,1]),1])*sum(1-af[!is.na(af[,1]),1])
  X<-sum(af[!is.na(af[,1]),1]^2)+sum((1-af[!is.na(af[,1]),1])^2)
  F0<-((2*M*X)-r)/(((2*M)-1)*r)
  F1<-((Y-X)/(r*(r-1)))
  HB<-1-F1
  betahat<-(F0-F1)/HB
  return(c(HB,betahat))
}

calc.weighted.fst<-function(df,i){
  df.split<-split(df[,i],df[,1])
  af<-do.call("rbind",lapply(df.split,calc.allele.freq))
  hexp<-apply(af,1,calc.exp.het)
  if(length(hexp[is.na(hexp)]) > 0){
    print(paste("WARNING:", length(hexp[is.na(hexp)]),"populations are missing locus",colnames(df)[i]))
  }
  ns<-unlist(lapply(df.split[which(!is.na(hexp))],length))
  hs<-sum(hexp[!is.na(hexp)]*ns)/sum(ns)
  ht<-calc.exp.het(calc.allele.freq(df[,i]))
  fst<-(ht-hs)/ht
  return(c(ht,fst))
}



fst.boot.onecol<-function(df,fst.choice){
  fst.options<-c("FST","Fst","fst","var","VAR","Var","theta","Theta","THETA",
                 "Betahat","BETAHAT","betahat")
  if(!(fst.choice %in% fst.options)) { stop("Fst choice not an option. Use fst.options.print() to see options.")}
  col<-sample(3:ncol(df),1)
  if(fst.choice %in% c("Fst","FST","fst")){
    ht.fst<-calc.fst(df,col)
  }
  if(fst.choice %in% c("VAR","var","Var")){
    ht.fst<-var.fst(df,col)
  }
  if(fst.choice %in% c("Theta","theta","THETA")){
    ht.fst<-calc.theta(df,col)
  }
  if(fst.choice %in% c("betahat","BETAHAT","Betahat")){
    ht.fst<-calc.betahat(df,col)
  }
  return(ht.fst)
}

fst.options.print<-function(){
  print("For Wright's Fst: fst, FST, Fst",quote=F)
  print("For a variance-based Fst (beta): var, VAR, Var",quote=F)
  print("For Cockerham and Weir's Theta: theta, Theta, THETA", quote=F)
  print("For Beta-hat (LOSITAN): Betahat, betahat, BETAHAT",quote=F)
}

make.bins<-function(fsts,num.breaks=25, Ht.name="Ht", Fst.name="Fst",min.per.bin=20)
{
  breaks<-hist(fsts[,Ht.name],breaks=(num.breaks/2),plot=F)$breaks
  low.breaks<-breaks[1:(length(breaks)-1)]
  upp.breaks<-breaks[2:length(breaks)]
  br.rate<-breaks[2]-breaks[1]
  newbreaks<-breaks-(br.rate/2)
  low.breaks<-c(low.breaks,newbreaks[1:(length(breaks)-1)])
  upp.breaks<-c(upp.breaks,newbreaks[2:length(breaks)])
  bins<-data.frame(low.breaks=sort(low.breaks),upp.breaks=sort(upp.breaks))
  bins<-bins[bins[,1]>=0 & bins[,2]>=0,]
  mids<-apply(bins,1,mean)
  #bin fsts
  bin.fst<-apply(bins, 1, function(x){ #this returns a list of Fst vectors
    out<-fsts[fsts[,Ht.name] > x[1] &	fsts[,Ht.name] < x[2],Fst.name] })
  names(bin.fst)<-bins$upp.breaks
  #merge those with too few with the next bin up
  rmvec<-NULL
  for(i in 1:(length(bin.fst)-1)){
    if(length(bin.fst[[i]]) < min.per.bin){	
      bin.fst[[(i+1)]]<-c(bin.fst[[i]],bin.fst[[(i+1)]])
      rmvec<-c(rmvec,i)
    }
    if((i+1)==length(bin.fst)){
      if(length(bin.fst[[i+1]]) < min.per.bin){	
        bin.fst[[i]]<-c(bin.fst[[i]],bin.fst[[(i+1)]])
        rmvec<-c(rmvec,(i+1))
      }
    }
  }
  #remove bins with too few
  if(!is.null(rmvec)){	bin.fst<-bin.fst[-rmvec] }
  bin.fst<-lapply(bin.fst,sort)
  for(i in 1:length(rmvec)){
    bins[(rmvec[i]+1),1]<-bins[rmvec[i],1]
  }
  if(!is.null(rmvec)){ bins<-bins[-rmvec,] }
  #final.bins<-bins[bins$breaks %in% as.numeric(names(bin.fst)),]#no good
  return(list(bins=bins,bin.fst=bin.fst))
}

find.quantiles<-function(bins, bin.fst, ci=0.05)
{
  fst.CI<-list()
  for(j in 1:length(ci)){
    ci.min<-(ci[j]/2)
    ci.max<-1-(ci[j]/2)
    fstCI<-lapply(bin.fst, function(x){
      keep.fst<-x[round(length(x)*ci.min):round(length(x)*ci.max)]
      if(length(keep.fst)>0){
        fst.thresh<-c(min(keep.fst),max(keep.fst))
        names(fst.thresh)<-c("Low","Upp")
      } else{
        fst.thresh<-c("","")
        names(fst.thresh)<-c("Low","Upp")
      }
      return(fst.thresh)
    })		
    ci.name<-paste("CI",(1-ci[j]),sep="")
    cis<-data.frame(do.call(rbind,fstCI))
    cis$UppHet<-as.numeric(rownames(cis))
    cis<-apply(cis,c(1,2),round,3)
    bins<-apply(bins,c(1,2),round,3)
    cis<-merge(cis,bins,by.x="UppHet",by.y="upp.breaks",keep=T)
    fst.CI[[j]]<-data.frame(Low=cis$Low,Upp=cis$Upp,LowHet=cis$low.breaks,UppHet=cis$UppHet)
    names(fst.CI)[j]<-ci.name
  }
  return(fst.CI)
}

fst.boot<-function(df, fst.choice="fst", ci=0.05,num.breaks=25,bootstrap=TRUE,min.per.bin=20){	
		#updated 2 Dec 2016
  #Fst options are Wright, WeirCockerham, or WeirCockerhamCorrected
  fst.options<-c("FST","Fst","fst","var","VAR","Var","theta","Theta","THETA",
                 "Betahat","BETAHAT","betahat")
  if(!(fst.choice %in% fst.options)) { stop("Fst choice not an option. Use fst.options.print() to see options.")}
	nloci<-(ncol(df)-2)
	if(bootstrap == TRUE)
	{
  	boot.out<-as.data.frame(t(replicate(nloci, fst.boot.onecol(df,fst.choice))))
  	colnames(boot.out)<-c("Ht","Fst")
  	boot.out$Fst[boot.out$Fst=="NaN"]<-0
  	print("Bootstrapping done. Now Calculating CIs")
	} else{
	  boot.out<-calc.actual.fst(df,fst.choice)
	  rownames(boot.out)<-boot.out$Locus
	  boot.out<-data.frame(cbind(as.numeric(boot.out$Ht),as.numeric(boot.out$Fst)))
	  colnames(boot.out)<-c("Ht","Fst")
	  boot.out$Fst[boot.out$Fst=="NaN"]<-0
	  print("Fsts calculated. Now Calculating CIs")
	}
	#order by het
	boot.out<-as.data.frame(boot.out[order(boot.out$Ht),])
	#create overlapping bins 
  bins<-make.bins(boot.out,num.breaks,min.per.bin=min.per.bin)
	#find quantile
  fst.CI<-find.quantiles(bins$bins,bins$bin.fst,ci)

	return(list(Fsts=boot.out,Bins=bins$bins,fst.CI))
}

fst.boot.means<-function(boot.out){ #boot.out[[1]]
		#updated 2 Dec 2016
	all.boot<-data.frame(
		Ht=unlist(lapply(boot.out$Fsts, 
			function(x) { out<-x$Ht; return(out) })),
		Fst=unlist(lapply(boot.out$Fsts, 
			function(x) { out<-x$Fst; return(out) })))
	breaks<-boot.out$Bins[[1]]
	bmu<-t(apply(breaks,1,function(x){
		x.ht<-all.boot[all.boot$Ht >= x[1] & 
			all.boot$Ht <= x[2],"Ht"]
		x.fst<-all.boot[all.boot$Ht >= x[1] & 
			all.boot$Ht <= x[2],"Fst"]
		x.fh<-c(mean(x.ht),mean(x.fst),length(x.ht))
		return(x.fh)		
	}))
	bmu<-data.frame(bmu[,1],bmu[,2],bmu[,3],breaks[,1],breaks[,2])
	colnames(bmu)<-c("Ht","Fst","Num","LowBin","UppBin")
	return(bmu)
}

p.boot<-function(actual.fsts, boot.out=NULL,boot.means=NULL){
		#updated 2 Dec 2016
	if(is.null(boot.out) & is.null(boot.means)) { 
		stop("You must provide either bootstrapping output or bootstrap means") }
  if(class(boot.out[[2]])=="data.frame"){
    stop("Can only calculate p-values if bootstrapping was run more than once.") }
	if(is.null(boot.means)){ #then we need to calculate means
		boot.means<-fst.boot.means(boot.out)				
	}#boot.means
	boot.means$real.means<-apply(boot.means,1,function(x){
		rmu<-mean(as.numeric(actual.fsts[
			as.numeric(actual.fsts$Ht) >= as.numeric(x[4]) & 
			as.numeric(actual.fsts$Ht) <= as.numeric(x[5]),"Fst"]),na.rm=T)
		return(rmu)
	})
	boot.means$unitsaway<-abs(boot.means$real.means - boot.means$Fst)
	boot.means$low<-boot.means$Fst-boot.means$unitsaway
	boot.means$upp<-boot.means$Fst+boot.means$unitsaway
	pvals<-unlist(apply(actual.fsts,1, function(x){
		bin<-boot.means[
			as.numeric(boot.means$LowBin) <= as.numeric(x["Ht"])
			& as.numeric(boot.means$UppBin) >= as.numeric(x["Ht"]),]
		fsts.in.bin<-apply(bin,1,function(y){
			actual.fsts[actual.fsts$Ht >= as.numeric(y["LowBin"] )
				& actual.fsts$Ht <= as.numeric(y["UppBin"]),]
		})
		
		unitsaway<-apply(bin,1,function(y){ 
			abs(as.numeric(x["Fst"]) - as.numeric(y["Fst"])) })
		low<-apply(bin,1,function(y){ 
			as.numeric(y["Fst"]) - as.numeric(unitsaway)  })
		upp<-apply(bin,1,function(y){ 
			as.numeric(y["Fst"]) + as.numeric(unitsaway)  })
		n<-lapply(fsts.in.bin, nrow)
		p<-NULL
		if(length(fsts.in.bin)>0){
  		for(i in 1:length(fsts.in.bin)){
  			p[i]<-(nrow(fsts.in.bin[[i]][as.numeric(fsts.in.bin[[i]]$Fst) < 
  				as.numeric(low[i]),]) + 
  				nrow(fsts.in.bin[[i]][as.numeric(fsts.in.bin[[i]]$Fst) > 
  				as.numeric(upp[i]),]))/as.numeric(n[i])
  		}
  		p<-max(p)	#Some of these loci are in multiple bins.
		}else{
		  p<-NA
		  print("No bins found. Were actual.fsts and boot.out/boot.means calculated with the same Fst method?")
		}
		return(p)
	
	}))
	names(pvals)<-actual.fsts$Locus
	return(pvals)
}

ci.means<-function(boot.out.list){ #boot.out[[3]]
	#updated 2 Dec 2016
	if(class(boot.out.list)=="list" & length(boot.out.list) > 1) {
		boot.ci<-as.data.frame(do.call(rbind,
			lapply(boot.out.list,function(x){
				y<-as.data.frame(x[[1]])
				return(y)
			})))
	} else {
		boot.ci<-as.data.frame(boot.out.list[[1]])
	}
  colnames(boot.ci)<-c("Low","Upp","LowHet","UppHet")
	avg.cil<-tapply(boot.ci[,1],boot.ci$UppHet,mean)
	avg.ciu<-tapply(boot.ci[,2],boot.ci$UppHet,mean)
	low.het<-tapply(boot.ci$LowHet,boot.ci$UppHet,mean)
	return(data.frame(Low=avg.cil,Upp=avg.ciu, LowHet=low.het,
		UppHet=as.numeric(levels(as.factor(boot.ci$UppHet)))))
}

plotting.cis<-function(df,boot.out=NULL,ci.df=NULL,sig.list=NULL,Ht.name="Ht",Fst.name="Fst",
	ci.col="red", pt.pch=1,file.name=NULL,sig.col=ci.col,make.file=TRUE) {
#This function takes a dataframe with empirical Fst and Ht measurements
#It must have at least two columns, one named "Ht" and one named "Fst"
#Or you must pass the column names to the function
#You must give it bootstrap output or a list of confidence interval values.
#updated 5 Dec 2016
	if(is.null(boot.out) & is.null(ci.df)){
		 stop("Must provide bootstrap output or a list of CI values") 
	} else if(is.null(ci.df)){
		avg.ci<-ci.means(boot.out[[3]])
	} else {
		avg.ci<-ci.df
	}
  #if(smooth.ci==TRUE){
   # avg.ci<-smooth.cis(avg.ci,smoothing.rate)
  #}
	#if(names(avg.ci[[1]])[1] != "0"){ #need to extend the CIs to 0
	#	avg.ci[[1]]<-c(0,avg.ci[[1]])
	#	avg.ci[[2]]<-c(0,avg.ci[[2]])
	#	names(avg.ci[[1]])[1]<-0
	#	names(avg.ci[[2]])[1]<-0
	#}
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
	mtext(expression(italic(F)["ST"]),2,line=2.5)
	mtext(expression(italic(H)["T"]),1,line=2.5)
	if(!is.null(sig.list)){
	  points(df[df[,1] %in% sig.list,Ht.name],df[df[,1] %in% sig.list,Fst.name],col=sig.col,pch=pt.pch)
	}
	points(avg.ci$LowHet,avg.ci$Low,type="l",col=ci.col)
	points(avg.ci$UppHet,avg.ci$Upp,type="l",col=ci.col)
	if(make.file==TRUE) dev.off()
	
}

find.outliers<-function(df,boot.out,ci.df=NULL,file.name=NULL){
#Updated 12 Dec 2016
#Need to give this function bootstrap output (or a list of CIs)
	if(is.null(boot.out) & is.null(ci.df)){
		stop("Must provide bootstrap output or a list of CI values") 
	} else if(is.null(ci.df)){ #need to calculate the mean ci values
		ci.df<-ci.means(boot.out[[3]])
	}
	if(class(boot.out[[2]])=="list"){
	  bin<-boot.out[[2]][[1]] }
	if(class(boot.out[[2]])=="data.frame"){
	  bin<-boot.out[[2]] }
	#match Fst and Ht
	#I want to find the outliers. 
	#For each df, find its closest cis
	actual.bin<-apply(ci.df, 1, function(x){ #this returns a list of Fst vectors
		out<-df[df$Ht >= x["LowHet"] &	df$Ht < x["UppHet"],] })
	out<-NULL

	for(i in 1:nrow(ci.df)){
		out<-rbind(out,actual.bin[[i]][
			as.numeric(actual.bin[[i]]$Fst) < as.numeric(ci.df[i,"Low"]) | 
			as.numeric(actual.bin[[i]]$Fst)> as.numeric(ci.df[i,"Upp"]),])
		
	}
	out<-out[!duplicated(out$Locus),]

#	outliers<-NULL
#	for(i in 1: nrow(df)){
#		x<-df[i,]
#		#get the closest ones to the ht
#		this.ci<-ci.df[which.min(abs(as.numeric(ci.df$Ht)-as.numeric(x["Ht"]))),]
#		if(x["Fst"] > this.ci$upp | x["Fst"] < this.ci$low){		
#			outliers<-rbind(outliers,x) }
#	}
	if(!is.null(file.name)){
		write.csv(out,paste(file.name,".csv",sep=""))
	}
	return(out)
}

calc.actual.fst<-function(df, fst.choice="fst"){
  fst.options<-c("FST","Fst","fst","var","VAR","Var","theta","Theta","THETA",
                 "Betahat","BETAHAT","betahat")
  if(!(fst.choice %in% fst.options)) { stop("Fst choice not an option. Use fst.options.print() to see options.")}
  fsts<-data.frame(Locus=character(),Ht=numeric(),Fst=numeric(),
                   stringsAsFactors=F)
  if(fst.choice %in% c("Fst","FST","fst")){
    for(i in 3:ncol(df)){
      fsts[i-2,]<-c(colnames(df)[i],calc.fst(df,i))
    }
  }
  if(fst.choice %in% c("VAR","var","Var")){
    for(i in 3:ncol(df)){
      fsts[i-2,]<-c(colnames(df)[i],var.fst(df,i))
    }
  }
  if(fst.choice %in% c("Theta","theta","THETA")){
    for(i in 3:ncol(df)){
      fsts[i-2,]<-c(colnames(df)[i],calc.theta(df,i))
    }
  }
  if(fst.choice %in% c("betahat","BETAHAT","Betahat")){
    for(i in 3:ncol(df)){
      fsts[i-2,]<-c(colnames(df)[i],calc.betahat(df,i))
    }
  }	
	
	fsts<-data.frame(Locus=as.character(fsts$Locus),Ht=as.numeric(fsts$Ht),Fst=as.numeric(fsts$Fst),
	                 stringsAsFactors=F)
	return(fsts)
}

fhetboot<-function(gpop, fst.choice="fst",alpha=0.05,nreps=10){
  fsts<-calc.actual.fst(gpop,fst.choice)
  boot.out<-as.data.frame(t(replicate(nreps, fst.boot(gpop,fst.choice))))
  boot.pvals<-p.boot(fsts,boot.out=boot.out)
  boot.cor.pvals<-p.adjust(boot.pvals,method="BH")
  boot.sig<-boot.cor.pvals[boot.cor.pvals <= alpha]
  plotting.cis(fsts,boot.out,make.file=F,sig.list=names(boot.sig),pt.pch = 19)
  outliers<-find.outliers(fsts,boot.out)
  fsts$P.value<-boot.pvals
  fsts$Corr.P.value<-boot.cor.pvals
  fsts$Outlier<-FALSE
  fsts[fsts$Locus %in% outliers$Locus,"Outlier"]<-TRUE
  return(fsts)
}

fsthet<-function(gpop, fst.choice="fst",alpha=0.05){
  fsts<-calc.actual.fst(gpop,fst.choice)
  quant.out<-as.data.frame(t(replicate(1, fst.boot(gpop,fst.choice=fst.choice,ci=alpha,bootstrap=FALSE))))
  plotting.cis(fsts,boot.out=quant.out,make.file=F,pt.pch = 19)
  outliers<-find.outliers(fsts,quant.out)
  fsts$Outlier<-FALSE
  fsts[fsts$Locus %in% outliers$Locus,"Outlier"]<-TRUE
  return(fsts)
}

