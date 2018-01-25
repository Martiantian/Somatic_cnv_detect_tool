args=commandArgs(TRUE);
alpha=as.numeric(args[4])
threshold=as.numeric(args[2])
binlength=as.numeric(args[6])
qn=qnorm(1-alpha/2)

detect_line=function(ww){
	a=ww[ww[,dim(ww)[2]]!="Normal",]
	if(dim(a)[1]!=0){
		W=NULL
		for(i in sort(unique(a[,1]))){
			b=a[a[,1]==i,]
			chr=b[1,1]
			start=b[1,2]
			end=start
			ye=b[1,4]
			nn=b[1,5]
			p=b[1,6]
			cnv=as.character(b[1,7])

			for(j in 1:dim(b)[1]){
				if(b[j,4]==ye){
					end=b[j,2]
				}else{
					tmp=c(chr,start*binlength,(end+1)*binlength,binlength*(end-start+1),ye,p,cnv,start,end,nn)
					W=rbind(W,tmp)
					start=b[j,2]
					end=start
					ye=b[j,4]
					nn=b[j,5]
					p=b[j,6]
					cnv=as.character(b[j,7])
				}
			}

			tmp=c(chr,start*binlength,(end+1)*binlength-1,binlength*(end-start+1),ye,p,cnv,start,end,nn)
			W=rbind(W,tmp)
		}
	}else{
		W=a
	}
	return(W)
}



checkplot=function(a){
	x=a[a[,1]<=22 & a[,1]!=19,]
	xx=abs(x[,4]-x[,3]+1)
	mean=mean(xx)
	sd=sd(xx)
	W=NULL
	P=NULL
	P1=NULL
	for(i in sort(unique(a[,1]))){
		b=a[a[,1]==i,]
		if(dim(b)[1]!=0){
			y1=b[,3]
			ye1=b[,4]
			tmpn=b[,5]
			dp=dl=ni=NULL
			dpp=dll=NULL
			for(j in 1:length(ye1)){
				if(! is.na(ye1[j])){
					if(i==19){
						bp=2*threshold
					}else{
						bp=threshold
					}
					if(ye1[j]>=mean){
						p=1-pnorm(ye1[j],mean*(1+bp),sd/sqrt(tmpn[j]))
					}else{
					p=pnorm(ye1[j],mean*(1-bp),sd/sqrt(tmpn[j]))
					}
				}else{
				p=1
				}
				p=signif(p,4)
				P=c(P,p)
				p1="Normal"
				if(ye1[j]>mean*(1+threshold) & p<alpha/2){
					p1="Dup"
				}
				if(ye1[j]<mean*(1-threshold) & p<alpha/2){
					p1="Del"
				}
				P1=c(P1,p1)
			}
		}
	}
	ww=cbind(a,P,as.matrix(P1)) 
	return(ww)
}
a=read.table(args[5])
ww=checkplot(a)
write.table(ww,args[7],quote=F,col.names=F,row.names=F,sep="\t")
s<-ww[ww[,7]=="Normal",]
sdr<-sd(s[,3])
write.table(sdr,args[9],sep="\t",quote=F,col.names=F,row.names=args[1],append=T)

w=detect_line(ww)
if(dim(w)[1]!=0){
	nname=rep(args[1],dim(w)[1])
	w=cbind(nname,w)
}
w.names=c("Sample","Chr","Start","End","Length","Estimate_copy_ratio","P_vlaue","CNV","Start_bin","End_bin","True_bin_number")
if(args[10]==1){
	w=rbind(w.names,w)
}
write.table(w,args[8],quote=F,col.names=F,row.names=F,sep="\t",append=T)
q()