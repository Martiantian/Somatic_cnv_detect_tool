args=commandArgs(TRUE);
	alpha=as.numeric(args[2])
	name=paste(args[4],"cv0.normal.txt",sep=".")

	qn=qnorm(1-alpha/2)

	normalize=function(a){
	x=a[a[,1]<=22,]
	xx=abs(x[,4]-x[,3]+1)
	mean=mean(xx)
	sd=sd(xx)
	dp=dl=ni=NULL
	ni_nXY=NULL
	for(i in 1:dim(a)[1]){
		if(! is.na(a[i,4])){
			if(a[i,4] > (mean*1.01+qn*sd/sqrt(a[i,5])) ){
				dp=c(dp,i)
			}else if(a[i,4] < (mean*0.99-qn*sd/sqrt(a[i,5])) ){
				dl=c(dl,i)
			}else{
				ni=c(ni,i)
				if(a[i,1]<=22){
					ni_nXY=c(ni_nXY,i)
				}
			}
		}
	}
##############calculate the cv0 for all chr exact chrX and chrY by the determined normal bin data#####################
		cv=sd(a[ni_nXY,3])/mean(a[ni_nXY,3])
		write.table(cv,name,quote=F,col.names=F,row.names=F,sep="\t")
######################################################################################################################
		normalized_y=median(a[ni,3])
		a[,3]=a[,3]/normalized_y
		return(a[1:3])
	}
	a=read.table(args[1])
	W=normalize(a)
	write.table(W,args[3],sep="\t",quote=F,col.names=F,row.names=F)
q()
