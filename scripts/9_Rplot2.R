#!/bin/perl -w
#Author: wangxiaofeng
#Email:wangxaiofeng@genomics.cn
#File Name:
#Description:
#	
#Edit History:
#2016-01-30 17:00:03	File created.
args<-commandArgs(TRUE)
	a<-read.table(args[1])
	a<-a[a[,1]<=22,]
	a[,3:4]<-a[,3:4]*2
	a[,8]<-rep(3,length(a[,1]))
	a[a[,7]=="Dup",8]<-2
	a[a[,7]=="Del",8]<-4
	for (i in 2:22) {
		a[a[,1]==i,2]=a[a[,1]==i,2]+max(a[a[,1]==i-1,2])+10
	}
	b<-c(1:22)
	b[1]<-0
	for (i in 1:22) {
		b[i+1]=max(a[a[,1]==i,2])+5
	}
	d<-c(1:22)
	d[1]=125
	for (i in 2:22) {
		d[i]=b[i]+(b[i+1]-b[i]-10)/2
	}
	ymax<-max(a[,3])
	if (ymax>4){
		yt=ymax
	}else{
		yt=4
	}
	name=paste(args[2],"cnv","tiff",sep=".")
	tiff(name,width =16000,height =4000,units = "px",compression = "lzw",res=600)
	par(mai=c(1.5,1.5,1,1))
	plot(a[,2],a[,3],pch=20,col=a[,8],cex=0.6,cex.axis=2,mgp=c(4.5,1.3,0),cex.lab=2,
		xlab="Chromosome",ylab="Copy Number",xlim =c(0,max(b)),ylim=c(0,yt),xaxt="n",main=args[3])
	abline(v=b,lty=2,col="grey")
	abline(h=2,lty=2,col="grey")
	lines(a[,2],a[,4],type="l",lwd=1)
	legend("topleft",legend=c("Normal","Duplication","Deletion"),col=c("green","red","blue"),cex=1.3,pch=20)
	axis(1,at=d,cex.lab=2,cex.axis=2,mgp=c(4.5,1.3,0),labels =c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
	dev.off()
