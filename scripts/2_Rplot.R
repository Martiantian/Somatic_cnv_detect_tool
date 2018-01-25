#Author: ***
#Email:***@genomics.org.cn
#File Name:
#Description:
# this program after extractInfo.pl chrno start length
#Edit History:
#2016-01-30 17:00:03  File created.
args<-commandArgs(TRUE)
	dir=paste(args[2],"tiff",sep=".")
	name=args[3]
	tiff(dir,width=2400,height=1800,res=600,units = "px",compression="lzw")
	a=read.table(args[1])
	a=a[a[,1]>=200 & a[,1]<=700,]
	x=a[,1]/10
	y=a[,2]/a[,3]
	par(mar=c(4.5,4.5,1.5,1.5))
	ymax=max(y)
	plot(0,cex=0,axes=F,xlab="GC(%)",xlim=c(15,75),ylim=c(0,ymax+0.005),ylab="RCr",main=name)
	axis(1)
	axis(2)
	lines(x,y,lwd=1.5)
	dev.off()
