#Author: wangsxiaofeng
#Email:wangxiaofeng@genomics.cn
#File Name:
#Description:
# 
#Edit History:
#2016-01-30 17:00:03  File created.
args=commandArgs(TRUE);
	a=read.table(args[1])
	n=dim(a)[2]
	b=a[a[,1]<=22,n] 
	y=median(b)  

	write.table(y,"normalized_nb",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
q()
