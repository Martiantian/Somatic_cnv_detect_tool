args=commandArgs(TRUE);
a=read.table(args[1])
mf=median(a[,3])
a=a[a[,1]<=22,]
a=a[(a[,3]>=mf*0.7 & a[,3]<=mf*1.4),]
if(dim(a)[2]>3){  
	 	b=a[,-c(1:2)]
		sd_a=apply(b,1,sd) 
		c=apply(b,1,median) 
		w=cbind(a[,c(1,2)],c)
	rs=which(sd_a>=mean(sd_a)+2*sd(sd_a))
	w=w[-rs,] 
	w=w[(w[,3]>0),]
}else{
	w=a
}

write.table(w,args[2],col.names=FALSE,row.names=FALSE,sep="\t")
q()


