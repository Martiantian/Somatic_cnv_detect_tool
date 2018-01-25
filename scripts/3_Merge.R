args=commandArgs(TRUE);
flist=read.table(args[1]);
		files=as.vector(flist[,1])
		n=length(files)

		names=NULL
		for(i in 1:n){
			a=read.table(files[i]);
			names=unique(c(names,as.matrix(paste(a[,1],a[,2],sep="\t")))) 
		}
			
		seq=names

		m=matrix(0,nrow=length(seq),ncol=n) 
		rownames(m)=seq
		for(i in 1:n){
			a=read.table(files[i]);
			rownames(a)=as.matrix(paste(a[,1],a[,2],sep="\t"))
			##filter the bin of large changes for one sample
			mf=median(a[,(dim(a)[2])])
			a=a[(a[,(dim(a)[2])]>=mf*0.7 & a[,(dim(a)[2])]<=mf*1.4),]
			m[as.matrix(paste(a[,1],a[,2],sep="\t")),i]=a[,(dim(a)[2])]
		}

		w=cbind(seq,m)
		write.table(w,args[2],sep="\t",quote=F,col.names=F,row.names=F)
		q()



