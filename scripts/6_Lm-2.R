#!/bin/perl -w
#Author: qianzhaoyang
#Email:qianzhaoyang@genomics.org.cn
#File Name:
#Description:

args<-commandArgs(TRUE)
a=read.table(args[1])
a=a[a[,1]<=22,]
k=ncol(a)
c=k
d=k-1
a=a[a[,c]>0.3,]
a=a[a[,c]<0.7,]
e=a[,c]^2
z=lm(a[,d]~1+a[,c]+e)
median(a[,d])
z$coefficients
