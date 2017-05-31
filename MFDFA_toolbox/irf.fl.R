library(reshape)
irf.fl<-function(tmpvar,exon=NULL,k,ecdet="none",r=1)
{
 ncol(tmpvar)->nlevels
 tmpseq<-seq(1,nlevels,1)
 for (i in 1:nlevels){
   print(tmpseq)
   ca.jo(tmpvar[,tmpseq],K=2,dumvar=exon,ecdet=ecdet)->cajo.out1
   summary(cajo.out1)
   vec2var(cajo.out1,r=r)->var2
   irf.data.underscore(irf(var2,n.ahead=10,runs=250))->var2.irf
   subset(var2.irf,entry.ord==max(var2.irf$entry.ord))->last.irf
   subset(var2.irf,entry.ord==min(var2.irf$entry.ord))->first.irf
   ifelse (i==1, last.irf->irf.last, rbind(irf.last,last.irf)->irf.last)
   ifelse (i==1, first.irf->irf.first, rbind(irf.first,first.irf)->irf.first)  
   tmpseq[1]->buff
   tmpseq[2:nlevels]->tmpseq[1:nlevels-1]
   buff->tmpseq[nlevels]
 }
 rbind(irf.last,irf.first)->irf.fl	
 return(irf.fl)	
}


#This function is called above. Assumes data are named "n_3" as in MFDFA output
irf.data.underscore<- function(irf.out)
{
 as.data.frame(irf.out$irf)->irf.dat
 as.data.frame(irf.out$Upper)->Upper.dat
 as.data.frame(irf.out$Lower)->Lower.dat
 names(irf.out$irf)->irf.names
 colsplit(irf.names,split=("_"),names=c("pn","q"))->tmpd
 ifelse(tmpd$pn=="n",(as.numeric(tmpd$q)*-1),tmpd$q)->qvalues
 irf.dat$step<-as.numeric(rownames(irf.dat))
 out<-data.frame(1)
 count<-0
 for (imp in 1:length(names(irf.out$irf))){
  for (resp in 1:length(names(irf.out$irf))){
   	count<-count+1
	for (step in 1:length(unique(irf.dat$step))){
		irf.dat[step,count]->out$x
		Upper.dat[step,count]->out$ux
		Lower.dat[step,count]->out$lx
		#print(count)
		qvalues[imp]->out$qimp
		qvalues[resp]->out$qresp
		step->out$step
		imp->out$entry.ord
		ifelse (imp==1&&resp==1&&step==1, out->irf.out2, rbind(irf.out2,out)->irf.out2)
		}
	}
 }
 irf.out2<-irf.out2[order(irf.out2$qimp,irf.out2$qresp,irf.out2$step),]
 print(qvalues)
 return(irf.out2)
}