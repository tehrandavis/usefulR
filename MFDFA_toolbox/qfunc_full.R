#save("qfunc_full" ,file="/Users/jdixon/Dropbox/R_Stable/qfunc_full.RData")
#This function does the operations on the output from flucts
qfunc_full<-function (rsd, qlo=-4,qhi=4,qinc=.5,onplot=1) 
{
   #This function, authored by Damian, computes the legendre transform to get the f(alpha) spectrum
   legendre<-function(Hurstq,qvec){
      #first calc tauq
      tauq<-Hurstq*qvec-1
      #next calc alpha
      alphaq<-lsfit(qvec[1:2],tauq[1:2])$coef[2]
      for (i in 2:(length(qvec)-1)) {
	     #print(i)
	     alphaq<-rbind(alphaq,lsfit(qvec[(i-1):(i+1)],tauq[(i-1):(i+1)])$coef[2])
    	}
      alphaq<-rbind(alphaq,lsfit(qvec[(length(qvec)-1):length(qvec)],tauq[(length(tauq)-1):length(tauq)])$coef[2])

      #next calc f
     f<-qvec*alphaq-tauq
     as.data.frame(cbind(tauq,alphaq,f),row.names=seq(1:length(f)))->out
     names(out)<-c("tauq","alphaq","f")
     return(out)
    }
	
	#TRAP FOR ZERO RESIDUALS
	#Catching zero residuals values and replacing with smallest nonzero residual
	 min(rsd[rsd$resid>0,]$resid)->minresid
     ifelse(rsd$resid==0,minresid,rsd$resid)->rsd$resid
	
 cntq=0
 qmatlen<-((qhi-qlo)/qinc)+2 #extra column for var controlling record type: q, H, fit
 qmatrix<-matrix(NA,3,qmatlen) 
 if (onplot==1) par(mfrow=c(2,2), cex=.5)
 sc.list<-sort(unique(rsd$sc))
 rsd$resid<-rsd$resid/rsd$sc
 lsc<-log(sc.list)
 
 #Pulling the Chbadra and Jensen method from the dicty work. 
 #Hence the renaming of the input object
 rsd$resid->values
 rsd$sc->indnum
 data.frame(values,indnum)->anch.long
 aggregate(anch.long$values,by=list(anch.long$indnum),sum)->tmp
  names(tmp)<-c("indnum","sum")
  summary(tmp)
  merge(anch.long,tmp,by="indnum")->anch.long1
  summary(anch.long1)
  anch.long1$pmass<-anch.long1$values/anch.long1$sum

  summary(anch.long1)

##for testing
#qlo<--2
#qhi<-2
#qinc<-.5

  qpass<-0
  for (q in (seq(qlo,qhi, by=qinc))){
 	qpass<-qpass+1
	if (q != 1) {
	#print(q)
	#q<-(-4)
	#Calculate pmass raised to current q
	anch.long1$pmass.q<-anch.long1$pmass^q
	
	#Find the sum of the pmass.q, merge back to the original file
	aggregate(anch.long1$pmass.q,by=list(indnum=anch.long1$indnum),sum)->sum.pmass.q
	names(sum.pmass.q)[2]<-"sum.pmass.q"
	merge(anch.long1,sum.pmass.q,by="indnum")->anch.long2
	#rm(anch.long2)
	#plot(log(sum.pmass.q$indnum),log(sum.pmass.q$sum.pmass.q))
	#Compute Dq
	fit<-lm(log(sum.pmass.q$sum.pmass.q)~log(sum.pmass.q$indnum))
    summary(fit)$r.squared->fit.Dq
    fit$coefficients[[2]]/(q-1)->Dq
    #Compute the Muq
    anch.long2$muq<-anch.long2$pmass.q/anch.long2$sum.pmass.q
    #Compute the Muq.lnpmass
    anch.long2$muq.lnpmass<-anch.long2$muq*(log(anch.long2$pmass))
    
    #Find the sum of the muq.lnpmass, merge back to the original file
	aggregate(anch.long2$muq.lnpmass,by=list(indnum=anch.long2$indnum),sum)->sum.muq.lnpmass
	names(sum.muq.lnpmass)[2]<-"sum.muq.lnpmass"
	merge(anch.long2,sum.muq.lnpmass,by="indnum")->anch.long3
	#Compute aq
	fit<-lm(sum.muq.lnpmass$sum.muq.lnpmass~log(sum.muq.lnpmass$indnum))
    summary(fit)$r.squared->fit.aq
	fit$coefficients[[2]]->aq
	
	#Compute aq higher polynomials
	fit<-lm(sum.muq.lnpmass$sum.muq.lnpmass~log(sum.muq.lnpmass$indnum)+I(log(sum.muq.lnpmass$indnum)^2))
    summary(fit)$r.squared->fit.aq2
	summary(fit)[[4]][3]->aq2.coef
	summary(fit)[[4]][9]->aq2.t
	
	fit<-lm(sum.muq.lnpmass$sum.muq.lnpmass~log(sum.muq.lnpmass$indnum)+I(log(sum.muq.lnpmass$indnum)^2)+I(log(sum.muq.lnpmass$indnum)^3))
    summary(fit)$r.squared->fit.aq3
	summary(fit)[[4]][4]->aq3.coef
	summary(fit)[[4]][12]->aq3.t

	#Compute teq and mean tq
	table(anch.long3$indnum)[]->nboxes
	nboxes^(-1)->nboxes1
	anch.long3$pmass.q.1<-anch.long3$pmass^(q-1)
	aggregate(anch.long3$pmass.q.1,by=list(indnum=anch.long3$indnum),sum)->sum.pmass.q.1
	names(sum.pmass.q.1)[2]<-"sum.pmass.q.1"
	sum.pmass.q.1$sum.pmass.q.1*nboxes1->teq
	#mean tq
	fit<-lm(log(teq)~log(sum.pmass.q.1$indnum))
    summary(fit)$r.squared->fit.mtq
	fit$coefficients[[2]]->mtq
	
	#Compute FEQ and then faq
    anch.long3$muq.lnmuq<-anch.long3$muq*(log(anch.long3$muq)) 
	aggregate(anch.long3$muq.lnmuq,by=list(indnum=anch.long3$indnum),sum)->sum.muq.lnmuq
 	sum.muq.lnmuq
	names(sum.muq.lnmuq)[2]<-"sum.muq.lnmuq"
	
	#Compute faq
	fit<-lm(sum.muq.lnmuq$sum.muq.lnmuq~log(sum.muq.lnmuq$indnum))
    summary(fit)$r.squared->fit.faq
	fit$coefficients[[2]]->faq
	
    #Compute faq higher polynomials
	fit<-lm(sum.muq.lnmuq$sum.muq.lnmuq~log(sum.muq.lnmuq$indnum)+I(log(sum.muq.lnmuq$indnum)^2))
    summary(fit)$r.squared->fit.faq2
	summary(fit)[[4]][3]->faq2.coef
	summary(fit)[[4]][9]->faq2.t
	
	fit<-lm(sum.muq.lnmuq$sum.muq.lnmuq~log(sum.muq.lnmuq$indnum)+I(log(sum.muq.lnmuq$indnum)^2)+I(log(sum.muq.lnmuq$indnum)^3))
    summary(fit)$r.squared->fit.faq3
	summary(fit)[[4]][4]->faq3.coef
	summary(fit)[[4]][12]->faq3.t

	#Compute tq
	tq<-(q*aq)-faq
	#Alternative Dq computation. 
	#Dqv2<-tq/(q-1)
	
	#Now from the old qfunc code, pulling the hurst and fit
	  #Steps below are done at the onset
#	  sc.list<-sort(unique(rsd$sc))
#	  print(sc.list)
#      rsd$resid<-rsd$resid/rsd$sc
#      print(summary(rsd$resid))
      if (q==0) q<-.01
      avg.resid<-tapply((rsd$resid^I(q/2)), rsd$sc,mean) 
      avg.resid<-avg.resid^(1/q)
      lresid<-log(avg.resid)
      lsc<-log(sc.list)
      logtmp<-as.data.frame(cbind(lsc,lresid))
      fit<-lm(logtmp$lresid~logtmp$lsc)
      if (onplot==1) {
  	    plot(logtmp$lsc,logtmp$lresid,main=q)
  	    lines(logtmp$lsc,fitted(fit),type="b",col="red") 
  	  }
      H <-fit$coefficients["logtmp$lsc"]
      fit.H<-summary.lm(fit)$r.square
      #H with higher polynomials
      fit<-lm(logtmp$lresid~logtmp$lsc+I((logtmp$lsc)^2))
      summary(fit)$r.squared->fit.H2
	  summary(fit)[[4]][3]->H2.coef
	  summary(fit)[[4]][12]->H2.t
	  fit<-lm(logtmp$lresid~logtmp$lsc+I((logtmp$lsc)^2)+I((logtmp$lsc)^3))
      summary(fit)$r.squared->fit.H3
	  summary(fit)[[4]][4]->H3.coef
	  summary(fit)[[4]][12]->H3.t
  	} 

	curr.out<-data.frame(q,Dq,tq,aq,mtq,faq,fit.Dq,fit.aq,fit.mtq,fit.faq,fit.aq2,aq2.coef,aq2.t,fit.aq3,aq3.coef,aq3.t,fit.faq2,faq2.coef,faq2.t,fit.faq3,faq3.coef,faq3.t,
	H,fit.H,fit.H2,H2.coef,H2.t,fit.H3,H3.coef,H3.t)
    ifelse (qpass==1, curr.out->full.out, rbind(full.out,curr.out)->full.out)
  }
  #print(full.out)
  
  #Computes the legendre transform based on H and q to give the mf spectrum
  legendre(full.out$H,full.out$q)->l.out
  #print(l.out)
  l.out$tauq->full.out$tau.lg
  l.out$alphaq->full.out$aq.lg
  l.out$f->full.out$f.lg
  return(full.out)
} 
