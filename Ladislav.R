require(fracdiff)
require(fArma)
require(tseries)

#First-differencing
firstdiff<-function(x) {
  x-c(0,x[1:(length(x)-1)])
}

#Moving-block bootstrap for bivariate series
MBB<-function(x,y,bs){
  xm<-matrix(x,nrow=bs)
  xmb<-matrix(x,nrow=bs)
  ym<-matrix(y,nrow=bs)
  ymb<-matrix(y,nrow=bs)
  order<-sample(1:(length(x)/bs))
  for(i in (1:(length(x)/bs))){
    xmb[,i]<-xm[,order[i]]
    ymb[,i]<-ym[,order[i]]
  }
  return(matrix(c(xmb,ymb),ncol=2))
}

#Cross-periodogram coherency with $D$ Daniell's spans
CPS_coherency<-function(x,y,D){
  s<-spectrum(ts.union(ts(x),ts(y)),method="pgram",plot=FALSE,spans=D)
  return(matrix(c(sqrt(s$coh),s$freq*2*pi),ncol=2))	
}

#Cross-periodogram with $D$ Daniell's spans
CPS<-function(x,y,D){
  s<-spectrum(ts.union(ts(x),ts(y)),method="pgram",plot=FALSE,spans=D)
  return(matrix(c(sqrt(s$coh*spectrum(x,method="pgram",plot=FALSE,spans=D)$spec*spectrum(y,method="pgram",plot=FALSE,spans=D)$spec),s$freq*2*pi),ncol=2))
}

#Raw cross-periodogram
CPS_raw<-function(x,y){
  s<-spectrum(ts.union(ts(x),ts(y)),method="pgram",plot=FALSE)
  return(matrix(c(sqrt(s$coh*spectrum(x,method="pgram",plot=FALSE)$spec*spectrum(y,method="pgram",plot=FALSE)$spec),s$freq*2*pi),ncol=2))
}

#Aggregate cross-correlation test (giving test statistics and $p$-value)
ACC<-function(x,y,lag,brep,bs){
  CCF<-runif(brep+1)
  CCF[1]<-sum(abs(ccf(x,y,lag.max=lag,type="correlation",plot=FALSE)$acf))
  for(i in 1:brep){
    bb<-MBB(x,y,bs)
    CCF[i+1]<-sum(abs(ccf(bb[,1],bb[,2],lag.max=lag,type="correlation",plot=FALSE)$acf))
  }
  return(c(CCF[1],1-(rank(CCF)[1]/(brep+1))))
}

#Partial sums covariance divergence test (giving test statistics and $p$-value)
CD<-function(x,y,brep,bs){
  CD<-runif(brep+1)
  CD[1]<-cov(cumsum(x),cumsum(y))/length(x)	
  for(i in 1:brep){
    bb<-MBB(x,y,bs)
    CD[i+1]<-cov(cumsum(bb[,1]),cumsum(bb[,2]))/length(bb[,1])	
  }
  return(c(CD[1],1-2*abs(rank(CD)[1]/(brep+1)-0.5)))
}

#DCCA-based LRCC test (giving the test statistic, standard error and $p$-value)
DCCA_LRCC<-function(x,y,s){
  xx<-cumsum(x)
  yy<-cumsum(y)
  t<-1:length(xx)
  F2sj_xy<-runif(floor(length(xx)/s))
  F2sj_xx<-F2sj_xy
  F2sj_yy<-F2sj_xy
  for(ss in seq(1,(floor(length(xx)/s)*s),by=s)){
    F2sj_xy[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_xx[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_yy[(ss-1)/s+1]<-sum((summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
  }
  rho<-mean(F2sj_xy)/sqrt(mean(F2sj_xx)*mean(F2sj_yy))
  return(c(rho,1/sqrt(length(xx)),1-pnorm(rho,mean=0,sd=1/sqrt(length(xx)))))
}

#Rescaled covariance test (giving the test statistic)
RCT<-function(x,y,q,d1,d2){
  xx<-cumsum(x-mean(x))
  yy<-cumsum(y-mean(x))
  M<-(q)^(d1+d2)*cov(xx,yy)/(length(xx)*sum((1-abs(-q:q)/(q+1))*ccf(x,y,lag.max=q,type="covariance",plot=FALSE)$acf))
  return(M)
}

#Simulates correlated noise (of length $n$ and with correlation coefficient $\rho$)
Noise_rho<-function(n,rho){
  e1<-rnorm(n)
  e<-rnorm(n)
  e2<-rho*e1+sqrt(1-rho^2)*e	
  return(matrix(c(e1,e2),ncol=2))
}

#Two ARFIMA processes with correlated innovations (of length $n$, fractionally integrated series with $d1$ and $d2$ and with correlation coefficient of innovations $rho$)
ARFIMA_rho<-function(n,d1,d2,rho){
  e1<-rnorm(n)
  e<-rnorm(n)
  e2<-rho*e1+sqrt(1-rho^2)*e
  x<-fracdiff.sim(n,d=d1,innov=e1)$series
  y<-fracdiff.sim(n,d=d2,innov=e2)$series
  return(matrix(c(x,y),ncol=2))
}

#ARFIMA and AR(1) processes with correlated innovations (of length $n$, fractionally integrated series with $d1$, AR(1) parameter $\theta$ and with correlation coefficient of innovations $rho$)
ARFIMA_AR<-function(n,d1,theta,rho){
  e1<-rnorm(n)
  e<-rnorm(n)
  e2<-rho*e1+sqrt(1-rho^2)*e
  x<-fracdiff.sim(n,d=d1,innov=e1)$series
  y<-arima.sim(list(order=c(1,0,0),ar=theta),n,innov=e2)	
  return(matrix(c(x,y),ncol=2))	
}

#Two AR(1) processes with correlated innovations (of length $n$, with AR(1) parameters $\theta1$ and $\theta2$ and with correlation coefficient of innovations $rho$)
AR_rho<-function(n,theta1,theta2,rho){
  e1<-rnorm(n)
  e<-rnorm(n)
  e2<-rho*e1+sqrt(1-rho^2)*e
  x<-arima.sim(list(order=c(1,0,0),ar=theta1),n,innov=e1)
  y<-arima.sim(list(order=c(1,0,0),ar=theta2),n,innov=e2)	
  return(matrix(c(x,y),ncol=2))	
}

#Mixed-correlated ARFIMA processes - uncorrelated (of length $n$, fractionally integrated of order $d1$ and $d2$, with correlation coefficient of innovations $rho$)
Mixed_ARFIMA_noise_rho<-function(n,d1,d2,rho){
  e<-Noise_rho(n,rho)
  e2<-e[,1]
  e3<-e[,2]
  x<-fracdiff.sim(n,d=d1)$series+e2
  y<-e3+fracdiff.sim(n,d=d2)$series
  return(matrix(c(x,y),ncol=2))
}

#Mixed-correlated ARFIMA processes - short-range cross-correlated (of length $n$, fractionally integrated of order $d1$ and $d2$, short-range cross-correlation of AR(1) with $theta$ parameter and with correlation coefficient of innovations $rho$)
Mixed_ARFIMA_AR_rho<-function(n,d1,d2,theta,rho){
  e<-AR_rho(n,theta,theta,rho)
  e2<-e[,1]
  e3<-e[,2]
  x<-fracdiff.sim(n,d=d1)$series+e2
  y<-e3+fracdiff.sim(n,d=d2)$series
  return(matrix(c(x,y),ncol=2))
}

#Mixed-correlated ARFIMA processes - long-range cross-correlated (of length $n$, fractionally integrated of order $d1$ and $d4$ for the uncorrelated parts and $d2$ and $d3$ for the correlated parts, with correlation coefficient of innovations $rho$)
Mixed_ARFIMA_ARFIMA_rho<-function(n,d1,d2,d3,d4,rho){
  e<-ARFIMA_rho(n,d2,d3,rho)
  e2<-e[,1]
  e3<-e[,2]
  x<-fracdiff.sim(n,d=d1)$series+e2
  y<-e3+fracdiff.sim(n,d=d4)$series
  return(matrix(c(x,y),ncol=2))
}

#Fractionally cointegrated series with correlated innovations of $x_t$ and $u_t$ (of length $n$, parameter $beta$, memory $d_x$ of $x_t$ and $y_t$ and memory $d_U$ of the innovations, with correlation coefficient between $x_t$ and $u_t$ -- $rho$)
FCI<-function(n,beta,d_x,d_U,rho){
  e1<-rnorm(n)
  e<-rnorm(n)
  e2<-rho*e1+sqrt(1-rho^2)*e
  x<-fracdiff.sim(n,d=d_x,innov=e1)$series
  u<-fracdiff.sim(n,d=d_U,innov=e2)$series
  y<-beta*x+u
  return(matrix(c(x,y),ncol=2))
}

#Height cross-correlation analysis, both specifications
HXA_both<-function(x,y,tmin,tmax){
  options(warn=-1)
  tau<-rnorm(tmax-tmin+1)
  covxy<-rnorm(tmax-tmin+1)
  covxy_abs<-rnorm(tmax-tmin+1)
  for(nu in tmin:tmax){
    t<-1:ceiling(length(x)/nu)
    xx<-summary(lm(c(matrix(cumsum(x),nrow=nu)[1,])~t))$res
    yy<-summary(lm(c(matrix(cumsum(y),nrow=nu)[1,])~t))$res
    cov<-(xx[2:length(xx)]-xx[1:(length(xx)-1)])*(yy[2:length(yy)]-yy[1:(length(yy)-1)])
    covxy[(nu-tmin+1)]<-mean(cov,na.rm=TRUE)
    covxy_abs[(nu-tmin+1)]<-mean(abs(cov))
    tau[(nu-tmin+1)]<-nu
  }
  H<-summary(lm(log10(covxy)~log10(tau)))$coefficients[2]/2
  H_abs<-summary(lm(log10(covxy_abs)~log10(tau)))$coefficients[2]/2
  return(c(H,H_abs))
}

#HXA jackknife estimator
HXA_jackknife_both<-function(x,y,tmax1,tmax2){
  HH<-rnorm(tmax2-tmax1+1)
  HH_abs<-rnorm(tmax2-tmax1+1)
  for(i in tmax1:tmax2){
    HH[i-tmax1+1]<-HXA_both(x,y,1,i)[1]
    HH_abs[i-tmax1+1]<-HXA_both(x,y,1,i)[2]
  }
  return(c(mean(HH),mean(HH_abs)))
}

#HXA-based estimator of power law coherency
HXA_rho<-function(x,y,tmin,tmax){
  options(warn=-1)
  tau<-rnorm(tmax-tmin+1)
  covxy<-rnorm(tmax-tmin+1)
  varxx<-rnorm(tmax-tmin+1)
  varyy<-rnorm(tmax-tmin+1)
  for(nu in tmin:tmax){
    t<-1:ceiling(length(x)/nu)
    xx<-summary(lm(c(matrix(cumsum(x),nrow=nu)[1,])~t))$res
    yy<-summary(lm(c(matrix(cumsum(y),nrow=nu)[1,])~t))$res
    cov<-(xx[2:length(xx)]-xx[1:(length(xx)-1)])*(yy[2:length(yy)]-yy[1:(length(yy)-1)])
    var_x<-(xx[2:length(xx)]-xx[1:(length(xx)-1)])^2
    var_y<-(yy[2:length(xx)]-yy[1:(length(xx)-1)])^2		
    covxy[(nu-tmin+1)]<-mean(cov,na.rm=TRUE)
    varxx[(nu-tmin+1)]<-mean(var_x)
    varyy[(nu-tmin+1)]<-mean(var_y)
    tau[(nu-tmin+1)]<-nu
  }
  H<-summary(lm(log10(covxy^2/(varxx*varyy))~log10(tau)))$coefficients[2]/4
  return(H)
}

#Detrended cross-correlation analysis, both specifications with step between used scales of $sstep$ between $smin$ and $smax$
DCCA_both_nonoverlap_step<-function(x,y,smin,smax,sstep){
  options(warn=-1)
  xx<-cumsum(x)
  yy<-cumsum(y)
  t<-1:length(xx)
  F2s<-rnorm((smax-smin)/sstep+1)
  F2s_abs<-rnorm((smax-smin)/sstep+1)
  sss<-rnorm((smax-smin)/sstep+1)
  for(s in seq(smin,smax,by=sstep)){
    F2sj<-rnorm(floor(length(xx)/s))
    F2sj_abs<-rnorm(floor(length(xx)/s))
    for(ss in seq(1,(floor(length(xx)/s)*s),by=s)){
      FF<-(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)
      F2sj[(ss-1)/s+1]<-sum(FF)/(s-1)
      F2sj_abs[(ss-1)/s+1]<-sum(abs(FF))/(s-1)
    }
    F2s[(s-smin)/sstep+1]<-mean(F2sj,na.rm=TRUE)
    F2s_abs[(s-smin)/sstep+1]<-mean(F2sj_abs)
    sss[(s-smin)/sstep+1]<-s
  }
  H_DCCA<-summary(lm(log10(F2s)~log10(sss)))$coefficients[2]/2	
  H_DCCA_abs<-summary(lm(log10(F2s_abs)~log10(sss)))$coefficients[2]/2	
  return(c(H_DCCA,H_DCCA_abs))
}

#DCCA-based estimator of the power law coherency
DCCA_rho<-function(x,y,smin,smax,sstep){
  options(warn=-1)
  xx<-cumsum(x)
  yy<-cumsum(y)
  t<-1:length(xx)
  F2s_xy<-rnorm((smax-smin)/sstep+1)
  F2s_xx<-rnorm((smax-smin)/sstep+1)
  F2s_yy<-rnorm((smax-smin)/sstep+1)
  sss<-rnorm((smax-smin)/sstep+1)
  for(s in seq(smin,smax,by=sstep)){
    F2sj_xy<-rnorm(floor(length(xx)/s))
    F2sj_xx<-rnorm(floor(length(xx)/s))
    F2sj_yy<-rnorm(floor(length(xx)/s))
    for(ss in seq(1,(floor(length(xx)/s)*s),by=s)){
      FF_xy<-(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)
      FF_xx<-(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)^2
      FF_yy<-(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)^2
      F2sj_xy[(ss-1)/s+1]<-sum(FF_xy)/(s-1)
      F2sj_xx[(ss-1)/s+1]<-sum(FF_xx)/(s-1)
      F2sj_yy[(ss-1)/s+1]<-sum(FF_yy)/(s-1)
      
    }
    F2s_xy[(s-smin)/sstep+1]<-mean(F2sj_xy)
    F2s_xx[(s-smin)/sstep+1]<-mean(F2sj_xx)
    F2s_yy[(s-smin)/sstep+1]<-mean(F2sj_yy)
    sss[(s-smin)/sstep+1]<-s
  }
  H_DCCA_abs<-summary(lm(log10(F2s_xy^2/(F2s_xx*F2s_yy))~log10(sss)))$coefficients[2]/4	  
  return(H_DCCA_abs)
}

#CCDMA, both specifications, centered MA, use odd $kmin$ and $kmax$
CCDMA_both<-function(x,y,kmin,kmax){
  options(warn=-1)
  xx<-cumsum(x)
  yy<-cumsum(y)
  F2k<-rnorm((kmax-kmin)/2+1)
  F2k_abs<-rnorm((kmax-kmin)/2+1)
  kkk<-seq(kmin,kmax,by=2)
  for(k in seq(kmin,kmax,by=2)){
    ma<-c(rep(1,k))/k
    X_MA<-filter(xx,ma)
    Y_MA<-filter(yy,ma)
    F2<-(xx-X_MA)[(1+floor(k/2)):(length(xx)-floor(k/2))]*(yy-Y_MA)[(1+floor(k/2)):(length(yy)-floor(k/2))]
    F2k[floor(k/2)]<-mean(F2,na.rm=TRUE)
    F2k_abs[floor(k/2)]<-mean(abs(F2))		
  }
  H<-summary(lm(log10(F2k)~log10(kkk)))$coefficients[2]/2	
  H_abs<-summary(lm(log10(F2k_abs)~log10(kkk)))$coefficients[2]/2	
  return(c(H,H_abs))	
}

#CCDMA-based estimator of the power law coherency
CCDMA_rho<-function(x,y,kmin,kmax){
  options(warn=-1)
  xx<-cumsum(x)
  yy<-cumsum(y)
  F2k_xy<-rnorm((kmax-kmin)/2+1)
  F2k_xx<-rnorm((kmax-kmin)/2+1)
  F2k_yy<-rnorm((kmax-kmin)/2+1)
  kkk<-seq(kmin,kmax,by=2)
  for(k in seq(kmin,kmax,by=2)){
    ma<-c(rep(1,k))/k
    X_MA<-filter(xx,ma)
    Y_MA<-filter(yy,ma)
    F2_xy<-(xx-X_MA)[(1+floor(k/2)):(length(xx)-floor(k/2))]*(yy-Y_MA)[(1+floor(k/2)):(length(yy)-floor(k/2))]
    F2_xx<-(xx-X_MA)[(1+floor(k/2)):(length(xx)-floor(k/2))]^2
    F2_yy<-(yy-Y_MA)[(1+floor(k/2)):(length(yy)-floor(k/2))]^2
    F2k_xy[floor(k/2)]<-mean(F2_xy)
    F2k_xx[floor(k/2)]<-mean(F2_xx)
    F2k_yy[floor(k/2)]<-mean(F2_yy)
  }
  H<-summary(lm(log10(F2k_xy^2/(F2k_xx*F2k_yy))~log10(kkk)))$coefficients[2]/4	
  return(H)	
}

#Cross-periodogram estimator (raw), using up to 2*pi*$part$ frequency
CPer<-function(x,y,part){
  options(warn=-1)
  C<-CPS_raw(x,y)
  lambda<-C[,2]
  f<-C[,1]
  H<-((-summary(lm(log10(f[1:(floor(length(x)*part))])~log10(lambda[1:(floor(length(x)*part))])))$coefficients[2])+1)/2
  return(H)
}

#Cross-periodogram estimator (using $Daniell_m$ spans), using up to 2*pi*$part$ frequency
CPer_smooth<-function(x,y,part,Daniell_m){
  options(warn=-1)
  C<-CPS(x,y,Daniell_m)
  lambda<-C[,2]
  f<-C[,1]
  H<-((-summary(lm(log10(f[1:(floor(length(x)*part))])~log10(lambda[1:(floor(length(x)*part))])))$coefficients[2])+1)/2
  return(H)
}

#XPE-based estimator of the power law coherency (using $Daniell_m$ spans), using up to 2*pi*$part$ frequency
CPer_smooth_rho<-function(x,y,part,Daniell_m){
  options(warn=-1)
  C_xy<-CPS(x,y,Daniell_m)
  lambda<-C_xy[,2]
  f_xy<-C_xy[,1]
  C_x<-CPS(x,x,Daniell_m)
  C_y<-CPS(y,y,Daniell_m)
  f_x<-C_x[,1]
  f_y<-C_y[,1]
  H_xy<-((-summary(lm(log10(f_xy[1:(floor(length(x)*part))])~log10(lambda[1:(floor(length(x)*part))])))$coefficients[2])+1)/2
  H_xx<-((-summary(lm(log10(f_x[1:(floor(length(x)*part))])~log10(lambda[1:(floor(length(x)*part))])))$coefficients[2])+1)/2
  H_yy<-((-summary(lm(log10(f_y[1:(floor(length(x)*part))])~log10(lambda[1:(floor(length(x)*part))])))$coefficients[2])+1)/2
  return(H_xy-(H_xx+H_xy)/2)
}

#Averaged periodogram estimator using $q$ (suggested to use q=0.5) and using up to 2*pi*$part$ frequency
APE<-function(x,y,part,q){
  options(warn=-1)
  C<-CPS_raw(x,y)
  f<-C[,1]	
  lambda<-C[,2]
  H<-1-log10(sum(f[1:(floor(length(x)*part*q))])/sum(f[1:(floor(length(x)*part))]))/(2*log10(q))
  return(H)
}

#Averaged periodogram estimator using $q$ (suggested to use q=0.5), with $D$ smoothing and using up to 2*pi*$part$ frequency
APE_smooth<-function(x,y,part,q,D){
  options(warn=-1)
  C<-CPS(x,y,D)
  f<-C[,1]	
  lambda<-C[,2]
  H<-1-log10(sum(f[1:(floor(length(x)*part*q))])/sum(f[1:(floor(length(x)*part))]))/(2*log10(q))
  return(H)
}

#APE-based estimator of the power law coherency using $q$ (suggested to use q=0.5) and using up to 2*pi*$part$ frequency
APE_smooth_rho<-function(x,y,part,q,D){
  options(warn=-1)
  C_xy<-CPS(x,y,D)
  f_xy<-C_xy[,1]	
  f_x<-CPS(x,x,D)[,1]
  f_y<-CPS(y,y,D)[,1]
  lambda<-C_xy[,2]
  H_rho<-(1-log10(sum(f_xy[1:(floor(length(x)*part*q))])/sum(f_xy[1:(floor(length(x)*part))]))/(2*log10(q)))-(1-(log10(sum(f_x[1:(floor(length(x)*part*q))])/sum(f_x[1:(floor(length(x)*part))]))/(2*log10(q)))+1-log10(sum(f_y[1:(floor(length(x)*part*q))])/sum(f_y[1:(floor(length(x)*part))]))/(2*log10(q)))/2
  return(H_rho)
}

#Local X-Whittle with using up to 2*pi*$part$ frequency
LXW<-function(x,y,part){
  C<-CPS_raw(x,y)
  lambda<-C[,2]
  f<-C[,1]
  R<-function(i) (log(mean(((lambda)^(2*i-1))[1:(floor(length(x)*part))]*f[1:(floor(length(x)*part))]))-(2*i-1)*mean(log(lambda)[1:(floor(length(x)*part))]))
  H<-nlm(R,0.5)$estimate
  return(H)
}

#Local X-Whittle with using up to 2*pi*$part$ frequency using $Daniell_m$ smoothing spans
LXW_smooth<-function(x,y,part,Daniell_m){
  C<-CPS(x,y,Daniell_m)
  lambda<-C[,2]
  f<-C[,1]
  R<-function(i) (log(mean(((lambda)^(2*i-1))[1:(floor(length(x)*part))]*f[1:(floor(length(x)*part))]))-(2*i-1)*mean(log(lambda)[1:(floor(length(x)*part))]))
  H<-nlm(R,0.5)$estimate
  return(H)
}

#LXW-based estimator of the power law coherency
LXW_coh<-function(x,y,part,Daniell_m){
  C<-CPS(x,y,Daniell_m)
  lambda<-C[,2]
  f_xy<-C[,1]
  C_x<-CPS(x,x,Daniell_m)
  C_y<-CPS(y,y,Daniell_m)
  f_x<-C_x[,1]
  f_y<-C_y[,1]
  coh<-f_xy^2/(f_x*f_y)
  R<-function(i) (log(mean(((lambda)^(4*i))[1:(floor(length(x)*part))]*coh[1:(floor(length(x)*part))]))-(4*i)*mean(log(lambda)[1:(floor(length(x)*part))]))
  rho<-nlm(R,0)$estimate 
  return(rho)
}