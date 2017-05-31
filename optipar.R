optiPar<-function (ts1, ts2, par, min.rec = 2, max.rec = 5) 
{
  maxEmbed = 5
  
  rescale = normalize = mindiagline = minvertline = lgM = fnnpercent = recpt = whiteline = tw = steps = radiusspan = radiussample = NULL
  if (is.vector(ts1) != TRUE) {
    ts1 = as.numeric(as.matrix(ts1))
  }
  if (is.vector(ts2) != TRUE) {
    ts2 = as.numeric(as.matrix(ts2))
  }
  for (v in 1:length(par)) assign(names(par)[v], par[[v]])
  if (radiusspan <= 1) {
    stop("Radius span too small, please choose a larger value")
  }
  
  mi1 = mutual(ts1, lag.max = lgM, plot = FALSE)
  mi2 = mutual(ts2, lag.max = lgM, plot = FALSE)
  
  ## iterate over a range of possible steps to find out the optimal one
  
  lag1 = lag2 = vector()
  for (s in steps){
    # print(s)
    
    for (l in 1:(lgM - s)){
      if (mi1[l] < mi1[l + s]){
        lg1 = which(mi1 == mi1[l])
        lag1 = c(lag1, lg1)
        break  }
    }
    
    for (l in 1:(lgM - s)){
      if (mi2[l] < mi2[l + s]){
        lg2 = which(mi2 == mi2[l])
        lag2 = c(lag2, lg2)
        break  }
    }
  }
  
  if (length(lag1) == 0 | length(lag2) == 0){
    stop("Please try varying maximum lag (and/or step size):
         minimal mutual information was not reached")
    del = 0
  }
  
  ## take unique lags, i.e., avoid repeat
  lag1 = unique(lag1); lag2 = unique(lag2)
  
  
  ## for each time-series take the lag where MI is minimal
  
  ix.lg1 = which( mi1[lag1] == min(mi1[lag1]) )
  if (length(ix.lg1) > 1){
    ix.lg1 = sample(ix.lg1, 1) ## if there are more than one minimum
    ## sample one 
    lg1 = lag1[ix.lg1]}
  
  ix.lg2 = which( mi2[lag2] == min(mi2[lag2]) )
  if (length(ix.lg2) > 1){
    ix.lg2 = sample(ix.lg2, 1)
    lg2 = lag2[ix.lg2]}
  
  del = as.numeric(max(lg1,lg2))
  if (length(del) > 1) {
    del = del[length(del)]
  }
  del = del - 1
  
  embdts1 = false.nearest(ts1, m = l, d = del, t = 0, rt = 10, 
                          eps = sd(ts1)/10)
  fnnfraction1 = embdts1[1, ]
  fnnfraction1 = fnnfraction1[which(is.na(fnnfraction1) == 
                                      FALSE)]
  emdthd1 = fnnfraction1[1]/fnnpercent
  emdix1 = which(diff(fnnfraction1) < -emdthd1)
  if (length(emdix1) == 1) {
    emdmints1 = as.numeric(emdix1) + 1
  }
  else if (length(emdix1) > 1) {
    emdmints1 = as.numeric(tail(emdix1, 1) + 1)
  }
  else {
    emdmints1 = 1
  }
  embdts2 = false.nearest(ts2, m = maxEmbed, d = del, t = 0, rt = 10, 
                          eps = sd(ts2)/10)
  fnnfraction2 = embdts2[1, ]
  fnnfraction2 = fnnfraction2[which(is.na(fnnfraction2) == 
                                      FALSE)]
  emdthd2 = fnnfraction2[1]/fnnpercent
  emdix2 = which(diff(fnnfraction2) < -emdthd2)
  if (length(emdix2) == 1) {
    emdmints2 = as.numeric(emdix2) + 1
  }
  else if (length(emdix2) > 1) {
    emdmints2 = as.numeric(tail(emdix2, 1) + 1)
  }
  else {
    emdmints2 = 1
  }
  if (length(emdmints1) > 1) {
    emdmints1 = emdmints1[1]
  }
  if (length(emdmints2) > 1) {
    emdmints2 = emdmints2[1]
  }
  embdim = max(c(emdmints1, emdmints2))
  dm = rdist(ts1, ts2)
  if (rescale > 0) {
    switch(rescale, {
      1
      rescaledist = mean(dm)
      dmrescale = (dm/rescaledist) * 100
    }, {
      2
      rescaledist = max(dm)
      dmrescale = (dm/rescaledist) * 100
    })
  }
  else {
    dmrescale = dm
  }
  combo = c(ts1, ts2)
  sdun = sd(dmrescale)
  mnun = median(dmrescale) * 2
  radi = seq(mnun, 0, -(sdun/radiusspan))
  radi = radi[(length(radi)/2):length(radi)]
  kpt = ceiling(length(radi)/radiussample)
  rsamples = sample(1:kpt, 1)
  syssamp = seq(rsamples, rsamples + kpt * (radiussample - 
                                              1), kpt)
  syssamp = syssamp[syssamp <= length(radi)]
  radi = radi[syssamp]
  delay = del
  embed = embdim
  optrad = vector()
  end.flag <- 0
  while (end.flag == 0) {
    hi.loc <- 1
    lo.loc <- length(radi)
    curr.loc <- round(length(radi)/2)
    r <- radi[curr.loc]
    radi
    res = crqa(ts1, ts2, delay, embed, rescale, r, normalize, 
               mindiagline, minvertline, tw, whiteline, recpt)
    if (res$RR >= min.rec & res$RR <= max.rec) {
      optrad = r
      end.flag <- 1
    }
    else {
      if (res$RR < min.rec) {
        lo.loc <- curr.loc
      }
      if (res$RR > max.rec) {
        hi.loc <- curr.loc
      }
      if ((lo.loc - hi.loc) < 2) {
        end.flag <- 1
        warning("Optimal Radius Not found: try again choosing a wider radius span and larger sample size")
        return(list(radius = 0, emddim = 0, delay = 0))
      }
    }
    radi <- radi[hi.loc:lo.loc]
  }
  if (length(optrad) == 0) {
    optrad = NA
    return(list(radius = 0, emddim = 0, delay = 0))
  }
  if (!is.na(optrad)) {
    return(list(radius = optrad, emddim = embdim, delay = del))
  }
}


