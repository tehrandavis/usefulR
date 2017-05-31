#library(fractal)
#library(multicore)
#save("flucts.mc.distrib" ,file="/Users/jdixon_old/Documents/R_Stuff/flucts.mc.distrib.RData")
#save("qfunc",file="qfunc.RData")
#
#    x<-rnorm(3000)
#	detrend<-"poly1"
#	sum.order<-1
#	overlap<-.9
#	scale.max<-trunc(length(x)/2)
#	scale.min<-NULL
#	scale.ratio<-.4
#	verbose<-FALSE
	
flucts.mc.distrib<-function (x, detrend = "poly1", sum.order = 0, overlap = 0, scale.max = trunc(length(x)/2), 
    scale.min = NULL, scale.ratio = 2, verbose = FALSE,framesize = 50) {  
    #This function builds a list of all binsizes and starting locations (for the time series)   
    #Each of which requires a regression to estimate the residual. It then farms out the elements of the list to the multiple cores
    # The framesize parameter controls how many list elements each core gets. 
    		
	array(data=NA,c(1,1))->tmp.array
    as.data.frame(tmp.array)->tmp
    polyfit.model <- function(polyfit.order) {
        if (polyfit.order < 0) 
            stop("Polynomial fit order must be positive")
        if (polyfit.order == 0) 
            return("x ~ 1")
        if (polyfit.order == 1) 
            return("x ~ 1 + t")
        else {
            poly.str <- paste("x ~ 1 + t +", paste("t^", seq(2, 
                polyfit.order), "+", sep = "", collapse = ""), 
                sep = "", collapse = "")
            return(substring(poly.str, 1, nchar(poly.str) - 1))
        }
    }
    regression.poly <- function(x, model = polyfit.model(1)) {
        #order <- length(attr(terms(formula(x ~ 1)), "order"))
        #print(order)
        #t <- x@positions
        x <- x@data
        t<-seq(length(x))
        if (polyfit.order == 0) return(sum((x - mean(x))^2))
        fit <- lm(model, data = data.frame(list(t = t, x = x)))
        return(sum(fit$residuals^2))
        #plot(t,x,main=q)
        #lines(xx,fitted(fit),type="b",col="red")
    }
    regression.bridge <- function(x, ...) {
        x <- x@data
        N <- length(x)
        bridge <- seq(x[1], x[N], length = N)
        x <- x - bridge
        return(sum(x^2))
    }
    regression.none <- function(x, ...) return(sum(x@data^2))
    checkScalarType(detrend, "character")
    checkScalarType(scale.ratio, "numeric")
    checkScalarType(verbose, "logical")
    checkScalarType(overlap, "numeric")
    data.name <- deparse(substitute(x))
    x <- create.signalSeries(x)
    if ((overlap < 0) | (overlap >= 1)) 
        stop("Overlap factor must be in the range [0,1)")
    detrend <- lowerCase(detrend)
    if (substring(detrend, 1, 4) == "poly") {
        polyfit.order <- as.integer(substring(detrend, 5))
        if (is.na(polyfit.order)) {
            stop("Improperly formed detrending s.tring")
        }
        else if (polyfit.order < 0) {
            stop("Polynomial fit order must be positive")
        }
        model <- formula(polyfit.model(polyfit.order))
        if (verbose) {
            cat("Detrending model: ")
            print(model)
        }
        regressor <- regression.poly
        modstr <- as.character(model)[2:3]
        regress.str <- paste(modstr, collapse = " ~ ")
    }  else if (charmatch(detrend, "bridge", nomatch = FALSE)) {
        regressor <- regression.bridge
        regress.str <- "Bridge detrended"
    } else if (charmatch(detrend, "none", nomatch = FALSE)) {
        regressor <- regression.none
        regress.str <- "None"
    } else stop("Detrending method is not supported")
    sum.order <- trunc(sum.order)
    if (sum.order > 0) {
        for (i in seq(sum.order)) x <- cumsum(x)
    }   else if (sum.order < 0) {
        for (i in seq(sum.order)) x <- diff(x)
    }
    N <- length(x)
        if (is.null(scale.min)) 
        scale.min <- ifelse(substring(detrend, 1, 4) == "poly", 
            2 * (polyfit.order + 1), min(N/4, 4))
    checkScalarType(scale.min, "numeric")
    checkScalarType(scale.max, "numeric")
    scale.min <- trunc(scale.min)
    scale.max <- trunc(scale.max)
    if (scale.min > scale.max) 
        stop("Scale minimum cannot exceed scale maximum")
    if (any(c(scale.min, scale.max) < 0)) 
        stop("Scale minimum and maximum must be positive")
    if (any(c(scale.min, scale.max) > N)) 
        stop(paste("Scale minimum and maximum must not exceed length", 
            "of (differenced/cummulatively summed) time series"))
    scale <- logScale(scale.min, scale.max, scale.ratio = scale.ratio, 
        coerce = trunc)
    scale <- scale[scale > 1]
    
    #Flips the ordering of the scale so very small bins are done first. 
    scale[order(scale)]->scale
    

   #This function creates the starting positions for bin size (scale) sc.  
     scale.distrib<-function(sc,x=x,overlap=overlap){
    	#print("Start scale.distrib")
    	#print(sc)
    	#overlap<-.3
        noverlap <- trunc(overlap * sc)
        noverlap
        sc-noverlap->movesize
        movesize
        trunc(N/movesize)->nmoves
        nmoves
        N
        N-(movesize*nmoves)
        (((nmoves-1)*movesize)+sc)
        (N-(((nmoves-1)*movesize)+sc) > 0)
        
    	
    	#sets the starting positions for two passess through the series. 
    	#First pass starts from the first value, second starting point determined by the end value if the number of moves does not go evenly into the time series
    	#Others starts at a small random value
    	if(N-(((nmoves-1)*movesize)+sc) > 0){
	     start2<-N-(((nmoves-1)*movesize)+sc)
	    } else {
         start2<-sample(10,1)
         #print(start2)
         (N-start2)/movesize->nmoves
	    }
    	seq(from=1,by=movesize,length=trunc(N/movesize))->start.list
        seq(from=start2,by=movesize,length=trunc((N-start2)/movesize))->start.list2
        c(start.list,start.list2)->start.all
        start.all
        data.frame(start.all,sc)->out
        return(out)
    } 

  
   #Calls scale distrib for all scales.   
   lapply(scale,scale.distrib,x,overlap)->out2
    do.call("rbind",out2)->out3
#    summary(out3)
#    print(out3)
    #Does same as above, requires reshape library
    #library(reshape)
    #rbind.fill(out2)->out3
    seq(nrow(out3))->rownum
    cbind(out3,rownum)->out4


nrow(out4)->totrow
#framesize<-50

#Compute the number of frames of size "framesize" will be used.
ifelse(((totrow/framesize)!=trunc(totrow/framesize)),(trunc(totrow/framesize)+1),(totrow/framesize))->nframes

#Create a list of rows with the proper starting positions 
#object out4 has the list of binsize and starting positions
seq(from=1,to=totrow,by=framesize)->frame.start.list

#frame.start.list
#totrow
#frame.start.list[1]
#
#frame.process(frame.start.list[1],framesize,totrow,x)->fr.out

#This function computes the regression using the model and method selected by the function call
 binresid<-function(rownum,x,model){
    	   #print("Start binresid")	
    	  # print("rownum")
#    	   print(rownum)           
    	   start<-out4[rownum,]$start.all
    	   sc<-out4[rownum,]$sc
    	       	   
    	   #print("start")
#    	   print(start)
#    	   print(sc)
	       end<-start+sc
	       #print("end")
#	       print(end)
	       t<-seq(length(x[start:end]))
	       index<-seq(start,end)
	#       x <- create.signalSeries(x)
#           x <- x@data
#           fit<-lm(x[index]~1+t)
#           resid<-(sum(fit$residuals^2))
           resid<-(regressor(x[index], model = model))
           #For the Hulber HRL videos I dumped the catch below. Now it is handled before submission to the qfunc_full function. 
           # The code to catch the zeros and replace with next smallest value is:
           #    min(flucts.out[flucts.out$resid>0,]$resid)->minresid
           #   ifelse(flucts.out$resid==0,minresid,flucts.out$resid)->flucts.out$resid
           
#           if (resid==0){
#        		backresid<-0
#        		forwardresid<-0
#        		backindex <- index - (sc - noverlap)
#        		forwardindex<-index + sc - noverlap
#        		if (count>1) {
#        	      backresid<-(regressor(x[backindex], model = model))
#        		}
#        		if (forwardindex[sc] <= N){
#        	      forwardresid<-(regressor(x[forwardindex], model = model))
#        		} 
#        		resid<-max(backresid,resid,forwardresid)
#        		}
           data.frame(resid,start,sc)->tmp
           return(tmp)
   }
    

#This function get passed to each core (using mclapply) and processes some number (set by framesize) of regressions. 

frame.process<-function(row.start,framesize,totrow,x){
	#print("row.start")
#	print(row.start)
    ifelse(row.start+framesize>totrow,totrow,(row.start+framesize-1))->row.end 
 #   print("row.end")
 #   print(row.end)
	for (i in row.start:row.end) {
		#print("i")
#		print(i)
		binresid(i,x,model)->tmp
		#print(tmp)
		ifelse(i==row.start,tmp->frame.out,rbind(frame.out,tmp)->frame.out)
		}
	return(frame.out) 
	}


mclapply(frame.start.list,frame.process,framesize,totrow,x)->fr.out
do.call("rbind",fr.out)->fr.out
#summary(fr.out)
 return(fr.out)
 }





    
    
    
    
    