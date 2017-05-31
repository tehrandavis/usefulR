#This is the command to submit this file to R in the terminal window : "R --vanilla <Particle_MFDFA_ManTrack_Terminal.R"
#The terminal and R.app are both pointed to the same directory. R.app is the "R" application that runs in the Mac standard interface.
#The functions called in the program reside in Dropbox/R_Stable as saved R objects that are loaded with the commands below.
#For this example, the terminal command to set the directory would be: cd /Users/jdixon_old/Documents/Hulber_HRL_Videos/Quads_Video1_HRL
# Make the edits to this file (saved under a new name) in R.app. Then submit in the terminal with the command above.



library(gtools)
library(gdata)
#library(gregmisc)
library(Hmisc)
library(lme4)
#library(plotrix)
library(lattice)
#library(latticeExtra)
library(reshape)
library(fractal)
library(multicore)

#Load the functions as R objects

load("/Users/jdixon_old/Dropbox/R_Stable/flucts.mc.distrib.RData")
load("/Users/jdixon_old/Dropbox/R_Stable/qfunc_full.RData")


#Sets the working directory
#Terminal will point to this directory too
setwd("/Users/jdixon_old/Documents/Hulber_HRL_Videos/Quads_Video1_HRL")
getwd() 

#This section is just makes a list of the files with that have the data. They all start with "Points_MTrackJ_".
dir()->file.names
file.names[]
length(file.names)
grep("Points_MTrackJ_",file.names)->listpoints
file.names[listpoints]->file.points
file.points

for (i in 1:length(file.points)){
read.table(file.points[i],header=FALSE,sep="\t",skip=1) -> x
summary(x)
x[1:400,]
names(x)
names(x)<-c("Frame","TID","PID","x","y","t","intensity","len","D2S","D2R","D2P","v")
cbind(x$x,x$y)->dist
summary(dist)
#print(summary(eu))

#This section computes the euclidean distance between temporally successive points
 data.matrix(dist, rownames.force=FALSE) ->y
 nrow(y)->rows.y
 yl1<-y[2:rows.y,2]
 yl0<-y[1:rows.y-1,2]
 xl0<-y[1:rows.y-1,1]
 xl1<-y[2:rows.y,1]
 yd2<-(yl0-yl1)^2
 xd2<-(xl0-xl1)^2
 eu<-sqrt((xd2+yd2))
#plot(eu)

#Begin Epoch analysis; epoching is optional, can be turned off by setting start to 1 and end to the length of eu :end <-length(eu); and 
# commenting out the start and end increments in the while loop
array(data=NA,c(1,1))->tmp.array
as.data.frame(tmp.array)->tmp
 epoch<-1
 start<-1
 end<-start+600         
 summary(eu[start:end])
 str(eu)
 while (end < length(eu)) {
 	#The flucts.mc.distrib function has quite a few options
 	  #sum.order controls integration or differentiation of the time series. 1 is first-order integration, -1 is first-order differentiation
 	  #overlap sets the proportion of overlap between adjacent bins 
 	  #scale.ratio controls the number of bins
 	  #framesize controls how many bins are distributed to each core (only affects processing speed)
    flucts.mc.distrib(eu[start:end],sum.order=1,overlap=.3,scale.ratio=.9,framesize=100)->flucts.out
    save(flucts.out,file="flucts.out.RData")
    print(summary(flucts.out))
    #Trap to replace zeros with next smallest value.
    #Zero residuals cannot be log transformed.
    min(flucts.out[flucts.out$resid>0,]$resid)->minresid
    ifelse(flucts.out$resid==0,minresid,flucts.out$resid)->flucts.out$resid
    #The qfunc_full function also has some options. When running in terminal leave onplot off (i.e., onplot=0). 
    #The qlo and qhi set the min and max values for q
    qfunc_full(flucts.out,onplot=0,qlo=-4,qhi=8)->DFA.eu
    #DFA.eu
	#DFA(eu[start:end],detrend="poly1",sum.order=1, overlap=.3,verbose=FALSE) ->DFA.eu
	#if (count < 30){
	#	plot(DFA.eu)}
	#count<-count+1	
    DFA.eu$epoch<-epoch
    DFA.eu$start<-start
    DFA.eu$end<-end
    DFA.eu$file<-i
    ifelse (start==1 & i==1, DFA.eu->out.DFA, rbind(out.DFA,DFA.eu)->out.DFA)
    start<-start+100
    end<-start+600
    print(epoch)	
    epoch<-epoch+1
 }
} 
 save(out.DFA,file="out.DFA.RData")
