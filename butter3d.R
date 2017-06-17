butter.3d<-function(bf,df){
   #The filter creates some transients at both ends, so trimming off 50 points
   x <- filtfilt(bf, df[,1])
   y <- filtfilt(bf, df[,2])
   z <- filtfilt(bf, df[,3])
   #filt.eu<-data.frame(b.X,b.Y,b.Z)
   # eu[51:(nrow(eu)-50),]$Marker->Marker
   out<-data.frame(x,y,z)
   return(out)
}
