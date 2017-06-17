downsample.3d <- function(ds.rate,df){
   x <- decimate(df[,1],q = ds.rate)
   y <- decimate(df[,2],q = ds.rate)
   z <- decimate(df[,3],q = ds.rate)
   out<-data.frame(x,y,z)
   return(out)
}
