removeParetheses <- function(df){
  # how many columns?
  lastColumn <- ncol(df)
  gsub(pattern="\\(", replacement="", df[,1]) -> df[,1]
  gsub(pattern="\\)", replacement="", df[,lastColumn]) -> df[,lastColumn]
  df[,1] <- as.numeric(df[,1])
  df[,lastColumn] <-as.numeric(df[,lastColumn])
  return(df)
  }