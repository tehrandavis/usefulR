DFA.cm <- function (x, detrend = "poly1", sum.order = 0, overlap = 0, scale.max = trunc(length(x)/2), 
          scale.min = NULL, scale.ratio = 2, verbose = FALSE) 
{
  "polyfit.model" <- function(polyfit.order) {
    if (polyfit.order < 0) 
      stop("Polynomial fit order must be positive")
    if (polyfit.order == 0) 
      return("x ~ 1")
    if (polyfit.order == 1) 
      return("x ~ 1 + t")
    else {
      return(paste(c("x ~ 1", "t", paste("I(t^", seq(2, 
                                                     polyfit.order), ")", collapse = " + ", sep = "")), 
                   collapse = " + ", sep = ""))
    }
  }
  "regression.poly" <- function(x, model = polyfit.model(1)) {
    order <- length(attr(terms(formula(x ~ 1)), "order"))
    t <- x@positions
    x <- x@data
    if (order == 0) 
      return(sum((x - mean(x))^2))
    if (order == 1) {
      x <- x - mean(x)
      t <- t - mean(t)
      slope <- sum(t * x)/sum(t^2)
      return(sum((x - slope * t)^2))
    }
    fit <- lm(model, data = data.frame(list(t = t, x = x)))
    return(sum(fit$residuals^2))
  }
  "regression.bridge" <- function(x, ...) {
    x <- x@data
    N <- length(x)
    bridge <- seq(x[1], x[N], length = N)
    x <- x - bridge
    return(sum(x^2))
  }
  "regression.none" <- function(x, ...) return(sum(x@data^2))
  checkScalarType(detrend, "character")
  checkScalarType(scale.ratio, "numeric")
  checkScalarType(verbose, "logical")
  checkScalarType(overlap, "numeric")
  data.name <- deparse(substitute(x))
  x <- wmtsa::create.signalSeries(x)
  if ((overlap < 0) | (overlap >= 1)) 
    stop("Overlap factor must be in the range [0,1)")
  detrend <- lowerCase(detrend)
  if (substring(detrend, 1, 4) == "poly") {
    polyfit.order <- as.integer(substring(detrend, 5))
    if (is.na(polyfit.order)) {
      stop("Improperly formed detrending string")
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
  }
  else if (charmatch(detrend, "bridge", nomatch = FALSE)) {
    regressor <- regression.bridge
    regress.str <- "Bridge detrended"
  }
  else if (charmatch(detrend, "none", nomatch = FALSE)) {
    regressor <- regression.none
    regress.str <- "None"
  }
  else stop("Detrending method is not supported")
  sum.order <- trunc(sum.order)
  if (sum.order > 0) {
    for (i in seq(sum.order)) x <- cumsum(x)
  }
  else if (sum.order < 0) {
    for (i in seq(sum.order)) x <- diff(x)
  }
  N <- length(x)
  if (is.null(scale.min)) 
    scale.min <- ifelse1(substring(detrend, 1, 4) == "poly", 
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
  rmse <- rep(0, length(scale))
  for (i in seq(along = scale)) {
    sc <- scale[i]
    noverlap <- trunc(overlap * sc)
    if (verbose) 
      cat("Processing scale", sc, "\n")
    cumsum <- 0
    count <- 0
    index <- seq(sc)
    while (index[sc] <= N) {
      rmse[i] <- rmse[i] + regressor(x[index], model = model)
      count <- count + 1
      index <- index + sc - noverlap
    }
    rmse[i] <- sqrt(rmse[i]/(count * sc))
  }
  xx <- log(scale)
  yy <- log(rmse)
  ww <- 1/seq(along = xx)
  w <- NULL
  logfit <- lm(y ~ 1 + x, data = data.frame(x = xx, y = yy, 
                                            w = ww), weights = w)
  exponent <- logfit$coefficients["x"]
  list(exponent,rmse)
  
  
}
