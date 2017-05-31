
 prepanel.ci <- function(x, y, ly, uy, subscripts, ...)
{
    y <- as.numeric(y)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    list(ylim = range(y, uy, ly, finite = TRUE))
}

panel.ci <- function(x, y, ly, uy, subscripts, pch = 16, ...)
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    panel.arrows(x, ly, x, uy, col = c('blue'),
                 length = .25, unit = "native",
                 angle = 90, code =7)
    panel.xyplot(x, y, pch = pch, ...)
}

