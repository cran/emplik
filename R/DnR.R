DnR <- function(x, d, w)
{
# inputs should be from  Wdataclean2()

allrisk <- rev(cumsum(rev(w)))
posi <- d == 1
uncenx <- x[posi]
uncenw <- w[posi]
uncenR <- allrisk[posi]

list( time = uncenx, n.risk = uncenR, n.event = uncenw )
}

