iter <- function(x, y, delta, beta)
{
# This function computes one iteration of the EM for 
# censored regression est. (Buckley-James) It first order the
# data according to the residuals, then 
# compute conditional expectations, then compute a new beta by lm(). 
#
# Input: 
# x is a matrix of N rows (the covariates), 
# y is the censored response, a vector of length N.
# delta is a vector of length N.
# (delta =1 means (y) is uncensored. delta = 0 means
# the  (y) is censored.)
# beta is the initial est. and is a vector of length = no. of column(x)
#
# Output: the new value of beta, 
#
        N <- length(delta)
        u <- x %*% beta
        res <- y - u
        niceorder <- order(res, - delta) # order the obs according to
        resorder <- res[niceorder]     # res, if tie then according to
        dorder <- delta[niceorder]  # delta value i.e. d=1 comes first
        uorder <- u[niceorder] 
        ystar <- y[niceorder]  # should I just let ystar <- delta ?
        xorder <- as.matrix(x[niceorder,])

temp <- WKM(x=resorder, d=dorder)

jifen <- rev( cumsum( rev(resorder * temp$jump)) )
Sresorder <- temp$surv

for (i in 1:N) if (dorder[i] == 0) {
           ystar[i] <- uorder[i] + jifen[i]/Sresorder[i]
}

return( lm( ystar ~ 0 + xorder )$coef )
}
