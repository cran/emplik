##################################################
### 2-sample right-censored data. Poisson EL   ###
### parameter = difference of 2 mean Haz       ###
##################################################
emplikH2P <- function(lambda, x1, d1, x2, d2, fun1, fun2, CIforTheta=FALSE) {
    x1 <- as.vector(x1)
    n1 <- length(x1)
    if (n1 <= 2)
        stop("Need more observations in x1")
    if (length(d1) != n1)
        stop("length of x1 and d1 must agree")
    if (any((d1 != 0) & (d1 != 1)))
        stop("d1 must be 0/1's for censor/not-censor")
    if (!is.numeric(x1))
        stop("x1 must be numeric -- observed times")
    if(!is.numeric(lambda)) stop("lambda must be a numeric value -- tilt")

    x2 <- as.vector(x2)
    n2 <- length(x2)
    if (n2 <= 2)
        stop("Need more observations for sample 2")
    if (length(d2) != n2)
        stop("length of x2 and d2 must agree")
    if (any((d2 != 0) & (d2 != 1)))
        stop("d2 must be 0/1's for censor/not-censor")
    if (!is.numeric(x2))
        stop("x2 must be numeric -- observed times")

    newdata1 <- Wdataclean2(z=x1, d=d1)
    temp1 <- DnR(newdata1$value, newdata1$dd, newdata1$weight)
    newdata2 <- Wdataclean2(z=x2, d=d2)
    temp2 <- DnR(newdata2$value, newdata2$dd, newdata2$weight)

    jump1 <- (temp1$n.event)/temp1$n.risk
    jump2 <- (temp2$n.event)/temp2$n.risk

    funtime11 <- fun1(temp1$times)
    funtime21 <- fun2(temp2$times)

    index1 <- (jump1 < 1)
    index2 <- (jump2 < 1)

    K12 <- 0
    tm11 <- temp1$times[!index1]
    if(length(tm11) > 1 ) stop("more than 1 place jump>=1 in x1?")
    if( length(tm11) > 0 ) {
         K12 <- K12 + fun1(tm11)
    }
    tm21 <- temp2$times[!index2]
    if(length(tm21) > 1 ) stop("more than 1 place jump>=1 in x2?")
    if( length(tm21) > 0 ) {
         K12 <- K12 - fun2(tm21)
    }

    eve1 <- temp1$n.event[index1]
    tm1 <- temp1$times[index1]
    rsk1 <- temp1$n.risk[index1]
    jmp1 <- jump1[index1]
    funtime1 <- fun1(tm1)
#######################################################
#### it seems I need to include the last point, even it is 1???
########Sample two########
    eve2 <- temp2$n.event[index2]
    tm2 <- temp2$times[index2]
    rsk2 <- temp2$n.risk[index2]
    jmp2 <- jump2[index2]
    funtime2 <- fun2(tm2)
##############################################################
  lam <- lambda
  N <- max( c(2*(n1+n2), 10000) )  ### for the llog function ##
###########################################################
lamfun1 <- funtime1 * lam
lamfun2 <- funtime2 * lam
onePlamh1 <- (rsk1 + lamfun1)/rsk1   ### this is 1 + lam Zi in Ref.
oneMlamh2 <- (rsk2 - lamfun2)/rsk2   ### this is 1 - lam Zi in Ref.

DHazPara <- NA
if(CIforTheta) {
DHazPara <- K12 + sum(funtime1*jmp1/onePlamh1) - sum(funtime2*jmp2/oneMlamh2)
}
###weights <- jump/onepluslamh
###need to change last jump to 1? NO. see above

loglik1 <- (sum( eve1*llog(onePlamh1, 1/N)) -
                 sum(eve1*(lamfun1)/(rsk1 + lamfun1)) )

loglik2 <- (sum( eve2*llog(oneMlamh2, 1/N)) -
                 sum(eve2*(-lamfun2)/(rsk2 - lamfun2)) )
loglik <- 2*(loglik1 + loglik2)
#?is that right? YES  see (3.2) in Ref. above. This ALR, or Poisson LR.
#last jump of 1 can be excluded when computing the EL ratio
list("-2LLR"=loglik, lambda=lam, "-2LLR(sample1)"=2*loglik1, HazDiff=DHazPara)
}