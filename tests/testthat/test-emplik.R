UNITEST_el.cen.EM <- function(){
    x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
    d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)
    result = el.cen.EM(x,d,mu=3.5)
    return(sprintf('%0.5f', result[["lam"]][1]));

}

UNITEST_findUL <- function(){
     x <- c(1, 1.5, 2, 3, 4, 5, 6, 5, 4, 1, 2, 4.5)
     d <- c(1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,  1)
     myfun6 <- function(theta, x, d) {
     el.cen.EM2(x, d, fun=function(t){t}, mu=theta)
     }
    result = findUL(step=0.2, fun=myfun6, MLE=4.0659, x=x, d=d)
    return(sprintf('%0.5f', result[["Low"]][1]));
}

UNITEST_el.cen.kmc1d <- function(){
     x <- c(1, 1.5, 2, 3, 4.2, 5, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3)
     d <- c(1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1)
     ff <- function(x) {
         x - 3.7
     }
     
    result = el.cen.kmc1d(x=x, d=d, fun=ff, mu=0)
    return(sprintf('%0.5f', result[["-2LLR"]][1]));
}



test_that("emplik works", {
  expect_equal(UNITEST_el.cen.EM(), "2.11514");
  expect_equal(UNITEST_findUL(), "3.04973");
  expect_equal(UNITEST_el.cen.kmc1d(), "0.61496")
  
})
