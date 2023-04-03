
#ifndef EMPLIK_FOO_H
#define EMPLIK_FOO_H

#include <R.h>
#include <Rinternals.h>

void cumsumsurv(const double *x, double *s, const int *LLL);
void eltestwt(double *x, double *wt, const double *mu1, const int *Lx1, double *prob, double *lamre);

#endif /* EMPLIK_FOO_H */
