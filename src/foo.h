
#ifndef EMPLIK_FOO_H
#define EMPLIK_FOO_H

#include <R.h>
#include <Rinternals.h>

void cumsumsurv(double *x, double *s, int *nin);
void eltestwt(double *x, double *w, double *lam, int *nin, double *mu, double *del);


#endif /* EMPLIK_FOO_H */
