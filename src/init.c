

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "foo.h"

static R_NativePrimitiveArgType cumsumsurv_types[3] = {REALSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType eltestwt_types[6] = {REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP};

static R_CMethodDef cMethods[] = {
    {"cumsumsurv", (DL_FUNC) &cumsumsurv, 3, cumsumsurv_types},
	{"eltestwt", (DL_FUNC) &eltestwt, 6, eltestwt_types},
    {NULL, NULL, 0, NULL}
};

void attribute_visible R_init_fooR(DllInfo *info)    /* change fooR to emplik? */
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}


