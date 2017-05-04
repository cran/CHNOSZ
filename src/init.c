#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(h2o92)(void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"h2o92", (DL_FUNC) &F77_NAME(h2o92), 4},
    {NULL, NULL, 0}
};

void R_init_CHNOSZ(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
