#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(prxgrd)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(prxgrdf)(void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"prxgrd",  (DL_FUNC) &F77_NAME(prxgrd),  7},
    {"prxgrdf", (DL_FUNC) &F77_NAME(prxgrdf), 7},
    {NULL, NULL, 0}
};

void R_init_covchol(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}