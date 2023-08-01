#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _PREGS_permutation_conformal_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PREGS_permutation_FL_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PREGS_permutation_simple_C(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_PREGS_permutation_conformal_C", (DL_FUNC) &_PREGS_permutation_conformal_C, 6},
  {"_PREGS_permutation_FL_C",        (DL_FUNC) &_PREGS_permutation_FL_C,        6},
  {"_PREGS_permutation_simple_C",    (DL_FUNC) &_PREGS_permutation_simple_C,    4},
  {NULL, NULL, 0}
};

void R_init_PREGS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}