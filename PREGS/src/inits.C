#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _PREGS_permutation_FL(SEXP, SEXP, SEXP, SEXP);
extern SEXP _PREGS_permutation_PREGSjoint(SEXP, SEXP, SEXP, SEXP);
extern SEXP _PREGS_permutation_PREGSseparate(SEXP, SEXP, SEXP, SEXP);
extern SEXP _PREGS_permutation_vanilla(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_PREGS_permutation_FL",            (DL_FUNC) &_PREGS_permutation_FL,            4},
  {"_PREGS_permutation_PREGSjoint",    (DL_FUNC) &_PREGS_permutation_PREGSjoint,    4},
  {"_PREGS_permutation_PREGSseparate", (DL_FUNC) &_PREGS_permutation_PREGSseparate, 4},
  {"_PREGS_permutation_vanilla",       (DL_FUNC) &_PREGS_permutation_vanilla,       4},
  {NULL, NULL, 0}
};

void R_init_PREGS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}