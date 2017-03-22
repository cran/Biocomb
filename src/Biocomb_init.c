#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .Call calls */
  extern SEXP Biocomb_CalcGene(SEXP, SEXP, SEXP, SEXP);
extern SEXP Biocomb_CalcROC(SEXP, SEXP, SEXP);
extern SEXP Biocomb_check_incons(SEXP, SEXP, SEXP);
extern SEXP Biocomb_forward_path(SEXP, SEXP);
extern SEXP Biocomb_fun1_chi(SEXP, SEXP);
extern SEXP Biocomb_fun2_chi(SEXP, SEXP);
extern SEXP Biocomb_fun3_chi(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Biocomb_fun4_chi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"Biocomb_CalcGene",     (DL_FUNC) &Biocomb_CalcGene,      4},
  {"Biocomb_CalcROC",      (DL_FUNC) &Biocomb_CalcROC,       3},
  {"Biocomb_check_incons", (DL_FUNC) &Biocomb_check_incons,  3},
  {"Biocomb_forward_path", (DL_FUNC) &Biocomb_forward_path,  2},
  {"Biocomb_fun1_chi",     (DL_FUNC) &Biocomb_fun1_chi,      2},
  {"Biocomb_fun2_chi",     (DL_FUNC) &Biocomb_fun2_chi,      2},
  {"Biocomb_fun3_chi",     (DL_FUNC) &Biocomb_fun3_chi,      5},
  {"Biocomb_fun4_chi",     (DL_FUNC) &Biocomb_fun4_chi,     13},
  {NULL, NULL, 0}
};

void R_init_Biocomb(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
