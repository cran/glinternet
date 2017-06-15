#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP R_compute_norms_cat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_compute_norms_cat_cat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_compute_norms_cat_cont(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_compute_norms_cont_cont(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_gl_solver(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_initialize_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_rescale_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_retrieve_beta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_x_times_rescaled_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_compute_norms_cat",       (DL_FUNC) &R_compute_norms_cat,        7},
    {"R_compute_norms_cat_cat",   (DL_FUNC) &R_compute_norms_cat_cat,    9},
    {"R_compute_norms_cat_cont",  (DL_FUNC) &R_compute_norms_cat_cont,  11},
    {"R_compute_norms_cont_cont", (DL_FUNC) &R_compute_norms_cont_cont,  9},
    {"R_gl_solver",               (DL_FUNC) &R_gl_solver,               18},
    {"R_initialize_beta",         (DL_FUNC) &R_initialize_beta,         16},
    {"R_rescale_beta",            (DL_FUNC) &R_rescale_beta,            13},
    {"R_retrieve_beta",           (DL_FUNC) &R_retrieve_beta,            5},
    {"R_x_times_rescaled_beta",   (DL_FUNC) &R_x_times_rescaled_beta,   12},
    {NULL, NULL, 0}
};

void R_init_glinternet(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
