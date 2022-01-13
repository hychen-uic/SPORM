#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>   // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h> // for F77_name

// .Fortran calls
extern void F77_NAME(analysiswmg)(void *, void *, void *, void *, void *, void *, void *, void *,
                                  void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(analysiswmgfix)(void *, void *, void *, void *, void *, void *, void *, void *,
                                     void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(networkselectbypm)(void *, void *, void *, void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pwmlecol)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pwmlerow)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pwmlecolfix)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pwmlerowfix)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(networkselectbypw)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormlecol)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormlerow)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormlecolfix)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormlerowfix)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(networkselectbysp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"analysiswmg", (DL_FUNC) &F77_NAME(analysiswmg), 17},
  {"analysiswmgfix", (DL_FUNC) &F77_NAME(analysiswmgfix), 19},
  {"networkselectbypm", (DL_FUNC) &F77_NAME(networkselectbypm), 14},
  {"pwmlecol", (DL_FUNC) &F77_NAME(pwmlecol), 10},
  {"pwmlerow", (DL_FUNC) &F77_NAME(pwmlerow), 10},
  {"pwmlecolfix", (DL_FUNC) &F77_NAME(pwmlecolfix), 12},
  {"pwmlerowfix", (DL_FUNC) &F77_NAME(pwmlerowfix), 12},
  {"networkselectbypw", (DL_FUNC) &F77_NAME(networkselectbypw), 10},
  {"ormlecol", (DL_FUNC) &F77_NAME(ormlecol), 10},
  {"ormlerow", (DL_FUNC) &F77_NAME(ormlerow), 10},
  {"ormlecolfix", (DL_FUNC) &F77_NAME(ormlecolfix), 12},
  {"ormlerowfix", (DL_FUNC) &F77_NAME(ormlerowfix), 12},
  {"networkselectbysp", (DL_FUNC) &F77_NAME(networkselectbysp), 10},
  {NULL, NULL, 0}
};

void R_init_Sporm(DllInfo *dll) { // R looks R_init_<pkg>() to register routines
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
