#ifndef __PYX_HAVE__dmaps__critical__bias
#define __PYX_HAVE__dmaps__critical__bias

struct BiasedMD;
struct DMSConfig;
struct Fit;
struct FEHist;

/* "dmaps/critical/bias.pxd":3
 * cimport numpy as np
 * 
 * cdef public struct BiasedMD:             # <<<<<<<<<<<<<<
 *     int natoms # number of atoms
 *     int nheavy_atoms # number of heavy atoms
 */
struct BiasedMD {
  int natoms;
  int nheavy_atoms;
  int step;
  int nsteps;
  int *heavy_atoms_idxs;
  double vbias;
  double *dcs;
  float *coord;
  float *force;
};

/* "dmaps/critical/bias.pxd":14
 *     float* force # force
 * 
 * cdef public struct DMSConfig:             # <<<<<<<<<<<<<<
 *     int isfirst # 1 if first dmaps iter, 0 if not
 *     int nstride # number of configs to save
 */
struct DMSConfig {
  int isfirst;
  int nstride;
  int ndcs;
  double fefrac;
  double kT;
};

/* "dmaps/critical/bias.pxd":21
 *     double kT # kT value
 * 
 * cdef public struct Fit:             # <<<<<<<<<<<<<<
 *     int npoints # number of points used to fit
 *     char* function # function name
 */
struct Fit {
  int npoints;
  char *function;
  char *metric;
  double *weights;
  double *sigma;
  double *coords;
};

/* "dmaps/critical/bias.pxd":29
 *     double* coords # coordinates of points used to fit
 * 
 * cdef public struct FEHist:             # <<<<<<<<<<<<<<
 *     int nbins # number of bins along each dimension
 *     int nnebins # number of non-empty bins (overall)
 */
struct FEHist {
  int nbins;
  int nnebins;
  int *nebins_idxs;
  int *nebins_idxs_s;
  double *steps;
  double *bins;
  double *values;
  double *gradient;
};

#ifndef __PYX_HAVE_API__dmaps__critical__bias

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(struct DMSConfig) *initDMSConfig(char const *);
__PYX_EXTERN_C DL_IMPORT(struct FEHist) *initFEHist(struct DMSConfig *, char const *);
__PYX_EXTERN_C DL_IMPORT(struct Fit) *initFit(struct DMSConfig *, char const *);
__PYX_EXTERN_C DL_IMPORT(struct BiasedMD) *initBiasedMD(struct DMSConfig *, char const *);
__PYX_EXTERN_C DL_IMPORT(int) do_biased_force(struct BiasedMD *, struct DMSConfig *, struct Fit *, struct FEHist *);

#endif /* !__PYX_HAVE_API__dmaps__critical__bias */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initbias(void);
#else
PyMODINIT_FUNC PyInit_bias(void);
#endif

#endif /* !__PYX_HAVE__dmaps__critical__bias */
