#include <math.h>
#include <R.h>

int sign_d(double x);

void lineAboveNoise(int * beginning,int * posLeft,int * posRight,double * intval,int * sizeMz,int * size_min,double * noise);

void findingPeakLimit(int * posLeft, int * posRight,double * intval,int *sizeMz);

void findLimDensity(double * seq,int * size,int * istart,int * linflex, int * rinflex,int * state);

void FindEqualGreaterM(const double *in, const int *size, const double *values,
                       const int *valsize, int *index);

void findLim(double *yseq,int *istartl,int *istartr,int *lengthmax, int *ileft, int *iright);

void findMax(double *yseq,int *istart,int *res);

void binarySearch(double *xseq, double *val, int *imin, int *imax, int *imid);

void replaceByZero(double *rseq,double *yseq,int* sizex);

void linearInterpolation(double *intseq,int *imin, int *imax, int *counter);

void checkIso(double *intensity, int *num_iso, int *valid, int *smax);
