#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, int* flagField, const double * const bddParams,int *xlength);

/** compute boundary helper function */
void compute_boundary(double *collideField, const double * const bddParams, int *flagField, int *xlength, const int const *iList, int *const coordinate, int * const freeSlipNeighbours, const int const * fsList, const double const * inflowFeq);
#endif

