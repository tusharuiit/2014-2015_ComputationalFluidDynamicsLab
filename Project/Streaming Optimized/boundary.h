#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, int* flagField, const double * const bddParams,int *xlength);

/** compute boundary helper function */
void compute_boundary(double *collideField, const double * const bddParams, int *flagField, int *xlength, const int directionToUpdate, int * const coordinate, int * const freeSlipNeighbours, const int freeSlipDirectionToCopy, const double * const inflowFeq);
#endif

