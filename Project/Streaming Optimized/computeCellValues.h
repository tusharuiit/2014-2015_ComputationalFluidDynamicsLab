#ifndef _COMPUTECELLVALUES_H_
#define _COMPUTECELLVALUES_H_
#include "LBDefinitions.h"
/** computes the density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void computeDensity(const double *const currentCell, double *density, const int tot_cells);
/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(const double *const currentCell, const double * const density, double *velocity, const int tot_cells);
/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double * const density, const double * const velocity, double *feq);
void computeFeqAVX(const double * const density, const double * const velocity, double *feq);
void computeFeqAVXv2(double density[4], double velocity[3][4], double feq[PARAMQ][4]);
#endif

