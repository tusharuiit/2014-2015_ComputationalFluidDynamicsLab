#include "collision.h"
#include "LBDefinitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq)
{
    int i;
    for (i = 0; i < PARAMQ; i++)
    {
        currentCell[i] = currentCell[i] * (1 - 1./ (*tau )) + 1./ (*tau)  * feq[i];
    }
}

void doCollision(double *collideField, int *flagField,const double * const tau,int *xlength)
{
    double velocity[3], density, feq[PARAMQ], *currentCell;
    int x, y, z;
    for (z = 1; z <= xlength[2]; z++)
        for (y = 1; y <= xlength[1] ; y++)
            for (x = 1; x <= xlength[0]; x++)
            {
                currentCell = collideField + PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x);
                density = 0.0;
                computeDensityAVX(currentCell, &density);
                if (__builtin_expect(fabs(1.0 - density) > 0.1, 0))
                    fprintf(stderr, "WARNING: Density is %.3f in cell (%d, %d, %d)\n", density, x, y, z);
                computeVelocityAVX(currentCell, &density, velocity);
                computeFeqAVX(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);

            }
}

