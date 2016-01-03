#include "collision.h"
#include "LBDefinitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq)
{
    for (int i = 0; i < PARAMQ; i++)
    {
        currentCell[i] = currentCell[i] * (1 - 1./ (*tau )) + 1./ (*tau)  * feq[i];
    }
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength)
{
    double velocity[3], density,  feq[PARAMQ], *currentCell;
    for (int z = 1; z <= xlength; z++)
        for (int y = 1; y <= xlength ; y++)
            for (int x = 1; x <= xlength; x++)
            {
                currentCell = collideField + PARAMQ * (z * (xlength + 2) * (xlength + 2) + y * (xlength + 2) + x);
                density = 0.0;
                computeDensitySSE(currentCell, &density);

                if (__builtin_expect(fabs(1.0 - density) > 0.1, 0))
                    fprintf(stderr, "WARNING: Density is %.3f in cell (%d, %d, %d)\n", density, x, y, z);

                computeVelocitySSE(currentCell, &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);

            }
}

