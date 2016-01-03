#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdlib.h>

/**
 *  The array bddParams has the following structure
 *  bddParams[0]   Inflow velocity in x direction
 *  bddParams[1]   Inflow velocity in y direction
 *  bddParams[2]   Inflow velocity in z direction
 *  bddParams[3]   Pressure surplus for PRESSURE_IN cells
 *  bddParams[4]   Moving wall velocity in x direction
 *  bddParams[5]   Moving wall velocity in y direction
 *  bddParams[6]   Moving wall velocity in z direction
 */
void treatBoundary(double *collideField, int* flagField, const double * const bddParams, int *xlength)
{
    int i, x, y, z, neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex;
    int coordinate[3];
    int freeSlipNeighbours[3];
    double referenceDensity = 1.0;
    double inflowFeq[PARAMQ];
    computeFeqAVX(&referenceDensity, bddParams, inflowFeq);

    /* top boundary
     * only need to set directions 0, 5, 6, 7, 14 */
    y = coordinate[1] = xlength[1] + 1;
    freeSlipNeighbours[1] = xlength[1];
    for (z = 0; z <= xlength[2] + 1; z++)
        for (x = 0; x <= xlength[0] + 1; x++)
        {
            coordinate[0] = freeSlipNeighbours[0] = x;
            coordinate[2] = freeSlipNeighbours[2] = z;
            compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESBOTTOMOUT, coordinate, freeSlipNeighbours, VELOCITIESTOPOUT, inflowFeq);
        }


    /* back boundary */
    // i = 14,15,16,17,18
    z = coordinate[2] = 0;
    freeSlipNeighbours[2] = 1;
    for (y = 0; y <= xlength[1] + 1; y++)
        for (x = 0; x <= xlength[0] + 1; x++)
        {
            coordinate[0] = freeSlipNeighbours[0] = x;
            coordinate[1] = freeSlipNeighbours[1] = y;
            compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESFRONTOUT, coordinate, freeSlipNeighbours, VELOCITIESBACKOUT, inflowFeq);

        }


    /* bottom boundary */
    // i = 4, 11, 12, 13, 18
    y = coordinate[1] = 0;
    freeSlipNeighbours[1] = 1;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (x = 0; x <= xlength[0] + 1; x++)
        {
            coordinate[0] = freeSlipNeighbours[0] = x;
            coordinate[2] = freeSlipNeighbours[2] = z;
            compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESTOPOUT, coordinate, freeSlipNeighbours, VELOCITIESBOTTOMOUT, inflowFeq);
        }

    /* left boundary */
    // i = 3, 7 ,10, 13, 17
    coordinate[0] = x = 0;
    freeSlipNeighbours[0] = 1;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++)
        {
            coordinate[1] = freeSlipNeighbours[1] = y;
            coordinate[2] = freeSlipNeighbours[2] = z;
            compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESRIGHTOUT, coordinate, freeSlipNeighbours, VELOCITIESLEFTOUT, inflowFeq);
        }

    /* right boundary */
    // i = 1,5,8,11,15
    coordinate[0] = x = xlength[0] + 1;
    freeSlipNeighbours[0] = xlength[0];
    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++)
        {
            coordinate[1] = freeSlipNeighbours[1] = y;
            coordinate[2] = freeSlipNeighbours[2] = z;
            compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESLEFTOUT, coordinate, freeSlipNeighbours, VELOCITIESRIGHTOUT, inflowFeq);
        }

    /* front boundary, i.e. z = xlength + 1 */
    // i=0,1,2,3,4
    z = coordinate[2] = xlength[2] + 1;
    freeSlipNeighbours[2] = xlength[2];
    for (y = 0; y <= xlength[1] + 1; y++)
        for (x = 0; x <= xlength[0] + 1; x++)
        {
            coordinate[0] = freeSlipNeighbours[0] = x;
            coordinate[1] = freeSlipNeighbours[1] = y;
            compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESBACKOUT, coordinate, freeSlipNeighbours, VELOCITIESFRONTOUT, inflowFeq);
        }

    /* inner boundary cells
     * assumes inner boundary cells can only be NO_SLIP */
    for (z = 1; z <= xlength[2]; z++)
        for (y = 1; y <= xlength[1]; y++)
            for(x = 1; x <= xlength[0]; x++)
                if (flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
                    for (i = 0; i < PARAMQ; i++)
                    {
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (!flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ])
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
}

/** boundary helper function */
void compute_boundary(double *collideField, const double * const bddParams, int *flagField, int *xlength, const int const *iList, int *const coordinate, int * const freeSlipNeighbours, const int const * fsList, const double const * inflowFeq)
{

    int i;
    int neighbourCoordX, neighbourCoordY, neighbourCoordZ;
    int currentCellIndex, neighbourCellIndex;
    double density;
    const double * const wallVelocity = bddParams + 4;
    double referenceDensity = 1.0;
    double feq[PARAMQ];
    double velocity[3];
    const double pressureIn = bddParams[3];

    int boundaryType = flagField[coordinate[0] + (xlength[0] + 2) * coordinate[1] + (xlength[0] + 2) * (xlength[1] + 2) * coordinate[2]];

    switch (boundaryType)
    {
    case NO_SLIP:
    {
        for(i = 0; i < 5; i++ )
        {
            neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
            neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
            neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
            {
                currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                collideField[currentCellIndex] = collideField[neighbourCellIndex];
            }
        }
        break;
    }

    case MOVING_WALL:
    {
        for(i = 0; i < 5; i++ )
        {
            neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
            neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
            neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
            {
                currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                computeDensityAVX(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                  + 2 * LATTICEWEIGHTS[iList[i]] * (LATTICEVELOCITIES[iList[i]][0] * wallVelocity[0]
                                                          + LATTICEVELOCITIES[iList[i]][1] * wallVelocity[1] + LATTICEVELOCITIES[iList[i]][2] * wallVelocity[2]) * density / (C_S * C_S);
            }
        }
        break;
    }

    case INFLOW:
    {
        for(i = 0; i < 5; i++ )
        {
            neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
            neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
            neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
            {
                currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                collideField[currentCellIndex] = inflowFeq[iList[i]];
            }
        }

        break;
    }

    case OUTFLOW:
    {
        for(i = 0; i < 5; i++ )
        {
            neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
            neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
            neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
            {
                currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                computeDensityAVX(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                computeVelocityAVX(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                computeFeqAVX(&referenceDensity, velocity, feq);
                collideField[currentCellIndex] = feq[PARAMQ - iList[i] - 1] + feq[iList[i]] - collideField[neighbourCellIndex];
            }
        }
        break;
    }

    case PRESSURE_IN:
    {
        for(i = 0; i < 5; i++ )
        {
            neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[iList[i]][0];
            neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[iList[i]][1];
            neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[iList[i]][2];

            if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                    && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
            {
                currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - iList[i]);
                density = 0.0;
                computeDensityAVX(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                computeVelocityAVX(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                density = referenceDensity + pressureIn;
                computeFeqAVX(&density, velocity, feq);
                collideField[currentCellIndex] = feq[PARAMQ - iList[i] - 1] + feq[iList[i]] - collideField[neighbourCellIndex];
            }
        }
        break;
    }
    case FREE_SLIP:
    {
        neighbourCoordX = freeSlipNeighbours[0];
        neighbourCoordY = freeSlipNeighbours[1];
        neighbourCoordZ = freeSlipNeighbours[2];
        if (!flagField[neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX])
        {
            for(i = 0; i < 5; i++)
            {
                currentCellIndex = PARAMQ * (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + iList[i];
                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + fsList[i];
                collideField[currentCellIndex] = collideField[neighbourCellIndex];
            }
        }
    }
    }

}
