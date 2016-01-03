#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

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
    int i, x, y, z;
    #ifdef _ARBITRARYGEOMETRY_
    int neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex;
    #endif // _ARBITRARYGEOMETRY_
    int coordinate[3];
    int freeSlipNeighbours[3];
    double referenceDensity = 1.0;
    double inflowFeq[PARAMQ] __attribute__((aligned(32)));
    computeFeq(&referenceDensity, bddParams, inflowFeq);

    #ifdef _ARBITRARYGEOMETRY_
    #pragma omp parallel private(neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex, x,y,z,i,coordinate,freeSlipNeighbours), firstprivate(referenceDensity), shared(collideField, flagField, xlength, inflowFeq)
    #else
    #pragma omp parallel private(x,y,z,i,coordinate,freeSlipNeighbours), firstprivate(referenceDensity), shared(collideField, flagField, xlength, inflowFeq)
    #endif // _ARBITRARYGEOMETRY_
    {
    y = coordinate[1] = xlength[1] + 1;
    freeSlipNeighbours[1] = xlength[1];
    #pragma omp for nowait schedule(static)
    for (i = 0; i < 5; i++)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
            {
                coordinate[0] = freeSlipNeighbours[0] = x;
                coordinate[2] = freeSlipNeighbours[2] = z;
                compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESBOTTOMOUT[i], coordinate, freeSlipNeighbours, VELOCITIESTOPOUT[i], inflowFeq);
            }

    /* back boundary */
    // i = 14,15,16,17,18
    z = coordinate[2] = 0;
    freeSlipNeighbours[2] = 1;
    #pragma omp for nowait schedule(static)
    for (i = 0; i < 5; i++)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
            {
                coordinate[0] = freeSlipNeighbours[0] = x;
                coordinate[1] = freeSlipNeighbours[1] = y;
                compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESFRONTOUT[i], coordinate, freeSlipNeighbours, VELOCITIESBACKOUT[i], inflowFeq);

            }


    /* bottom boundary */
    // i = 4, 11, 12, 13, 18
    y = coordinate[1] = 0;
    freeSlipNeighbours[1] = 1;
    #pragma omp for nowait schedule(static)
    for (i = 0; i < 5; i++)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
            {
                coordinate[0] = freeSlipNeighbours[0] = x;
                coordinate[2] = freeSlipNeighbours[2] = z;
                compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESTOPOUT[i], coordinate, freeSlipNeighbours, VELOCITIESBOTTOMOUT[i], inflowFeq);
            }

    /* left boundary */
    // i = 3, 7 ,10, 13, 17
    coordinate[0] = x = 0;
    freeSlipNeighbours[0] = 1;
    #pragma omp for nowait schedule(static)
    for (i = 0; i < 5; i++)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
            {
                coordinate[1] = freeSlipNeighbours[1] = y;
                coordinate[2] = freeSlipNeighbours[2] = z;
                compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESRIGHTOUT[i], coordinate, freeSlipNeighbours, VELOCITIESLEFTOUT[i], inflowFeq);
            }

    /* right boundary */
    // i = 1,5,8,11,15
    coordinate[0] = x = xlength[0] + 1;
    freeSlipNeighbours[0] = xlength[0];
    #pragma omp for nowait schedule(static)
    for (i = 0; i < 5; i++)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
            {
                coordinate[1] = freeSlipNeighbours[1] = y;
                coordinate[2] = freeSlipNeighbours[2] = z;
                compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESLEFTOUT[i], coordinate, freeSlipNeighbours, VELOCITIESRIGHTOUT[i], inflowFeq);
            }

    /* front boundary, i.e. z = xlength + 1 */
    // i=0,1,2,3,4
    z = coordinate[2] = xlength[2] + 1;
    freeSlipNeighbours[2] = xlength[2];
    #pragma omp for nowait schedule(static)
    for (i = 0; i < 5; i++)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
            {
                coordinate[0] = freeSlipNeighbours[0] = x;
                coordinate[1] = freeSlipNeighbours[1] = y;
                compute_boundary(collideField, bddParams, flagField, xlength, VELOCITIESBACKOUT[i], coordinate, freeSlipNeighbours, VELOCITIESFRONTOUT[i], inflowFeq);
            }

    /* inner boundary cells
     * assumes inner boundary cells can only be NO_SLIP */
     #ifdef _ARBITRARYGEOMETRY_
     #pragma omp for nowait schedule(static)
    for (i = 0; i < PARAMQ; i++)
        for (z = 1; z <= xlength[2]; z++)
            for (y = 1; y <= xlength[1]; y++)
                for(x = 1; x <= xlength[0]; x++)
                    if (flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
                    {
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (!flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ])
                            {
                                currentCellIndex = z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x + i * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
                                neighbourCellIndex = neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX + (PARAMQ - 1 - i) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
    #endif // _ARBITRARYGEOMETRY_
}
}

/** boundary helper function */
void compute_boundary(double *collideField, const double * const bddParams, int *flagField, int *xlength, const int directionToUpdate, int * const coordinate, int * const freeSlipNeighbours, const int freeSlipDirectionToCopy, const double * const inflowFeq)
{
    int neighbourCoordX, neighbourCoordY, neighbourCoordZ;
    int currentCellIndex, neighbourCellIndex;
    double density;
    const double * const wallVelocity = bddParams + 4;
    double referenceDensity = 1.0;
    double feq[PARAMQ] __attribute__((aligned(32)));
    double velocity[3];
    const double pressureIn = bddParams[3];

    int boundaryType = flagField[coordinate[0] + (xlength[0] + 2) * coordinate[1] + (xlength[0] + 2) * (xlength[1] + 2) * coordinate[2]];

    switch (boundaryType)
    {
    case NO_SLIP:
        neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[directionToUpdate][0];
        neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[directionToUpdate][1];
        neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[directionToUpdate][2];

        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
        {
            currentCellIndex = (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + directionToUpdate * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            neighbourCellIndex = (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - directionToUpdate) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            collideField[currentCellIndex] = collideField[neighbourCellIndex];
        }
        break;
    case MOVING_WALL:
        neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[directionToUpdate][0];
        neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[directionToUpdate][1];
        neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[directionToUpdate][2];

        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
        {
            currentCellIndex = (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + directionToUpdate * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            neighbourCellIndex = (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - directionToUpdate) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            computeDensity(collideField + (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));

            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                              + 2 * LATTICEWEIGHTS[directionToUpdate] * (LATTICEVELOCITIES[directionToUpdate][0] * wallVelocity[0]
                                                      + LATTICEVELOCITIES[directionToUpdate][1] * wallVelocity[1] + LATTICEVELOCITIES[directionToUpdate][2] * wallVelocity[2]) * density / (C_S * C_S);
        }
        break;
    case INFLOW:
        neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[directionToUpdate][0];
        neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[directionToUpdate][1];
        neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[directionToUpdate][2];

        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
        {
            currentCellIndex = (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + directionToUpdate * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            neighbourCellIndex = (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - directionToUpdate) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            collideField[currentCellIndex] = inflowFeq[directionToUpdate];
        }

        break;

    case OUTFLOW:
        neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[directionToUpdate][0];
        neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[directionToUpdate][1];
        neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[directionToUpdate][2];

        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
        {
            currentCellIndex = (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + directionToUpdate * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            neighbourCellIndex = (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - directionToUpdate) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            computeDensity(collideField + (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
            computeVelocity(collideField + (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity, (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
            computeFeq(&referenceDensity, velocity, feq);
            collideField[currentCellIndex] = feq[PARAMQ - directionToUpdate - 1] + feq[directionToUpdate] - collideField[neighbourCellIndex];
        }
        break;
    case PRESSURE_IN:
        neighbourCoordX = coordinate[0] + LATTICEVELOCITIES[directionToUpdate][0];
        neighbourCoordY = coordinate[1] + LATTICEVELOCITIES[directionToUpdate][1];
        neighbourCoordZ = coordinate[2] + LATTICEVELOCITIES[directionToUpdate][2];

        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1
                && flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
        {
            currentCellIndex = (coordinate[2] * (xlength[0] + 2 ) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + directionToUpdate * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            neighbourCellIndex = (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - directionToUpdate) * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            computeDensity(collideField +  (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
            computeVelocity(collideField +  (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity, (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));
            density = referenceDensity + pressureIn;
            computeFeq(&density, velocity, feq);
            collideField[currentCellIndex] = feq[PARAMQ - directionToUpdate - 1] + feq[directionToUpdate] - collideField[neighbourCellIndex];
        }
        break;
    case FREE_SLIP:

        neighbourCoordX = freeSlipNeighbours[0];
        neighbourCoordY = freeSlipNeighbours[1];
        neighbourCoordZ = freeSlipNeighbours[2];
        if (!flagField[neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX])
        {
            currentCellIndex = (coordinate[2] * (xlength[0] + 2) * (xlength[1] + 2) + coordinate[1] * (xlength[0] + 2) + coordinate[0]) + directionToUpdate * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            neighbourCellIndex = (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + freeSlipDirectionToCopy * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
            collideField[currentCellIndex] = collideField[neighbourCellIndex];

        }

    }

}

