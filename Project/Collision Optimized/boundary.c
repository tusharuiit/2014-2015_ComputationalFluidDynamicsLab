#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <omp.h>
#include "communication.h"

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

void treatBoundary(double *collideField, int* flagField, const double * const bddParams, int *xlength, const int iProc, const int jProc, const int kProc)
{
    int i, x, y, z, neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex;
    const double * const wallVelocity = bddParams + 4;
    double density = 0.0;
    int iCoord, jCoord, kCoord;

    computePosition(iProc, jProc, kProc, &iCoord, &jCoord, &kCoord);

#ifdef _ARBITRARYGEOMETRIES_
    const double pressureIn = bddParams[3];
    double velocity[3];
    double feq[PARAMQ];
    double referenceDensity = 1.0;
    double inflowFeq[PARAMQ];
    computeFeq(&referenceDensity, bddParams, inflowFeq);
#endif // _ARBITRARYGEOMETRIES_

#ifdef _ARBITRARYGEOMETRIES_
    #pragma omp parallel private(i, x, y, z, neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex, density, velocity, feq) shared(inflowFeq, referenceDensity, collideField, flagField, xlength, iCoord, jCoord, kCoord)
#else
    #pragma omp parallel private(i, x, y, z, neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex, density) shared(collideField, flagField, xlength, iCoord, jCoord, kCoord)
#endif // _ARBITRARYGEOMETRIES_
    {
#ifndef _ARBITRARYGEOMETRIES_
    if (iCoord == 0)
#endif // _ARBITRARYGEOMETRIES_
    {

        /* left boundary */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
#ifdef _ARBITRARYGEOMETRIES_
                switch(flagField[(xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
                {
                case NO_SLIP:
#endif // _ARBITRARYGEOMETRIES_
                {
                    i = 3;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 10;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
#ifdef _ARBITRARYGEOMETRIES_
                    break;
#endif // _ARBITRARYGEOMETRIES_
                }
#ifdef _ARBITRARYGEOMETRIES_
                case MOVING_WALL:
                {
                    i = 3;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 7;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 10;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 13;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 17;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }
                    break;
                }
                case INFLOW:
                {
                    i = 3;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 7;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 10;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 13;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 17;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 3;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 10;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 3;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 10;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = 1;
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + 1])
                    {
                        neighbourCoordX = 1;
                        neighbourCoordY = y;
                        neighbourCoordZ = z;

                        i = 3;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 1;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 7;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 5;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 10;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 8;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 13;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 11;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 17;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 15;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + 1] == NO_SLIP)
                    {
                        // here the cell has to behave like a NO_SLIP cell
                        i = 3;
                        neighbourCoordX = 1;
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 7;
                        neighbourCoordX = 1;
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        // i = 10 is not needed since we already know that there is a bdd and not a fluid cell
                        i = 13;
                        neighbourCoordX = 1;
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 17;
                        neighbourCoordX = 1;
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if (neighbourCoordY >= 0 && neighbourCoordZ >= 0 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2)) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
                    break;
                }
                }
#endif // _ARBITRARYGEOMETRIES_
    }

#ifndef _ARBITRARYGEOMETRIES_
    if (iCoord == iProc - 1)
#endif // _ARBITRARYGEOMETRIES_
    {
        /* right boundary */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (y = 0; y <= xlength[1] + 1; y++)
#ifdef _ARBITRARYGEOMETRIES_
                switch(flagField[(xlength[0] + 1) + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
                {
                case NO_SLIP:
#endif // _ARBITRARYGEOMETRIES_
                {
                    i = 1;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 8;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
#ifdef _ARBITRARYGEOMETRIES_
                    break;
#endif // _ARBITRARYGEOMETRIES_
                }
#ifdef _ARBITRARYGEOMETRIES_
                case MOVING_WALL:
                {
                    i = 1;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 5;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 8;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 11;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 15;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }
                    break;
                }
                case INFLOW:
                {
                    i = 1;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 5;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 8;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 11;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 15;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 1;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 8;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 1;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 8;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = xlength[0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0]])
                    {
                        neighbourCoordX = xlength[0];
                        neighbourCoordY = y;
                        neighbourCoordZ = z;

                        i = 1;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 3;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 5;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 7;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 8;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 10;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 11;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 13;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 15;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + (xlength[0] + 1)) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 17;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    break;
                }
                }
#endif // _ARBITRARYGEOMETRIES_
    }

#ifndef _ARBITRARYGEOMETRIES_
    if (kCoord == kProc - 1)
#endif // _ARBITRARYGEOMETRIES_
    {
        /* front boundary, i.e. z = xlength + 1 */
        #pragma omp for nowait schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
#ifdef _ARBITRARYGEOMETRIES_
                switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 1)])
                {
                case NO_SLIP:
#endif // _ARBITRARYGEOMETRIES_
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
#ifdef _ARBITRARYGEOMETRIES_
                    break;
#endif // _ARBITRARYGEOMETRIES_
                }
#ifdef _ARBITRARYGEOMETRIES_
                case MOVING_WALL:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }
                    break;
                }
                case INFLOW:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = xlength[2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[xlength[2] * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y ;
                        neighbourCoordZ = xlength[2];

                        i = 0;
                        currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 14;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 1;
                        currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 15;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 2;
                        currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 16;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 3;
                        currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 17;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 4;
                        currentCellIndex = PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 18;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    break;
                }
                }
#endif // _ARBITRARYGEOMETRIES_
    }

#ifndef _ARBITRARYGEOMETRIES_
    if (kCoord == 0)
#endif // _ARBITRARYGEOMETRIES_
    {
        /* back boundary */
        z = 0;
        #pragma omp for nowait schedule(static)
        for (y = 0; y <= xlength[1] + 1; y++)
            for (x = 0; x <= xlength[0] + 1; x++)
#ifdef _ARBITRARYGEOMETRIES_
                switch(flagField[x + (xlength[0] + 2) * y])
                {
                case NO_SLIP:
#endif // _ARBITRARYGEOMETRIES_
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
#ifdef _ARBITRARYGEOMETRIES_
                    break;
#endif // _ARBITRARYGEOMETRIES_
                }
#ifdef _ARBITRARYGEOMETRIES_
                case MOVING_WALL:
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }
                    break;
                }
                case INFLOW:
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = 1;
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[(xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y;
                        neighbourCoordZ = 1;

                        i = 14;
                        currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 0;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 15;
                        currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 1;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 16;
                        currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 2;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 17;
                        currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 3;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 18;
                        currentCellIndex = PARAMQ * (y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 4;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    break;
                }
                }
#endif // _ARBITRARYGEOMETRIES_
    }
#ifndef _ARBITRARYGEOMETRIES_
    if (jCoord == 0)
#endif // _ARBITRARYGEOMETRIES_
    {
        /* bottom boundary */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
#ifdef _ARBITRARYGEOMETRIES_
                switch(flagField[x + (xlength[0] + 2) * (xlength[1] + 2) * z])
                {
                case NO_SLIP:
#endif // _ARBITRARYGEOMETRIES_
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY =  LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
#ifdef _ARBITRARYGEOMETRIES_
                    break;
#endif // _ARBITRARYGEOMETRIES_
                }
#ifdef _ARBITRARYGEOMETRIES_
                case MOVING_WALL:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }
                    break;
                }
                case INFLOW:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = 1;
                        neighbourCoordZ = z;

                        i = 4;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 0;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 11;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 5;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 12;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 6;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 13;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 7;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 18;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 14;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    break;
                }
                }
#endif // _ARBITRARYGEOMETRIES_
    }

#ifndef _ARBITRARYGEOMETRIES_
    if (jCoord == jProc - 1)
#endif // _ARBITRARYGEOMETRIES_
    {
        /* top boundary
                * only need to set directions 0, 5, 6, 7, 14 */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= xlength[2] + 1; z++)
            for (x = 0; x <= xlength[0] + 1; x++)
#ifdef _ARBITRARYGEOMETRIES_
                switch(flagField[x + (xlength[0] + 2) * (xlength[1] + 1) + (xlength[0] + 2) * (xlength[1] + 2) * z])
                {
                case NO_SLIP:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = (xlength[1]);
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = (xlength[1]);
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    break;
                }
                case MOVING_WALL:
#endif // _ARBITRARYGEOMETRIES_
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
#ifdef _ARBITRARYGEOMETRIES_
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
#endif // _ARBITRARYGEOMETRIES_
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }
#ifdef _ARBITRARYGEOMETRIES_
                    break;
#endif // _ARBITRARYGEOMETRIES_
                }
#ifdef _ARBITRARYGEOMETRIES_
                case INFLOW:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = inflowFeq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = xlength[1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1  && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocity(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + xlength[1] * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = xlength[1];
                        neighbourCoordZ = z;

                        i = 0;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 4;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 5;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 11;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 6;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 12;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 7;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 13;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 14;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 18;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    break;
                }
                }
#endif // _ARBITRARYGEOMETRIES_
    }

        /* inner boundary cells
         * assumes inner boundary cells can only be NO_SLIP */
#ifdef _ARBITRARYGEOMETRIES_
        #pragma omp for nowait schedule(static)
        for (z = 1; z <= xlength[2]; z++)
            for (y = 1; y <= xlength[1]; y++)
                for(x = 1; x <= xlength[0]; x++)
                    if (flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
                        for (i = 0; i < PARAMQ; i++)
                        {
                            neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                            neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                            neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                            if (!flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ])
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        }
#endif
    } // end parallel region
}
