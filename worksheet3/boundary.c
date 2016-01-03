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
    int i, x, y, z, neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex;
    const double * const wallVelocity = bddParams + 4;
    const double pressureIn = bddParams[3];
    double density = 0.0;
    double velocity[3];
    double feq[PARAMQ];
    double referenceDensity = 1.0;

    /* top boundary
     * only need to set directions 0, 5, 6, 7, 14 */

    y = xlength[1] + 1;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (x = 0; x <= xlength[0] + 1; x++)
            switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
            {
                case NO_SLIP:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    break;
                }
                case MOVING_WALL:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
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
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
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
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 6;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y - 1) * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y - 1;
                        neighbourCoordZ = z;

                        i = 0;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 4;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 5;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 11;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 6;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 12;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 7;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 13;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 14;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 18;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y - 1) * (xlength[0] + 2) + x] == NO_SLIP)
                    {
                        // here the cell has to behave like a NO_SLIP cell
                        i = 0;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 5;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        // i = 6 is not needed since we already know that there is a bdd and not a fluid cell
                        i = 7;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 14;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
                    break;
                }
            }

    /* back boundary */
    z = 0;
    for (y = 0; y <= xlength[1] + 1; y++)
        for (x = 0; x <= xlength[0] + 1; x++)
            switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
            {
                case NO_SLIP:
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    break;
                }
                case MOVING_WALL:
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
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
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 14;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
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
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 16;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[(z + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y;
                        neighbourCoordZ = z + 1;

                        i = 14;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 0;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 15;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 1;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 16;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 2;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 17;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 3;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 18;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 4;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    else if (flagField[(z + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] == NO_SLIP)
                    {
                        // here the cell has to behave like a NO_SLIP cell
                        i = 14;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 15;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        // i = 16 is not needed since we already know that there is a bdd and not a fluid cell
                        i = 17;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 18;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
                    break;
                }
            }


    /* bottom boundary */
    y = 0;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (x = 0; x <= xlength[0] + 1; x++)
            switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
            {
                case NO_SLIP:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    break;
                }
                case MOVING_WALL:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
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
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 12;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 18;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y + 1) * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y + 1;
                        neighbourCoordZ = z;

                        i = 4;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 0;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 11;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 5;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 12;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 6;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 13;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 7;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 18;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 14;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (y + 1) * (xlength[0] + 2) + x] == NO_SLIP)
                    {
                        // here the cell has to behave like a NO_SLIP cell
                        i = 4;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 11;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        // i = 12 is not needed since we already know that there is a bdd and not a fluid cell
                        i = 13;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 18;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
                    break;
                }
            }


    /* left boundary */
    x = 0;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++)
            switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
            {
                case NO_SLIP:
                {
                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 10;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    break;
                }
                case MOVING_WALL:
                {
                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 10;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
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
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 10;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 10;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 7;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 10;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 13;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 17;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x + 1])
                    {
                        neighbourCoordX = x + 1;
                        neighbourCoordY = y;
                        neighbourCoordZ = z;

                        i = 3;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 1;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 7;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 5;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 10;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 8;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 13;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 11;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 17;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 15;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x + 1] == NO_SLIP)
                    {
                        // here the cell has to behave like a NO_SLIP cell
                        i = 3;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 7;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        // i = 10 is not needed since we already know that there is a bdd and not a fluid cell
                        i = 13;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 17;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
                    break;
                }
            }

    /* right boundary */
    x = xlength[0] + 1;
    for (z = 0; z <= xlength[2] + 1; z++)
        for (y = 0; y <= xlength[1] + 1; y++)
            switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
            {
                case NO_SLIP:
                {
                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 8;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    break;
                }
                case MOVING_WALL:
                {
                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 8;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
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
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 8;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 8;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case PRESSURE_IN:
                {
                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 5;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 8;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 11;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 15;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x - 1])
                    {
                        neighbourCoordX = x - 1;
                        neighbourCoordY = y;
                        neighbourCoordZ = z;

                        i = 1;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 3;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 5;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 7;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 8;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 10;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 11;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 13;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 15;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 17;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    else if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x - 1] == NO_SLIP)
                    {
                        // here the cell has to behave like a NO_SLIP cell
                        i = 1;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 5;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        // i = 8 is not needed since we already know that there is a bdd and not a fluid cell
                        i = 11;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 15;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
                    break;
                }
            }

    /* front boundary, i.e. z = xlength + 1 */
    z = xlength[2] + 1;
    for (y = 0; y <= xlength[1] + 1; y++)
        for (x = 0; x <= xlength[0] + 1; x++)
            switch(flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
            {
                case NO_SLIP:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            collideField[currentCellIndex] = collideField[neighbourCellIndex];
                        }
                    break;
                }
                case MOVING_WALL:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                              + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
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
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            computeFeq(&referenceDensity, bddParams, feq);
                            collideField[currentCellIndex] = feq[i];
                        }
                    break;
                }
                case OUTFLOW:
                {
                    i = 0;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            computeFeq(&referenceDensity, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
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
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 1;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 2;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 3;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }

                    i = 4;
                    neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                    neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                    neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                    if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                        if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                        {
                            currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                            neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                            density = 0.0;
                            computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
                            computeVelocitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density, velocity);
                            density = referenceDensity + pressureIn;
                            computeFeq(&density, velocity, feq);
                            collideField[currentCellIndex] = feq[PARAMQ - i - 1] + feq[i] - collideField[neighbourCellIndex];
                        }
                    break;
                }
                case FREE_SLIP:
                {
                    if (!flagField[(z - 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x])
                    {
                        neighbourCoordX = x;
                        neighbourCoordY = y ;
                        neighbourCoordZ = z - 1;

                        i = 0;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 14;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 1;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 15;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 2;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 16;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 3;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 17;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                        i = 4;
                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + 18;
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];

                    }
                    else if (flagField[(z - 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] == NO_SLIP)
                    {
                        // here the cell has to behave like a NO_SLIP cell
                        i = 0;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 1;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        // i = 2 is not needed since we already know that there is a bdd and not a fluid cell
                        i = 3;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }

                        i = 4;
                        neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                        neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                        neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                        if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
                            if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
                            {
                                currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
                                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                                collideField[currentCellIndex] = collideField[neighbourCellIndex];
                            }
                    }
                    break;
                }
            }


    /* inner boundary cells
     * assumes inner boundary cells can only be NO_SLIP */
    for (z = 1; z <= xlength[2]; z++)
        for (y = 1; y <= xlength[1]; y++)
            for(x = 1; x <= xlength[0]; x++)
                if (!!flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z])
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





//
//
//
//
//    /* top boundary */
//    y = xlength[1] + 1;
//    for (z = 0; z <= xlength[2] + 1; z++)
//        for (x = 0; x <= xlength[0] + 1; x++)
//            for (i = 0; i < PARAMQ; i++)
//            {
//                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
//                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
//                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
//                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength[0] + 1 && neighbourCoordY <= xlength[1] + 1 && neighbourCoordZ <= xlength[2] + 1)
//                    if (flagField[neighbourCoordX + (xlength[0] + 2) * neighbourCoordY + (xlength[0] + 2) * (xlength[1] + 2) * neighbourCoordZ] == 0)
//                    {
//                        currentCellIndex = PARAMQ * (z * (xlength[0] + 2 ) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
//                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
//                        density = 0.0;
//                        computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength[0] + 2) * (xlength[1] + 2) + neighbourCoordY * (xlength[0] + 2) + neighbourCoordX), &density);
//                        collideField[currentCellIndex] =  collideField[neighbourCellIndex]
//                                                          + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
//                    }
//            }



}

