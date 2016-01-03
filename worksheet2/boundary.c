#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength)
{
    int i, x, y, z, neighbourCoordX, neighbourCoordY, neighbourCoordZ, currentCellIndex, neighbourCellIndex;

    /* top boundary */
    y = xlength + 1;
    for (z = 0; z <= xlength + 1; z++)
        for (x = 0; x <= xlength + 1; x++)
            for (i = 0; i < PARAMQ; i++)
            {
                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                currentCellIndex = PARAMQ * (z * (xlength + 2 ) * (xlength + 2) + y * (xlength + 2) + x) + i;
                neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength + 2) * (xlength + 2) + neighbourCoordY * (xlength + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength + 1 && neighbourCoordY <= xlength + 1 && neighbourCoordZ <= xlength + 1)
                    if (flagField[neighbourCoordX + (xlength + 2) * neighbourCoordY + (xlength + 2) * (xlength + 2) * neighbourCoordZ] == 0)
                    {
                        double density = 0.0;
                        computeDensitySSE(collideField + PARAMQ * (neighbourCoordZ * (xlength + 2) * (xlength + 2) + neighbourCoordY * (xlength + 2) + neighbourCoordX), &density);
                        collideField[currentCellIndex] =  collideField[neighbourCellIndex]
                                                          + 2 * LATTICEWEIGHTS[i] * (LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2]) * density / (C_S * C_S);
                    }
            }
    /* back boundary */
    z = 0;
    for (y = 0; y <= xlength + 1; y++)
        for (x = 0; x <= xlength + 1; x++)
            for (i = 0; i < PARAMQ; i++)
            {
                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength + 1 && neighbourCoordY <= xlength + 1 && neighbourCoordZ <= xlength + 1)
                    if (flagField[neighbourCoordX + (xlength + 2) * neighbourCoordY + (xlength + 2) * (xlength + 2) * neighbourCoordZ] == 0)
                    {
                        currentCellIndex = PARAMQ * (z * (xlength + 2 ) * (xlength + 2) + y * (xlength + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength + 2) * (xlength + 2) + neighbourCoordY * (xlength + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
                    }
            }


    /* bottom boundary */
    y = 0;
    for (z = 0; z <= xlength + 1; z++)
        for (x = 0; x <= xlength + 1; x++)
            for (i = 0; i < PARAMQ; i++)
            {
                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength + 1 && neighbourCoordY <= xlength + 1 && neighbourCoordZ <= xlength + 1)
                    if (flagField[neighbourCoordX + (xlength + 2) * neighbourCoordY + (xlength + 2) * (xlength + 2) * neighbourCoordZ] == 0)
                    {
                        currentCellIndex = PARAMQ * (z * (xlength + 2 ) * (xlength + 2) + y * (xlength + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength + 2) * (xlength + 2) + neighbourCoordY * (xlength + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
                    }
            }

    /* left boundary */
    x = 0;
    for (z = 0; z <= xlength + 1; z++)
        for (y = 0; y <= xlength + 1; y++)
            for (i = 0; i < PARAMQ; i++)
            {
                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength + 1 && neighbourCoordY <= xlength + 1 && neighbourCoordZ <= xlength + 1)
                    if (flagField[neighbourCoordX + (xlength + 2) * neighbourCoordY + (xlength + 2) * (xlength + 2) * neighbourCoordZ] == 0)
                    {
                        currentCellIndex = PARAMQ * (z * (xlength + 2 ) * (xlength + 2) + y * (xlength + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength + 2) * (xlength + 2) + neighbourCoordY * (xlength + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
                    }
            }

    /* right boundary */
    x = xlength + 1;
    for (z = 0; z <= xlength + 1; z++)
        for (y = 0; y <= xlength + 1; y++)
            for (i = 0; i < PARAMQ; i++)
            {
                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];
                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength + 1 && neighbourCoordY <= xlength + 1 && neighbourCoordZ <= xlength + 1)
                    if (flagField[neighbourCoordX + (xlength + 2) * neighbourCoordY + (xlength + 2) * (xlength + 2) * neighbourCoordZ] == 0)
                    {
                        currentCellIndex = PARAMQ * (z * (xlength + 2 ) * (xlength + 2) + y * (xlength + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength + 2) * (xlength + 2) + neighbourCoordY * (xlength + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
                    }
            }

    /* front boundary, i.e. z = xlength + 1 */
    z = xlength + 1;
    for (y = 0; y <= xlength + 1; y++)
        for (x = 0; x <= xlength + 1; x++)
            for (i = 0; i < PARAMQ; i++)
            {
                neighbourCoordX = x + LATTICEVELOCITIES[i][0];
                neighbourCoordY = y + LATTICEVELOCITIES[i][1];
                neighbourCoordZ = z + LATTICEVELOCITIES[i][2];

                if ((neighbourCoordX >= 0) && (neighbourCoordY >= 0) && (neighbourCoordZ >= 0) && neighbourCoordX <= xlength + 1 && neighbourCoordY <= xlength + 1 && neighbourCoordZ <= xlength + 1)
                    if (flagField[neighbourCoordX + (xlength + 2) * neighbourCoordY + (xlength + 2) * (xlength + 2) * neighbourCoordZ] == 0)
                    {
                        currentCellIndex = PARAMQ * (z * (xlength + 2 ) * (xlength + 2) + y * (xlength + 2) + x) + i;
                        neighbourCellIndex = PARAMQ * (neighbourCoordZ * (xlength + 2) * (xlength + 2) + neighbourCoordY * (xlength + 2) + neighbourCoordX) + (PARAMQ - 1 - i);
                        collideField[currentCellIndex] = collideField[neighbourCellIndex];
                    }
            }



}

