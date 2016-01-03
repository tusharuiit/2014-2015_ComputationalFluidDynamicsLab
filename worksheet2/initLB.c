#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[])
{
    /* TODO */
    if (argc == 2)
    {
        READ_INT(argv[1], *xlength);
        READ_DOUBLE(argv[1], *tau);
        read_double(argv[1], "velocityWallX", velocityWall);
        read_double(argv[1], "velocityWallY", velocityWall + 1);
        read_double(argv[1], "velocityWallZ", velocityWall + 2);
        READ_INT(argv[1], *timesteps);
        READ_INT(argv[1], *timestepsPerPlotting);
    }
    else

    {

        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return -1;
    }

    return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength)
{
    int i, x, y, z;

    for (z = 0; z <= xlength + 1; z++)
        for (y = 0; y <= xlength + 1; y++)
            for (x = 0; x <= xlength + 1; x++)
            {
                for (i = 0; i < PARAMQ; i++)
                {
                    collideField[PARAMQ * (z * (xlength + 2) * (xlength + 2) + y * (xlength + 2) + x) + i] = LATTICEWEIGHTS[i];
                    streamField[PARAMQ * (z * (xlength + 2) * (xlength + 2) + y * (xlength + 2) + x) + i] = LATTICEWEIGHTS[i];
                }
                flagField[z * (xlength + 2) * (xlength + 2) + y * (xlength + 2) + x] = FLUID;
            }

    /* back boundary */
    z = 0;
    for (y = 0; y <= xlength + 1; y++)
        for (x = 0; x <= xlength + 1; x++)
            flagField[z * (xlength + 2) * (xlength + 2) + (xlength + 2) * y + x] = NO_SLIP;

    /* bottom boundary */
    y = 0;
    for (z = 0; z <= xlength + 1; z++)
        for (x = 0; x <= xlength + 1; x++)
            flagField[z * (xlength + 2) * (xlength + 2) + (xlength + 2) * y + x] = NO_SLIP;

    /* left boundary */
    x = 0;
    for (z = 0; z <= xlength + 1; z++)
        for (y = 0; y <= xlength + 1; y++)
            flagField[z * (xlength + 2) * (xlength + 2) + (xlength + 2) * y + x] = NO_SLIP;

    /* right boundary */
    x = xlength + 1;
    for (z = 0; z <= xlength + 1; z++)
        for (y = 0; y <= xlength + 1; y++)
            flagField[z * (xlength + 2) * (xlength + 2) + (xlength + 2) * y + x] = NO_SLIP;

    /* front boundary, i.e. z = xlength + 1 */
    for (y = 0; y <= xlength + 1; y++)
        for (x = 0; x <= xlength + 1; x++)
            flagField[(xlength + 1) * (xlength + 2) * (xlength + 2) + y * (xlength + 2) + x] = NO_SLIP;

    /* top boundary */
    y = xlength + 1;
    for (z = 0; z <= xlength + 1; z++)
        for (x = 0; x <= xlength + 1; x++)
            flagField[z * (xlength + 2) * (xlength + 2)+ (xlength + 2) * y + x] = MOVING_WALL;

}

