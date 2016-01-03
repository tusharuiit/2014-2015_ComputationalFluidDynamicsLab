#include "initLB.h"
#include "LBDefinitions.h"
#include "communication.h"
#include <omp.h>

int readParameters( int *xlength, double *tau, double *bddParams, int *iProc, int *jProc, int *kProc,  int *timesteps, int *timestepsPerPlotting, char *problem, char *pgmInput, int argc, char *argv[] )
{
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
    if (argc == 2)
    {
        read_string(argv[1], "problem", problem);
        read_int(argv[1], "xlength", xlength);
        read_int(argv[1], "ylength", xlength + 1);
        read_int(argv[1], "zlength", xlength + 2);
        read_int(argv[1], "iProc", iProc);
        read_int(argv[1], "jProc", jProc);
        read_int(argv[1], "kProc" , kProc);
        READ_DOUBLE(argv[1], *tau);

        if(!strcmp(problem, "drivenCavity"))
        {
            read_double(argv[1], "velocityWallX", bddParams + 4);
            read_double(argv[1], "velocityWallY", bddParams + 5);
            read_double(argv[1], "velocityWallZ", bddParams + 6);
        }
        else
        {
            bddParams[4] = 0.0;
            bddParams[5] = 0.0;
            bddParams[6] = 0.0;
        }

        if(!strcmp(problem, "tiltedPlate"))
        {
            read_double(argv[1], "velocityInflowX", bddParams);
            read_double(argv[1], "velocityInflowY", bddParams + 1);
            read_double(argv[1], "velocityInflowZ", bddParams + 2);
            read_string(argv[1], "pgmInput", pgmInput);
        }
        else
        {
            bddParams[0] = 0.0;
            bddParams[1] = 0.0;
            bddParams[2] = 0.0;
            strcpy(pgmInput, "");
        }

        if (!strcmp(problem, "flowStep"))
        {
            read_double(argv[1], "velocityInflowX", bddParams);
            read_double(argv[1], "velocityInflowY", bddParams + 1);
            read_double(argv[1], "velocityInflowZ", bddParams + 2);
        }

        if(!strcmp(problem, "planeShearFlow"))
        {
            read_double(argv[1], "pressureIn", bddParams + 3);
        }
        else
        {
            bddParams[3] = 0.0;
        }
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


void initialiseFields( double *collideField, double *streamField, int *flagField, int * xlength, int *local_xlength, char *problem, char* pgmInput, int rank, int iProc, int jProc, int kProc )
{
    int iCoord, jCoord, kCoord, x, y, z, i;
#ifdef _ARBITRARYGEOMETRIES_
    int** pgmImage;
#endif // _ARBITRARYGEOMETRY_
    computePosition(iProc, jProc, kProc, &iCoord, &jCoord, &kCoord);
#ifdef _ARBITRARYGEOMETRIES_
    #pragma omp parallel private(x, y, z, i) shared(iCoord, jCoord, kCoord, collideField, flagField, streamField, xlength, problem, pgmInput, pgmImage, local_xlength, iProc, jProc, kProc)
#else
    #pragma omp parallel private(x, y, z, i) shared(iCoord, jCoord, kCoord, collideField, flagField, streamField, xlength, local_xlength, iProc, jProc, kProc)
#endif // _ARBITRARYGEOMETRY_
    {
        #pragma omp for schedule(static)
        for(x = 0; x < (local_xlength[0] + 2) * (local_xlength[1] + 2) * (local_xlength[2] + 2); x++)
        {
            flagField[x] = FLUID;
            for(i = 0; i < PARAMQ; i++)
                collideField[PARAMQ * x + i] = streamField[PARAMQ * x + i] = LATTICEWEIGHTS[i];
        }

        // independend of the problem first initialize all ghost cells as PARALLEL_BOUNDARY
        /* top boundary */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] = PARALLEL_BOUNDARY;

        /* back boundary */
        #pragma omp for nowait schedule(static)
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                flagField[y * (local_xlength[0] + 2) + x] = PARALLEL_BOUNDARY;

        /* bottom boundary */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] = PARALLEL_BOUNDARY;

        /* left boundary */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] = PARALLEL_BOUNDARY;

        /* right boundary */
        #pragma omp for nowait schedule(static)
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] = PARALLEL_BOUNDARY;

        /* front boundary, i.e. z = local_xlength + 1 */
        #pragma omp for schedule(static)
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                flagField[(local_xlength[2] + 1) * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = PARALLEL_BOUNDARY;

        /** initialization of different scenarios */
#ifdef _ARBITRARYGEOMETRIES_
        if (!strcmp(problem, "drivenCavity"))
        {
#endif // _ARBITRARYGEOMETRY_
            /* front boundary, i.e. z = local_xlength + 1 */
            if (kCoord == kProc - 1)
            {
                #pragma omp for nowait schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[(local_xlength[2] + 1) * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = NO_SLIP;
            }

            /* back boundary */
            if (kCoord == 0)
            {
                #pragma omp for schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[y * (local_xlength[0] + 2) + x] = NO_SLIP;
            }

            /* left boundary */
            if (iCoord == 0)
            {
                #pragma omp for nowait schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] = NO_SLIP;
            }
            /* right boundary */
            if (iCoord == iProc - 1)
            {
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] = NO_SLIP;
            }
            /* bottom boundary */
            if (jCoord == 0)
            {
                #pragma omp for nowait schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)

                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] = NO_SLIP;
            }

            /* top boundary */
            if (jCoord == jProc - 1)
            {
                #pragma omp for nowait schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] = MOVING_WALL;
            }
#ifdef _ARBITRARYGEOMETRIES_
        }
        if (!strcmp(problem, "tiltedPlate"))
        {
            #pragma omp single
            {
                pgmImage = read_pgm(pgmInput);
            }
            #pragma omp for nowait schedule(static)
            for (z = 1; z <= local_xlength[2]; z++)
                for (y = 1; y <= local_xlength[1]; y++)
                    for (x = 1; x <= local_xlength[0]; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = !!pgmImage[(xlength[0] / iProc) * iCoord + x][(xlength[1] / jProc) * jCoord + y];


            /* front boundary, i.e. z = local_xlength + 1 */
            if (kCoord == kProc - 1)
                #pragma omp for nowait schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                    {
                        // check for obstacle in the cell in the inner part of the domain adjacent to the boundary (i.e. NO_SLIP)
                        // if there is an obstacle, make the boundary NO_SLIP instead of FREE_SLIP as well
                        if (pgmImage[(xlength[0] / iProc) * iCoord + x][(xlength[1] / jProc) * jCoord + y])
                            flagField[(local_xlength[2] + 1) * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = NO_SLIP;
                        else
                            flagField[(local_xlength[2] + 1) * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = FREE_SLIP;
                    }

            /* back boundary */
            if (kCoord == 0)
                #pragma omp for schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                    {
                        // check for obstacle in the cell in the inner part of the domain adjacent to the boundary
                        // if there is an obstacle, make the boundary NO_SLIP instead of FREE_SLIP as well
                        if (pgmImage[(xlength[0] / iProc) * iCoord + x][(xlength[1] / jProc) * jCoord + y])
                            flagField[y * (local_xlength[0] + 2) + x] = NO_SLIP;
                        else
                            flagField[y * (local_xlength[0] + 2) + x] = FREE_SLIP;
                    }


            /* left boundary */
            if (iCoord == 0)
                #pragma omp for nowait schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] = INFLOW;

            /* right boundary */
            if (iCoord == iProc - 1)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] = OUTFLOW;


            /* top boundary */
            if (jCoord == jProc - 1)
                #pragma omp for nowait schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] = NO_SLIP;

            /* bottom boundary */
            if (jCoord == 0)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] = NO_SLIP;


            // Adjust PARALLEL_BOUNDARY if the "adjacent" cell in the domain of the neighbouring process is an obstacle cell
            /* top boundary */
            #pragma omp for nowait schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (x = 0; x <= local_xlength[0] + 1; x++)
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY && pgmImage[(xlength[0] / iProc) * iCoord + x][(xlength[1] / jProc) * (jCoord + 1) +  1])
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] = NO_SLIP;

            /* bottom boundary */
            #pragma omp for schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (x = 0; x <= local_xlength[0] + 1; x++)
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] == PARALLEL_BOUNDARY && pgmImage[(xlength[0] / iProc) * iCoord + x][(xlength[1] / jProc) * jCoord])
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] = NO_SLIP;

            /* left boundary */
            #pragma omp for nowait schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] == PARALLEL_BOUNDARY && pgmImage[(xlength[0] / iProc) * iCoord][(xlength[1] / jProc) * jCoord + y])
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] = NO_SLIP;

            /* right boundary */
            #pragma omp for schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] == PARALLEL_BOUNDARY && pgmImage[(xlength[0] / iProc) * (iCoord + 1) + 1][(xlength[1] / jProc) * jCoord + y])
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] = NO_SLIP;
            #pragma omp single
            {
                free_imatrix(pgmImage, 0, xlength[0] + 2, 0, xlength[1] + 2);
            }
        }
        if (!strcmp(problem, "flowStep"))
        {
            /* front boundary, i.e. z = local_xlength + 1 */
            if (kCoord == kProc - 1)
                #pragma omp for schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[(local_xlength[2] + 1) * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = NO_SLIP;

            /* back boundary */
            if (kCoord == 0)
                #pragma omp for schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[y * (local_xlength[0] + 2) + x] = NO_SLIP;

            /* left boundary */
            if (iCoord == 0)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] = INFLOW;

            /* right boundary */
            if (iCoord == iProc - 1)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] = OUTFLOW;

            /* top boundary */
            if (jCoord == jProc - 1)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] = NO_SLIP;

            /* bottom boundary */
            if (jCoord == 0)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] = NO_SLIP;


            /* step */
            // step height & width: jProc * local_xlength[1] / 2
            #pragma omp for schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (y = 1; y <= min(xlength[1] / 2 - jCoord * (xlength[1] / jProc), local_xlength[1] ) ; y++) /* integer division on purpose, half of the channel is blocked by step */
                    for (x = 1; x <= min(xlength[1] / 2 - iCoord * (xlength[0] / iProc), local_xlength[0] ); x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = NO_SLIP;

            // adjust the PARALLEL_BOUNDARY where "adjacent" cells in the neighbouring domain are NO_SLIP
            /* top boundary */
            #pragma omp for schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (x = 0; x <= local_xlength[0] + 1; x++)
                    if (x <= xlength[1] / 2 - iCoord * (xlength[0] / iProc) && xlength[1] / 2 - (jCoord + 1) * (xlength[1] / jProc) >= 1 && flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] = NO_SLIP;


            /* bottom boundary */
            #pragma omp for schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (x = 0; x <= local_xlength[0] + 1; x++)
                    if( x <= xlength[1] / 2 - iCoord * (xlength[0] / iProc) && 0 <= xlength[1] / 2 - jCoord * (xlength[1] / jProc) && flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] == PARALLEL_BOUNDARY)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] = NO_SLIP;

            /* left boundary */
            #pragma omp for schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    if( 0 <= xlength[1] / 2 - iCoord * (xlength[0] / iProc) && y <= xlength[1] / 2 - jCoord * (xlength[1] / jProc)  && flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] == PARALLEL_BOUNDARY)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] = NO_SLIP;

            /* right boundary */
            #pragma omp for schedule(static)
            for (z = 0; z <= local_xlength[2] + 1; z++)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    if( 1 <= xlength[1] / 2 - (iCoord + 1) * (xlength[0] / iProc) && y <= xlength[1] / 2 - jCoord * (xlength[1] / jProc)  && flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] == PARALLEL_BOUNDARY)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] = NO_SLIP;
        }
        if (!strcmp(problem, "planeShearFlow"))
        {

            /* front boundary, i.e. z = xlength + 1 */
            if (kCoord == kProc - 1)
                #pragma omp for schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[(local_xlength[2] + 1) * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] = FREE_SLIP;

            /* back boundary */
            if (kCoord == 0)
                #pragma omp for schedule(static)
                for (y = 0; y <= local_xlength[1] + 1; y++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[y * (local_xlength[0] + 2) + x] = FREE_SLIP;

            /* left boundary */
            if (iCoord == 0)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2)] = PRESSURE_IN;

            /* right boundary */
            if (iCoord == iProc - 1)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (y = 0; y <= local_xlength[1] + 1; y++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + local_xlength[0] + 1] = OUTFLOW;

            /* top boundary */
            if (jCoord == jProc - 1)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + (local_xlength[1] + 1) * (local_xlength[0] + 2) + x] = NO_SLIP;

            /* bottom boundary */
            if (jCoord == 0)
                #pragma omp for schedule(static)
                for (z = 0; z <= local_xlength[2] + 1; z++)
                    for (x = 0; x <= local_xlength[0] + 1; x++)
                        flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + x] = NO_SLIP;
        }
#endif // _ARBITRARYGEOMETRY_
    } // end of parallel region

    /** debugging code: checking the flagField */
#ifdef DEBUG
    int * exactFlagField;
    exactFlagField = (int *) malloc( (size_t) sizeof( int ) * (xlength[0] + 2) *  (xlength[1] + 2) * (xlength[2] + 2));
    FILE *fp2 = NULL;
    unsigned int line2 = 0;
    int noOfReadEntries;
    int error2 = 0;
    char szFileName2[80];
    computePosition(iProc, jProc, kProc, &iCoord, &jCoord, &kCoord);
    sprintf( szFileName2, "Testdata/%s/flagField.dat", problem );
    fp2 = fopen(szFileName2,"r");
    if (fp2 != NULL)
    {
        for (line2 = 0; line2 < (xlength[0] + 2) *  (xlength[1] + 2) * (xlength[2] + 2); line2++)
        {
            noOfReadEntries = fscanf(fp2,"%d",&exactFlagField[line2]);
            if (noOfReadEntries != 1)
                continue;
        }

    }
    fclose(fp2);
    for (z = 1; z <= local_xlength[2]; z++)
        for (y = 1; y <= local_xlength[1]; y++)
            for(x = 1; x <= local_xlength[0]; x++)
                for (i = 0; i < PARAMQ; i++)
                    if (flagField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1] / jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc)))])
                        error2 = 1;
    if (error2)
        printf("ERROR: Process %d has a different flagField at inner nodes.\n",rank);

    error2 = 0;
    // Check global boundaries as well
    if (iCoord == 0)
    {
        // Left boundary
        x = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                if (flagField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1] / jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc)))])
                    error2 = 1;

    }
    if (iCoord == iProc - 1)
    {
        // Right boundary
        x = local_xlength[0] + 1;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                if (flagField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1] / jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc)))])
                    error2 = 1;
    }
    if (jCoord == 0)
    {
        // Bottom boundary
        y = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (flagField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1] / jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc)))])
                    error2 = 1;

    }
    if (jCoord == jProc - 1)
    {
        // Top boundary
        y = local_xlength[1] + 1;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (flagField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1] / jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc)))])
                    error2 = 1;
    }
    if (kCoord == 0)
    {
        // back boundary
        z = 0;
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (flagField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1] / jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc)))])
                    error2 = 1;

    }
    if (kCoord == kProc - 1)
    {
        // front boundary
        z = local_xlength[2] + 1;
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (flagField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1] / jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc)))])
                    error2 = 1;
    }

    if (error2)
        printf("ERROR: Process %d has a different flagField at global boundary.\n",rank);
    free(exactFlagField);
#endif
    /** debugging code end */
}



