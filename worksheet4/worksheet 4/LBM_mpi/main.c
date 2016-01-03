#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include <time.h>
#include <unistd.h>
#include <mpi.h>
#include "parallel.h"

int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    char problem[100];
    char pgmInput[1000];
    int *flagField = NULL;
    clock_t begin, end;
    double time_spent;

    int xlength[3], timesteps, timestepsPerPlotting;
    double tau, bddParams[7];

    int rank;
    int number_of_ranks;
    int iProc, jProc, kProc;
    int leftNeighbour, rightNeighbour, bottomNeighbour, topNeighbour, backNeighbour, frontNeighbour;
    leftNeighbour = rightNeighbour = bottomNeighbour = topNeighbour = backNeighbour = frontNeighbour = MPI_PROC_NULL;
    MPI_Status mpistatus;



    int x, y, z, i;

    // send and read buffers for all possible directions :
    // [0: left, 1: right, 2: top, 3: bottom, 4: front, 5: back]
    double *sendBuffer[6];
    double *readBuffer[6];

    sendBuffer[0] = sendBuffer[1] = sendBuffer[2] = sendBuffer[3] = sendBuffer[4] = sendBuffer[5] = NULL;
    readBuffer[0] = readBuffer[1] = readBuffer[2] = readBuffer[3] = readBuffer[4] = readBuffer[5] = NULL;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );      /* asking for the number of processes  */
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );                 /* asking for the local process id   */

    if(readParameters(xlength, &tau, bddParams, &iProc, &jProc, &kProc, &timesteps, &timestepsPerPlotting, problem, pgmInput, argc, argv) == 0)
    {
        if (number_of_ranks != iProc * jProc * kProc)
        {
            if (rank == 0)
                printf("ERROR: number of processes started does not match the number specified in the input file!\n");
            MPI_Barrier( MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, 0);
        }
        begin = clock();

        time_t start = time(NULL);
        collideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        streamField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));
        flagField = (int *) malloc((size_t) sizeof (int) * (xlength[0] + 2)*(xlength[1] + 2)*(xlength[2] + 2));

        int iCoord, jCoord, kCoord; // position of the domain in the decomposition; needed to compute neighbouring processes and for vtk output
        iCoord = rank % iProc;
        jCoord = ((rank - iCoord) / iProc) % jProc;
        kCoord = (rank - iCoord - iProc * jCoord) / (iProc * jProc);

        if (iCoord > 0)
        {
            leftNeighbour = kCoord * iProc * jProc + jCoord * iProc + iCoord - 1;
            sendBuffer[0] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[1] + 2) * (xlength[2] + 2));
            readBuffer[1] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[1] + 2) * (xlength[2] + 2));
        }
        if (iCoord < iProc - 1)
        {
            rightNeighbour = kCoord * iProc * jProc + jCoord * iProc + iCoord + 1;
            sendBuffer[1] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[1] + 2) * (xlength[2] + 2));
            readBuffer[0] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[1] + 2) * (xlength[2] + 2));
        }
        if (kCoord < kProc - 1)
        {
            frontNeighbour = (kCoord + 1) * iProc * jProc + jCoord * iProc + iCoord;
            sendBuffer[4] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[1] + 2));
            readBuffer[5] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[1] + 2));
        }
        if (kCoord > 0)
        {
            backNeighbour = (kCoord - 1) * iProc * jProc + jCoord * iProc + iCoord;
            sendBuffer[5] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[1] + 2));
            readBuffer[4] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[1] + 2));
        }
        if (jCoord < jProc - 1)
        {
            topNeighbour = kCoord * iProc * jProc + (jCoord + 1) * iProc + iCoord;
            sendBuffer[2] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[2] + 2));
            readBuffer[3] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[2] + 2));
        }
        if (jCoord > 0)
        {
            bottomNeighbour = kCoord * iProc * jProc + (jCoord - 1) * iProc + iCoord;
            sendBuffer[3] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[2] + 2));
            readBuffer[2] = (double*) malloc((size_t) sizeof(double) * 5 * (xlength[0] + 2) * (xlength[2] + 2));
        }

        /* Debug only
        MPI_Barrier(MPI_COMM_WORLD);
        printf("DEBUG - I am process %d, my coordinates are (%d, %d, %d)\n", rank, iCoord, jCoord, kCoord);
        MPI_Barrier(MPI_COMM_WORLD);
        printf("DEBUG - I am process %d, my neighbours are %d, %d, %d, %d, %d, %d\n", rank, leftNeighbour, rightNeighbour, topNeighbour, bottomNeighbour, frontNeighbour, backNeighbour);
        MPI_Barrier(MPI_COMM_WORLD); */


        initialiseFields(collideField, streamField, flagField, xlength, problem, pgmInput, rank, iProc, jProc, kProc);

        /** debugging code, checking the flagField */
//        int exactFlagField[(iProc * xlength[0] + 2) *  (jProc * xlength[1] + 2) * (kProc * xlength[2] + 2)];
//        FILE *fp2 = NULL;
//        unsigned int line2 = 0;
//        int error2 = 0;
//        char szFileName2[80];
//        sprintf( szFileName2, "Testdata/%s/flagField.dat", problem );
//        fp2 = fopen(szFileName2,"r");
//        if (fp2 != NULL)
//        {
//            for (line2 = 0; line2 < (iProc * xlength[0] + 2) *  (jProc * xlength[1] + 2) * (kProc * xlength[2] + 2); line2++)
//                fscanf(fp2,"%d",&exactFlagField[line2]);
//        }
//        fclose(fp2);
//        for (z = 1; z <= xlength[2]; z++)
//            for (y = 1; y <= xlength[1]; y++)
//                for(x = 1; x <= xlength[0]; x++)
//                    for (i = 0; i < PARAMQ; i++)
//                        if (flagField[(z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x)] != exactFlagField[((z + kCoord * xlength[2]) * (xlength[0] * iProc + 2) * (xlength[1] * jProc + 2) + (y + jCoord * xlength[1]) * (xlength[0] * iProc + 2) + (x + iCoord * xlength[0]))])
//                            error2 = 1;
//
//
//        if (error2)
//            printf("ERROR: Process %d has a different flagField\n",rank);
        /** debugging code end */

        if (!rank)
            printf("Progress:     ");
        for(int t = 0; t < timesteps; t++)
        {
            double *swap = NULL;

            if (rightNeighbour >= 0 && leftNeighbour >= 0)
            {
                // both left and right neighbours
                // Extraction at right boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        for (i = 0; i < 5; i++)
                            sendBuffer[1][5 * ((xlength[1] + 2) * z + y) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0]) + VELOCITIESRIGHTOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[1], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, rightNeighbour, 1, readBuffer[1], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, leftNeighbour, 1, MPI_COMM_WORLD, &mpistatus );
                // Injection at left boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)) + VELOCITIESRIGHTOUT[i]] = readBuffer[1][5 * ((xlength[1] + 2) * z + y) + i];
                // Extraction at left boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        for (i = 0; i < 5; i++)
                            sendBuffer[0][5 * ((xlength[1] + 2) * z + y) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + 1) + VELOCITIESLEFTOUT[i]];
                MPI_Sendrecv( sendBuffer[0], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, leftNeighbour, 0, readBuffer[0], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, rightNeighbour, 0, MPI_COMM_WORLD, &mpistatus );
                // Injection at right boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1) + VELOCITIESLEFTOUT[i]] = readBuffer[0][5 * ((xlength[1] + 2) * z + y) + i];

            }
            else if (rightNeighbour >= 0)
            {
                // no left neighbour
                // Extraction at right boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        for (i = 0; i < 5; i++)
                            sendBuffer[1][5 * ((xlength[1] + 2) * z + y) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0]) + VELOCITIESRIGHTOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[1], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, rightNeighbour, 1, readBuffer[0], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, rightNeighbour, 0, MPI_COMM_WORLD, &mpistatus );
                // Injection at right boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + xlength[0] + 1) + VELOCITIESLEFTOUT[i]] = readBuffer[0][5 * ((xlength[1] + 2) * z + y) + i];
            }
            else if (leftNeighbour >= 0)
            {
                // no right neighbour
                // Extraction at left boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        for (i = 0; i < 5; i++)
                            sendBuffer[0][5 * ((xlength[1] + 2) * z + y) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + 1) + VELOCITIESLEFTOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[0], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, leftNeighbour, 0, readBuffer[1], 5 * (xlength[1] + 2) * (xlength[2] + 2), MPI_DOUBLE, leftNeighbour, 1, MPI_COMM_WORLD, &mpistatus );
                // Injection at left boundary
                for (z = 1; z <= xlength[2]; z++)
                    for (y = 1; y <= xlength[1]; y++)
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2)) + VELOCITIESRIGHTOUT[i]] = readBuffer[1][5 * ((xlength[1] + 2) * z + y) + i];
            }

            if (frontNeighbour >= 0 && backNeighbour >= 0)
            {
                // Extraction at front boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for (i = 0; i < 5; i++)
                            sendBuffer[4][5 * ((xlength[0] + 2) * y + x) + i] = collideField[PARAMQ * (xlength[2] * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[4], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, frontNeighbour, 4, readBuffer[4], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, backNeighbour, 4,  MPI_COMM_WORLD, &mpistatus );
                // Injection at back boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[y * (xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * (y * (xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i]] = readBuffer[4][5 * ((xlength[0] + 2) * y + x) + i];
                // Extraction at back boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for (i = 0; i < 5; i++)
                            sendBuffer[5][5 * ((xlength[0] + 2) * y + x) + i] = collideField[PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + VELOCITIESBACKOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[5], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, backNeighbour, 5, readBuffer[5], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, frontNeighbour, 5,  MPI_COMM_WORLD, &mpistatus );
                // Injection at front boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[(xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + VELOCITIESBACKOUT[i]] = readBuffer[5][5 * ((xlength[0] + 2) * y + x) + i];

            }
            else if (frontNeighbour >= 0)
            {
                // Extraction at front boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for (i = 0; i < 5; i++)
                            sendBuffer[4][5 * ((xlength[0] + 2) * y + x) + i] = collideField[PARAMQ * (xlength[2] * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[4], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, frontNeighbour, 4, readBuffer[5], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, frontNeighbour, 5,  MPI_COMM_WORLD, &mpistatus );
                // Injection at front boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[(xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * ((xlength[2] + 1) * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + VELOCITIESBACKOUT[i]] = readBuffer[5][5 * ((xlength[0] + 2) * y + x) + i];
            }
            else if (backNeighbour >= 0)
            {
                // Extraction at back boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for (i = 0; i < 5; i++)
                            sendBuffer[5][5 * ((xlength[0] + 2) * y + x) + i] = collideField[PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + VELOCITIESBACKOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[5], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, backNeighbour, 5, readBuffer[4], 5 * (xlength[0] + 2) * (xlength[1] + 2), MPI_DOUBLE, backNeighbour, 4,  MPI_COMM_WORLD, &mpistatus );
                // Injection at back boundary
                for ( y = 1; y <= xlength[1]; y++)
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[y * (xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                            for (i = 0; i < 5; i++)
                                collideField[PARAMQ * (y * (xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i]] = readBuffer[4][5 * ((xlength[0] + 2) * y + x) + i];
            }
//
            if (topNeighbour >= 0 && bottomNeighbour >= 0)
            {
                // Extraction at top boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for ( i = 0; i < 5; i++ )
                            sendBuffer[2][5 * ( z * (xlength[0] + 2) + x) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + xlength[1] * (xlength[0] + 2) + x) + VELOCITIESTOPOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[2], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, topNeighbour, 2, readBuffer[2], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, bottomNeighbour, 2, MPI_COMM_WORLD, &mpistatus );
                // Injection at bottom boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + x] == PARALLEL_BOUNDARY)
                            for ( i = 0; i < 5; i++ )
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + x) + VELOCITIESTOPOUT[i]] = readBuffer[2][5 * ( z * (xlength[0] + 2) + x) + i];
                // Extraction at bottom boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for ( i = 0; i < 5; i++ )
                            sendBuffer[3][5 * ( z * (xlength[0] + 2) + x) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[0] + 2) + x) + VELOCITIESBOTTOMOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[3], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, bottomNeighbour, 3, readBuffer[3], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, topNeighbour, 3, MPI_COMM_WORLD, &mpistatus );
                // Injection at top boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                            for ( i = 0; i < 5; i++ )
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + VELOCITIESBOTTOMOUT[i]] = readBuffer[3][5 * ( z * (xlength[0] + 2) + x) + i];

            }
            else if (topNeighbour >= 0)
            {
                // Extraction at top boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for ( i = 0; i < 5; i++ )
                            sendBuffer[2][5 * ( z * (xlength[0] + 2) + x) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + xlength[1] * (xlength[0] + 2) + x) + VELOCITIESTOPOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[2], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, topNeighbour, 2, readBuffer[3], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, topNeighbour, 3, MPI_COMM_WORLD, &mpistatus );
                // Injection at top boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                            for ( i = 0; i < 5; i++ )
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[1] + 1) * (xlength[0] + 2) + x) + VELOCITIESBOTTOMOUT[i]] = readBuffer[3][5 * ( z * (xlength[0] + 2) + x) + i];

            }
            else if (bottomNeighbour >= 0)
            {
                // Extraction at bottom boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        for ( i = 0; i < 5; i++ )
                            sendBuffer[3][5 * ( z * (xlength[0] + 2) + x) + i] = collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[0] + 2) + x) + VELOCITIESBOTTOMOUT[i]];
                // Swap
                MPI_Sendrecv( sendBuffer[3], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, bottomNeighbour, 3, readBuffer[2], 5 * (xlength[0] + 2) * (xlength[2] + 2), MPI_DOUBLE, bottomNeighbour, 2, MPI_COMM_WORLD, &mpistatus );
                // Injection at bottom boundary
                for (z = 0; z <= xlength[2] + 1; z++ )
                    for (x = 0; x <= xlength[0] + 1; x++ )
                        if (flagField[z * (xlength[0] + 2) * (xlength[1] + 2) + x] == PARALLEL_BOUNDARY)
                            for ( i = 0; i < 5; i++ )
                                collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + x) + VELOCITIESTOPOUT[i]] = readBuffer[2][5 * ( z * (xlength[0] + 2) + x) + i];

            }


            doStreaming(collideField, streamField, flagField, xlength);
            swap = collideField;
            collideField = streamField;
            streamField = swap;

            doCollision(collideField, flagField, &tau, xlength);
            treatBoundary(collideField, flagField, bddParams, xlength);

            if (t % timestepsPerPlotting == 0)
                writeVtkOutput(collideField, flagField, "./Paraview/output", (unsigned int) t / timestepsPerPlotting, xlength, rank, iCoord, jCoord, kCoord);

            if (!rank)
            {
                int pct = ((float) t / timesteps) * 100;
                printf("\b\b\b%02d%%", pct);
                fflush(stdout);
            }

            /** debugging code */
            /* check correctness of collideField with reference data */
//            if (t % timestepsPerPlotting == 0)
//            {
//                double exactCollideField[PARAMQ * (iProc * xlength[0] + 2) *  (jProc * xlength[1] + 2) * (kProc * xlength[2] + 2)];
//                FILE *fp = NULL;
//                unsigned int line = 0;
//                int error = 0;
//                char szFileName[80];
//                sprintf( szFileName, "Testdata/%s/%i.dat", problem, t / timestepsPerPlotting );
//                fp = fopen(szFileName,"r");
//                if (fp != NULL)
//                {
//                    for (line = 0; line < PARAMQ * (iProc * xlength[0] + 2) *  (jProc * xlength[1] + 2) * (kProc * xlength[2] + 2); line++)
//                        fscanf(fp,"%lf",&exactCollideField[line]);
//                }
//                fclose(fp);
//                for (z = 1; z <= xlength[2]; z++)
//                    for (y = 1; y <= xlength[1]; y++)
//                        for(x = 1; x <= xlength[0]; x++)
//                            for (i = 0; i < PARAMQ; i++)
//                                if (fabs(collideField[PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i] - exactCollideField[PARAMQ * ((z + kCoord * xlength[2]) * (xlength[0] * iProc + 2) * (xlength[1] * jProc + 2) + (y + jCoord * xlength[1]) * (xlength[0] * iProc + 2) + (x + iCoord * xlength[0])) + i]) > 1e-5)
//                                    error = 1;
//                if (error)
//                    printf("ERROR: Process %d has a different collideField in timestep %d\n",rank, t);
//            }
            /** end of debugging code */


        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (!rank)
        {
            printf("\b\b\b\b100%%\n");
            end = clock();
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

            printf("Running time (CPU time): %.2fs\n", time_spent);
            printf("Running time (Wall clock): %2.fs\n", (double)(time(NULL) - start) );
            printf("MLUPS: %.3f\n", ((double) (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2) * timesteps * number_of_ranks) / (1000000.0 * time_spent));
        }

        free(collideField);
        free(streamField);
        free(flagField);
        free(sendBuffer[0]);
        free(sendBuffer[1]);
        free(sendBuffer[2]);
        free(sendBuffer[3]);
        free(sendBuffer[4]);
        free(sendBuffer[5]);
        free(readBuffer[0]);
        free(readBuffer[1]);
        free(readBuffer[2]);
        free(readBuffer[3]);
        free(readBuffer[4]);
        free(readBuffer[5]);
    }

    MPI_Finalize();
    return 0;
}

#endif

