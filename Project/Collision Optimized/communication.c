#include "communication.h"
#include <mpi.h>
#include <stdlib.h>
#include "LBDefinitions.h"
#include <stdio.h>


void communicateBoundaryValues ( MPI_Datatype * MPISendTypes, MPI_Datatype * MPIRecvTypes, double * collideField, int * neighbours, MPI_Request * MPISendReq, MPI_Request * MPIRecvReq, int axis)
{
    // values of axis: 0 -> x direction, 1 -> y direction, 2 -> z direction
    if (neighbours[2 * axis ] >= 0)
    {
        MPI_Isend( collideField, 1, MPISendTypes[2 * axis], neighbours[2 * axis], 2 * axis + 1, MPI_COMM_WORLD, &MPISendReq[0]);
        MPI_Irecv( collideField, 1, MPIRecvTypes[2 * axis], neighbours[2 * axis], 2 * axis, MPI_COMM_WORLD, &MPIRecvReq[0]);
    }
    if (neighbours[2 * axis + 1] >= 0)
    {
        MPI_Isend( collideField, 1, MPISendTypes[2 * axis + 1], neighbours[2 * axis + 1], 2 * axis , MPI_COMM_WORLD, &MPISendReq[1]);
        MPI_Irecv( collideField, 1, MPIRecvTypes[2 * axis + 1], neighbours[2 * axis + 1], 2 * axis + 1, MPI_COMM_WORLD, &MPIRecvReq[1]);
    }
}
void computeLocalDomainSize(int * xlength, int * local_xlength, int iProc, int jProc, int kProc )
{
    int iCoord, jCoord, kCoord;
    computePosition( iProc, jProc, kProc, &iCoord, &jCoord, &kCoord);

    // The last process in each direction gets the largest chunk of the domain if it cannot be split evenly
    if (iCoord != iProc - 1)
        local_xlength[0] = xlength[0] / iProc;
    else
        local_xlength[0] = xlength[0] - (iProc - 1) * (xlength[0] / iProc);
    if (jCoord != jProc - 1)
        local_xlength[1] = xlength[1] / jProc;
    else
        local_xlength[1] = xlength[1] - (jProc - 1) * (xlength[1] / jProc);
    if (kCoord != kProc - 1)
        local_xlength[2] = xlength[2] / kProc;
    else
        local_xlength[2] = xlength[2] - (kProc - 1) * (xlength[2] / kProc);

    /** debugging code: local dimensions */
//    MPI_Barrier(MPI_COMM_WORLD);
//    printf( "Process %d has local dimension (%d, %d, %d)\n", rank, local_xlength[0], local_xlength[1], local_xlength[2] );
    /** end of debugging code */
}

void computeNeighbours( int iProc, int jProc, int kProc, int * neighbours )
{
    int iCoord, jCoord, kCoord, rank;
    //neighbours and buffers, [0: left, 1: right, 2: top, 3: bottom, 4: front, 5: back]
    neighbours[0] = neighbours[1] = neighbours[2] = neighbours[3] = neighbours[4] = neighbours[5] = MPI_PROC_NULL;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    computePosition(iProc, jProc, kProc, &iCoord, &jCoord, &kCoord);

    if (iCoord > 0)
        neighbours[0] = kCoord * iProc * jProc + jCoord * iProc + iCoord - 1;
    if (iCoord < iProc - 1)
        neighbours[1] = kCoord * iProc * jProc + jCoord * iProc + iCoord + 1;
    if (kCoord < kProc - 1)
        neighbours[4] = (kCoord + 1) * iProc * jProc + jCoord * iProc + iCoord;
    if (kCoord > 0)
        neighbours[5] = (kCoord - 1) * iProc * jProc + jCoord * iProc + iCoord;
    if (jCoord < jProc - 1)
        neighbours[2] = kCoord * iProc * jProc + (jCoord + 1) * iProc + iCoord;
    if (jCoord > 0)
        neighbours[3] = kCoord * iProc * jProc + (jCoord - 1) * iProc + iCoord;


    /** debugging code: correct neighbours */
//    MPI_Barrier( MPI_COMM_WORLD );
//    printf("DEBUG - I am process %d, my neighbours are %d, %d, %d, %d, %d, %d\n", rank, neighbours[0], neighbours[1], neighbours[2], neighbours[3], neighbours[4], neighbours[5]);
    /** end of debugging code */
}
void computePosition( int iProc, int jProc, int kProc, int * iCoord, int * jCoord, int * kCoord )
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    *iCoord = rank % iProc;
    *jCoord = ((rank - *iCoord) / iProc) % jProc;
    *kCoord = (rank - *iCoord - iProc **jCoord) / (iProc * jProc);

    /** debugging code: correct positions */
//    MPI_Barrier(MPI_COMM_WORLD);
//    printf("DEBUG - I am process %d, my coordinates are (%d, %d, %d)\n", rank, *iCoord, *jCoord, *kCoord);
    /** end of debugging code */
}

void getSendRecvCount( const int * const local_xlength, const int * const neighbours, const int * const flagField, int * sendCount, int * recvCount)
{
    int x, y, z, x_send, x_recv, y_send, y_recv, z_send, z_recv;
    sendCount[0] = sendCount[1] = sendCount[2] = sendCount[3] = sendCount[4] = sendCount[5] = 0;
    recvCount[0] = recvCount[1] = recvCount[2] = recvCount[3] = recvCount[4] = recvCount[5] = 0;

    if (neighbours[0] >= 0)
    {
        x_send = 1;
        x_recv = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                if (y == 0 || y == local_xlength[1] + 1 || z == 0 || z == local_xlength[2] + 1)
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FREE_SLIP)
                        sendCount[0] += 5;
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == FREE_SLIP)
                        recvCount[0] += 5;
                }
                else
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FLUID)
                        sendCount[0] += 5;
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == PARALLEL_BOUNDARY)
                        recvCount[0] += 5;
                }
    }

    if (neighbours[1] >= 0)
    {
        x_send = local_xlength[0];
        x_recv = local_xlength[0] + 1;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                if (y == 0 || y == local_xlength[1] + 1 || z == 0 || z == local_xlength[2] + 1)
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FREE_SLIP)
                        sendCount[1] += 5;
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == FREE_SLIP)
                        recvCount[1] += 5;
                }
                else
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FLUID)
                        sendCount[1] += 5;
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == PARALLEL_BOUNDARY)
                        recvCount[1] += 5;
                }
    }

    if (neighbours[2] >= 0)
    {
        y_send = local_xlength[1];
        y_recv = local_xlength[1] + 1;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
            {
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FLUID || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FREE_SLIP || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                    sendCount[2] += 5;
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == FREE_SLIP )
                    recvCount[2] += 5;
            }
    }

    if (neighbours[3] >= 0)
    {
        y_send = 1;
        y_recv = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
            {
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FLUID || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FREE_SLIP || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                    sendCount[3] += 5;
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == FREE_SLIP )
                    recvCount[3] += 5;
            }
    }

    if (neighbours[4] >= 0)
    {
        z_send = local_xlength[2];
        z_recv = local_xlength[2] + 1;
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (y == 0 || y == local_xlength[1] + 1)
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        sendCount[4] += 5;
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        recvCount[4] += 5;
                }
                else
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FLUID || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        sendCount[4] += 5;
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        recvCount[4] += 5;
                }
    }
    if (neighbours[5] >= 0)
    {
        z_send = 1;
        z_recv = 0;
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (y == 0 || y == local_xlength[1] + 1)
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        sendCount[5] += 5;
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        recvCount[5] += 5;
                }
                else
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FLUID || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        sendCount[5] += 5;
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        recvCount[5] += 5;
                }
    }
}

void getSendRecvIndices( const int * const local_xlength, const int * const neighbours, const int * const flagField, int ** sendIndices, int ** recvIndices )
{
    int x, y, z, i, x_send, x_recv, y_send, y_recv, z_send, z_recv, sendcounter, recvcounter;

    if (neighbours[0] >= 0)
    {
        x_send = 1;
        x_recv = 0;
        sendcounter = recvcounter = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                if (y == 0 || y == local_xlength[1] + 1 || z == 0 || z == local_xlength[2] + 1)
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[0][sendcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send) + VELOCITIESLEFTOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[0][recvcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv) + VELOCITIESRIGHTOUT[i];
                            recvcounter++;
                        }
                }
                else
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FLUID)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[0][sendcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send) + VELOCITIESLEFTOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == PARALLEL_BOUNDARY)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[0][recvcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv) + VELOCITIESRIGHTOUT[i];
                            recvcounter++;
                        }
                }
    }

    if (neighbours[1] >= 0)
    {
        x_send = local_xlength[0];
        x_recv = local_xlength[0] + 1;
        sendcounter = recvcounter = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (y = 0; y <= local_xlength[1] + 1; y++)
                if (y == 0 || y == local_xlength[1] + 1 || z == 0 || z == local_xlength[2] + 1)
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[1][sendcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send) + VELOCITIESRIGHTOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[1][recvcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv) + VELOCITIESLEFTOUT[i];
                            recvcounter++;
                        }
                }
                else
                {
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send] == FLUID)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[1][sendcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_send) + VELOCITIESRIGHTOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv] == PARALLEL_BOUNDARY)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[1][recvcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x_recv) + VELOCITIESLEFTOUT[i];
                            recvcounter++;
                        }
                }
    }

    if (neighbours[2] >= 0)
    {
        y_send = local_xlength[1];
        y_recv = local_xlength[1] + 1;
        sendcounter = recvcounter = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
            {
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FLUID || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FREE_SLIP || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                    for (i = 0; i < 5; i++)
                    {
                        sendIndices[2][sendcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x) + VELOCITIESTOPOUT[i];
                        sendcounter++;
                    }
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY ||
                    flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == FREE_SLIP )
                    for (i = 0; i < 5; i++)
                    {
                        recvIndices[2][recvcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x) + VELOCITIESBOTTOMOUT[i];
                        recvcounter++;
                    }
            }
    }

    if (neighbours[3] >= 0)
    {
        y_send = 1;
        y_recv = 0;
        sendcounter = recvcounter = 0;
        for (z = 0; z <= local_xlength[2] + 1; z++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
            {
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FLUID || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == FREE_SLIP || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY)
                    for (i = 0; i < 5; i++)
                    {
                        sendIndices[3][sendcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_send * (local_xlength[0] + 2) + x) + VELOCITIESBOTTOMOUT[i];
                        sendcounter++;
                    }
                if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x] == FREE_SLIP )
                    for (i = 0; i < 5; i++)
                    {
                        recvIndices[3][recvcounter] = PARAMQ * (z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y_recv * (local_xlength[0] + 2) + x) + VELOCITIESTOPOUT[i];
                        recvcounter++;
                    }
            }
    }

    if (neighbours[4] >= 0)
    {
        z_send = local_xlength[2];
        z_recv = local_xlength[2] + 1;
        sendcounter = recvcounter = 0;
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (y == 0 || y == local_xlength[1] + 1)
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[4][sendcounter] = PARAMQ * (z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[4][recvcounter] = PARAMQ * (z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESBACKOUT[i];
                            recvcounter++;
                        }
                }
                else
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FLUID || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[4][sendcounter] = PARAMQ * (z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[4][recvcounter] = PARAMQ * (z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESBACKOUT[i];
                            recvcounter++;
                        }
                }
    }
    if (neighbours[5] >= 0)
    {
        z_send = 1;
        z_recv = 0;
        sendcounter = recvcounter = 0;
        for (y = 0; y <= local_xlength[1] + 1; y++)
            for (x = 0; x <= local_xlength[0] + 1; x++)
                if (y == 0 || y == local_xlength[1] + 1)
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[5][sendcounter] = PARAMQ * (z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESBACKOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[5][recvcounter] = PARAMQ * (z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i];
                            recvcounter++;
                        }
                }
                else
                {
                    if (flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FLUID || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            sendIndices[5][sendcounter] = PARAMQ * (z_send * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESBACKOUT[i];
                            sendcounter++;
                        }
                    if (flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == PARALLEL_BOUNDARY || flagField[z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FREE_SLIP)
                        for (i = 0; i < 5; i++)
                        {
                            recvIndices[5][recvcounter] = PARAMQ * (z_recv * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + VELOCITIESFRONTOUT[i];
                            recvcounter++;
                        }
                }
    }

}

void checkRequestCompletion( MPI_Request * MPISendReq, MPI_Request * MPIRecvReq, int axis, int * neighbours)
{

    if (neighbours[2 * axis] >= 0)
    {
        MPI_Wait(&MPISendReq[0], MPI_STATUS_IGNORE);
        MPI_Wait(&MPIRecvReq[0], MPI_STATUS_IGNORE);
    }
    if (neighbours[2 * axis + 1] >= 0)
    {
        MPI_Wait(&MPISendReq[1], MPI_STATUS_IGNORE);
        MPI_Wait(&MPIRecvReq[1], MPI_STATUS_IGNORE);
    }

}

void initialiseMPITypes( const int * const local_xlength, const int * const neighbours, const int * const flagField, MPI_Datatype * MPISendTypes, MPI_Datatype * MPIRecvTypes )
{
    int sendCount[6], recvCount[6];
    int *sendIndices[6], *recvIndices[6];
    int i;

    getSendRecvCount(local_xlength, neighbours, flagField, sendCount, recvCount);

    for (i = 0; i < 6; i++)
    {
        if(neighbours[i] >= 0)
        {
            sendIndices[i] = (int*) malloc((size_t) sendCount[i] * sizeof(int));
            recvIndices[i] = (int*) malloc((size_t) recvCount[i] * sizeof(int));
        }
    }

    getSendRecvIndices(local_xlength, neighbours, flagField, sendIndices, recvIndices);
    // build self defined types
    for(int i = 0; i < 6; i++)
        if(neighbours[i] >= 0)
        {
            int * blocksizeSend = (int *) malloc( (size_t) sendCount[i] * sizeof( int ) );
            for (int j = 0; j < sendCount[i]; j++)
                blocksizeSend[j] = 1;
            int * blocksizeRecv = (int *) malloc ( (size_t) recvCount[i] * sizeof( int ) );
            for (int j = 0; j < recvCount[i]; j++)
                blocksizeRecv[j] = 1;

            MPI_Type_indexed(sendCount[i], blocksizeSend, sendIndices[i], MPI_DOUBLE, &MPISendTypes[i]);
            MPI_Type_commit( &MPISendTypes[i] );

            MPI_Type_indexed(recvCount[i], blocksizeRecv, recvIndices[i], MPI_DOUBLE, &MPIRecvTypes[i]);
            MPI_Type_commit( &MPIRecvTypes[i] );

            free(blocksizeSend);
            free(blocksizeRecv);
            free(sendIndices[i]);
            free(recvIndices[i]);
        }
}

