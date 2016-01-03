#ifndef _MAIN_C_
#define _MAIN_C_

#include <mpi.h>
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include "communication.h"
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <sys/time.h>
//#include <papi.h>

int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    char problem[100];
    char pgmInput[1000];
    int *flagField = NULL;
    clock_t begin, end;
    double time_spent;
    struct timeval time_start, time_end;

//    long long counters[3];
//    int PAPI_events[] =
//    {
//        PAPI_TOT_CYC,
//        PAPI_L2_DCM,
//        PAPI_L2_DCA
//    };
//
//    PAPI_library_init(PAPI_VER_CURRENT);

    int xlength[3], local_xlength[3], timesteps, timestepsPerPlotting;
    double tau, bddParams[7];

    int rank;
    int number_of_ranks;
    int iProc, jProc, kProc;

    int neighbours[6];  // [0: left, 1: right, 2: top, 3: bottom, 4: front, 5: back]

    // send and read MPI datatypes for all possible directions :
    // [0: left, 1: right, 2: top, 3: bottom, 4: front, 5: back]
    MPI_Datatype MPISendTypes[6], MPIRecvTypes[6];
    MPI_Request MPISendReq[2], MPIRecvReq[2];

#ifdef _DEBUG_
    double * exactCollideField; // for debugging only
#endif



    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );      /* asking for the number of processes  */
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );                 /* asking for the local process id   */

    MPI_File fh;
    int err = MPI_File_open( MPI_COMM_WORLD, "VTK.bin.data", MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
    if (err) {
        fprintf(stderr, "Can't open file for writing.\n");
        MPI_Abort( MPI_COMM_WORLD, 911 );
        MPI_Finalize();
    }

    if(readParameters( xlength, &tau, bddParams, &iProc, &jProc, &kProc, &timesteps, &timestepsPerPlotting, problem, pgmInput, argc, argv ) == 0)
    {
        if (number_of_ranks != iProc * jProc * kProc)
        {
            if (rank == 0)
                printf("ERROR: number of processes started does not match the number specified in the input file!\n");
            MPI_Barrier( MPI_COMM_WORLD);
            MPI_Finalize();
            return -1;
        }
        begin = clock();
        gettimeofday(&time_start, NULL);


        int iCoord, jCoord, kCoord; // position of the domain in the decomposition; needed to compute neighbouring processes and for vtk output
        computePosition( iProc, jProc, kProc, &iCoord, &jCoord, &kCoord );
        computeLocalDomainSize( xlength, local_xlength, iProc, jProc, kProc );

        collideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (local_xlength[0] + 2)*(local_xlength[1] + 2)*(local_xlength[2] + 2));
        streamField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (local_xlength[0] + 2)*(local_xlength[1] + 2)*(local_xlength[2] + 2));
        flagField = (int *) malloc((size_t) sizeof (int) * (local_xlength[0] + 2)*(local_xlength[1] + 2)*(local_xlength[2] + 2));

        initialiseFields( collideField, streamField, flagField, xlength, local_xlength, problem, pgmInput, rank, iProc, jProc, kProc );
        computeNeighbours( iProc, jProc, kProc, neighbours );
        initialiseMPITypes(local_xlength, neighbours, flagField, MPISendTypes, MPIRecvTypes);


#ifndef _NOPROGRESS_
        if (!rank)
            printf("Progress:     ");
#endif // _NOPROGRESS_
//        PAPI_start_counters( PAPI_events, 3 );
        for(int t = 0; t < timesteps; t++)
        {
            double *swap = NULL;


            communicateBoundaryValues( MPISendTypes, MPIRecvTypes, collideField, neighbours, MPISendReq, MPIRecvReq, 0); // communicate along x - axis (right to left - left to right)
#ifdef _AVX_
            doStreamingAndCollisionAVX(collideField, streamField, flagField, local_xlength, tau, 1);
#else
            doStreamingAndCollision(collideField, streamField, flagField, local_xlength, tau, 1);
#endif // _AVX_
            checkRequestCompletion(MPISendReq, MPIRecvReq, 0, neighbours);

            communicateBoundaryValues( MPISendTypes, MPIRecvTypes, collideField, neighbours, MPISendReq, MPIRecvReq, 2); // communicate along z - axis (back to front - front to back)
#ifdef _AVX_
            doStreamingAndCollisionAVX(collideField, streamField, flagField, local_xlength, tau, 2);
#else
            doStreamingAndCollision(collideField, streamField, flagField, local_xlength, tau, 2);
#endif // _AVX_
            checkRequestCompletion(MPISendReq, MPIRecvReq, 2, neighbours);

            communicateBoundaryValues( MPISendTypes, MPIRecvTypes, collideField, neighbours, MPISendReq, MPIRecvReq, 1); // communicate along y - axis (bottom to top - top to bottom)
#ifdef _AVX_
            doStreamingAndCollisionAVX(collideField, streamField, flagField, local_xlength, tau, 3);
#else
            doStreamingAndCollision(collideField, streamField, flagField, local_xlength, tau, 3);
#endif // _AVX_
            checkRequestCompletion(MPISendReq, MPIRecvReq, 1, neighbours);

#ifdef _AVX_
            doStreamingAndCollisionAVX(collideField, streamField, flagField, local_xlength, tau, 4);
#else
            doStreamingAndCollision(collideField, streamField, flagField, local_xlength, tau, 4);
#endif // _AVX_

            swap = collideField;
            collideField = streamField;
            streamField = swap;

            treatBoundary(collideField, flagField, bddParams, local_xlength);
#ifdef _VTK_
            if (t % timestepsPerPlotting == 0)
#ifdef _MPI_VTK_
                MPI_writeVtkOutput( &fh, collideField, flagField, local_xlength, iCoord, jCoord, kCoord, iProc, jProc, kProc, xlength);
#else
                writeVtkOutput(collideField, flagField, "./Paraview/output", (unsigned int) t / timestepsPerPlotting, xlength, local_xlength, rank, iCoord, jCoord, kCoord, iProc, jProc, kProc);
#endif // _MPI_VTK_
#endif // _VTK_

#ifndef _NOPROGRESS_
            if (!rank)
            {
                int pct = ((float) t / timesteps) * 100;
                printf("\b\b\b%02d%%", pct);
                fflush(stdout);
            }
#endif // _NOPROGRESS_

            /** debugging code: check collideField */
#ifdef _DEBUG_
            // check correctness of collideField with reference data
            if (t % timestepsPerPlotting == 0)
            {

                exactCollideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2));

                int x, y, z, i;
                FILE *fp = NULL;
                unsigned int line = 0;
                int error = 0;
                char szFileName[1200];
                sprintf( szFileName, "Testdata/%s/%i.dat", problem, t / timestepsPerPlotting );
                fp = fopen(szFileName,"r");
                if (fp != NULL)
                {
                    for (line = 0; line < PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2); line++)
                        fscanf(fp,"%lf",&exactCollideField[line]);

                    for (z = 1; z <= local_xlength[2]; z++)
                        for (y = 1; y <= local_xlength[1]; y++)
                            for(x = 1; x <= local_xlength[0]; x++)
                                for (i = 0; i < PARAMQ; i++)
                                    if (flagField[z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x] == FLUID)
                                        if (fabs(collideField[(z * (local_xlength[0] + 2) * (local_xlength[1] + 2) + y * (local_xlength[0] + 2) + x) + i * (local_xlength[0] + 2) * (local_xlength[1] + 2) * (local_xlength[2] + 2)] - exactCollideField[PARAMQ * ((z + kCoord * (xlength[2] / kProc)) * (xlength[0] + 2) * (xlength[1] + 2) + (y + jCoord * (xlength[1]/jProc)) * (xlength[0] + 2) + (x + iCoord * (xlength[0] / iProc) )) + i]) > 1e-4)
                                            error = 1;
                    if (error)
                        printf("ERROR: Process %d has a different collideField in timestep %d\n",rank, t);

                    fclose(fp);

                }
                else
                    printf("ERROR: Process %d cannot read file %s\n", rank, szFileName);

                free(exactCollideField);

            }
#endif
            /** end of debugging code */


        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (!rank)
        {
//            PAPI_read_counters( counters, 3 );
            printf("\b\b\b\b100%%\n");
            end = clock();
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
            gettimeofday(&time_end, NULL);

            printf("Running time (CPU time): %.2fs\n", time_spent);
            printf("Running time (Wall clock): %.2fs\n", ( (double) (( time_end.tv_sec - time_start.tv_sec) * 1000000u + time_end.tv_usec - time_start.tv_usec) )/ 1e6);
            printf("MLUPS: %.3f\n", ((double) (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2) * timesteps) / (1000000.0 * ((time_end.tv_sec - time_start.tv_sec) * 1000000u + time_end.tv_usec - time_start.tv_usec) / 1e6));

//            printf("%lld L2 cache misses (%.3lf%% misses) in %lld cycles\n",
//                   counters[1],
//                   (double)counters[1] / (double)counters[2] * 100,
//                   counters[0] );
        }

        free(collideField);
        free(streamField);
        free(flagField);
    }
    MPI_File_close(&fh);
    MPI_Finalize();
    return 0;
}

#endif
