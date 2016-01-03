#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include <time.h>

int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    int *flagField = NULL;
    clock_t begin, end;
    double time_spent;

    int xlength, timesteps, timestepsPerPlotting;
    double tau, velocityWall[3];


    if(readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv) == 0)
    {
        begin = clock();
        collideField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength + 2)*(xlength + 2)*(xlength + 2));
        streamField = (double*) malloc((size_t) sizeof(double) * PARAMQ * (xlength + 2)*(xlength + 2)*(xlength + 2));
        flagField = (int *) malloc((size_t) sizeof (int) * (xlength + 2)*(xlength + 2)*(xlength + 2));
        initialiseFields(collideField, streamField, flagField, xlength);

        writeVtkOutput(streamField, flagField, "./Paraview/output", 0, xlength);

        for(int t = 0; t < timesteps; t++)
        {
            double *swap = NULL;

            doStreaming(collideField, streamField, flagField, xlength);
            swap = collideField;
            collideField = streamField;
            streamField = swap;

            doCollision(collideField, flagField, &tau, xlength);

            treatBoundary(collideField, flagField, velocityWall, xlength);

            if (t % timestepsPerPlotting == 0)
                writeVtkOutput(collideField, flagField, "./Paraview/output", (unsigned int) t / timestepsPerPlotting, xlength);

        }
        end = clock();
        time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

        printf("Running time: %.2f\n", time_spent);
        printf("MLUPS: %.3f\n", ((double) (xlength + 2) * (xlength + 2) * (xlength + 2) * timesteps) / (1000000.0 * time_spent));

        free(collideField);
        free(streamField);
        free(flagField);

    }
    return 0;
}

#endif

