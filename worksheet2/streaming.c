#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength)
{
    int streamCellIndex, neighbourCellIndex;
    for (int z = 1; z <= xlength ; z++)
        for (int y = 1; y <= xlength ; y++)
            for (int x = 1; x <= xlength ; x++)
                for (int i = 0; i < PARAMQ; i++)
                {
                    streamCellIndex = PARAMQ * (x + y * (xlength + 2) + z * (xlength + 2) * (xlength + 2)) + i;
                    neighbourCellIndex = streamCellIndex - PARAMQ * (LATTICEVELOCITIES[i][0] + LATTICEVELOCITIES[i][1] * (xlength + 2) + LATTICEVELOCITIES[i][2] * (xlength + 2) * (xlength + 2));
                    streamField[streamCellIndex] = collideField[neighbourCellIndex];

                }

}

