#include "streaming.h"
#include "LBDefinitions.h"
#include "helper.h"
#ifdef _AVX_
#include <immintrin.h>
#endif // _AVX_
#include <omp.h>
#include <unistd.h>

void doStreamingAndCollision(double *collideField, double *streamField,int *flagField,int *xlength, const double tau, int part)
{
    int x ,y ,z, j;
    double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18;
    double density, velocity[3];
    double tmpFeq0, tmpFeq1, tmpFeq2, dotProduct;
    #ifdef _ARBITRARYGEOMETRY_
    int streamMask, copyMask;
    #endif // _ARBITRARYGEOMETRY_

    // Philip Neumann approach
    int numberOfCells = (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);

    // During communication the innermost parts of the domain can already be processed. After sending in direction of the x-axis
    // we process a third of the domain (part == 1). After that we send in direction of the z-axis and then continue with the second
    // third of the domain (part == 2). Finally we send in direction of the y-axis and process the last third of the domain (part == 3).
    // Once we recieved all data from neighbours in y-direction we can stream from the local boundary cells as well (part == 4).

    switch(part)
    {
    case 1:
        #ifdef _ARBITRARYGEOMETRY__
        #pragma omp parallel for private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField) schedule(static)
        #else
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField) schedule(static)
        #endif
        for (z = 2; z <= xlength[2] / 3; z++)
        {
            j = 2 + 2 * (xlength[0] + 2) + (xlength[0] + 2) * (xlength[1] + 2) * z;
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 2; x <= xlength[0] - 1; x++)
                {
                    int counter = x + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2);
                    if (counter != j)
                        printf("ERROR: Counters differ!\n");
                    density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                    velocity[0] = velocity[0] / density;
                    velocity[1] = velocity[1] / density;
                    velocity[2] = velocity[2] / density;

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProduct = - velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                    streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j++;
                }
                j += 4; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += 4 * xlength[0] + 8;
        }
        break;
    case 2:
        #ifdef _ARBITRARYGEOMETRY__
        #pragma omp parallel for private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField) schedule(static)
        #else
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField) schedule(static)
        #endif
        for (z = xlength[2] / 3 + 1; z <= 2 * xlength[2] / 3; z++)
        {
            j = 1 + 2 * (xlength[0] + 2) + (xlength[0] + 2) * (xlength[1] + 2) * z;
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 1; x <= xlength[0]; x++)
                {
                    int counter = x + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2);
                    if (counter != j)
                        printf("ERROR: Counters differ!\n");

                   density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                    velocity[0] = velocity[0] / density;
                    velocity[1] = velocity[1] / density;
                    velocity[2] = velocity[2] / density;

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProduct = - velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                    streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j++;
                }
                j += 2; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += 4 * xlength[0] + 8;
        }
        break;
    case 3:
        #ifdef _ARBITRARYGEOMETRY__
        #pragma omp parallel for private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField) schedule(static)
        #else
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField) schedule(static)
        #endif
        for (z = 2 * xlength[2] / 3 + 1; z <= xlength[2]; z++)
        {
            j = 1 + 2 * (xlength[0] + 2) + (xlength[0] + 2) * (xlength[1] + 2) * z;
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 1; x <= xlength[0]; x++)
                {
                    int counter = x + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2);
                    if (counter != j)
                        printf("ERROR: Counters differ!\n");
                    density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                    velocity[0] = velocity[0] / density;
                    velocity[1] = velocity[1] / density;
                    velocity[2] = velocity[2] / density;

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProduct = - velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                    streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j++;
                }
                j += 2; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += 4 * xlength[0] + 8;
        }
        break;
    case 4:
        #ifdef _ARBITRARYGEOMETRY__
        #pragma omp parallel private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField)
        #else
        #pragma omp parallel private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, j), firstprivate(numberOfCells, tau), shared(xlength, streamField, collideField)
        #endif
        {
        /** z = 1 */
        #pragma omp for schedule(static)
        for (y = 1; y <= xlength[1]; y++)
        {
            j = 1 + (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2);
            for (x = 1; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;
                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2; // skip the last two cells of the current row and the first two cell of the next row
        }

        /** 1 < z <= xlength[2] / 3 **/
        #pragma omp for schedule(static)
        for (z = 2; z <= xlength[2] / 3; z++)
        {
            j = z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[0] + 2) + 1;
            // y = 1
            for (x = 1; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2;
            // 1 < y < xlength[1]
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                // x = 1
                int counter = 1 + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2);
                if (counter != j)
                    printf("ERROR: Counters differ!\n");
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;
                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;

                j += xlength[0] - 1;

                // x = xlength[0]
                counter = xlength[0] + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2);
                if (counter != j)
                    printf("ERROR: Counters differ!\n");
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += 3;
            }
            // y = xlength[1]
            for (x = 1; x <= xlength[0]; x++)
            {
                int counter = x + (xlength[0] + 2) * (xlength[1]) + z * (xlength[0] + 2) * (xlength[1] + 2);
                if (counter != j)
                    printf("ERROR: Counters differ!\n");
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;
                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2 * xlength[0] + 6;
        }

        /** xlength[2] / 3 < z <= xlength[2] */
        #pragma omp for schedule(static)
        for (z = xlength[2] / 3 + 1; z <= xlength[2]; z++)
        {
            // y = 1
            j = (xlength[0] + 2) * (xlength[1] + 2) * z + (xlength[0] + 2) + 1;
            for (x = 1; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }

            j += (xlength[0] + 2) * (xlength[1] - 1) - (xlength[0]);
            // y = xlength[1]
            for (x = 1; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                    #ifdef _ARBITRARYGEOMETRY_
                    streamMask = !flagField[j];
                    copyMask = !streamMask;
                    #endif // _ARBITRARYGEOMETRY_
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j];
                    #else
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells];
                    #else
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 2];
                    #else
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 3];
                    #else
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 4];
                    #else
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    #endif

                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3] * streamMask + copyMask * collideField[j + numberOfCells * 5];
                    #else
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    #endif

                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2] * streamMask + copyMask * collideField[j + numberOfCells * 6];
                    #else
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    #endif

                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1] * streamMask + copyMask * collideField[j + numberOfCells * 7];
                    #else
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    #endif

                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f8 = collideField[8 * numberOfCells + j + 1] * streamMask + copyMask * collideField[j + numberOfCells * 8];
                    #else
                    f8 = collideField[8 * numberOfCells + j + 1];
                    #endif

                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f10 = collideField[10 * numberOfCells + j - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f10 = collideField[10 * numberOfCells + j - 1];
                    #endif

                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1] * streamMask + copyMask * collideField[j + numberOfCells * 11];
                    #else
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    #endif

                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2] * streamMask + copyMask * collideField[j + numberOfCells * 12];
                    #else
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    #endif

                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    #ifdef _ARBITRARYGEOMETRY_
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3] * streamMask + copyMask * collideField[j + numberOfCells * 13];
                    #else
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    #endif

                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)] * streamMask + copyMask * collideField[j + numberOfCells * 14];
                    #else
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    #endif

                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1] * streamMask + copyMask * collideField[j + numberOfCells * 15];
                    #else
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    #endif

                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)] * streamMask + copyMask * collideField[j + numberOfCells * 16];
                    #else
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    #endif

                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1] * streamMask + copyMask * collideField[j + numberOfCells * 17];
                    #else
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    #endif

                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    #ifdef _ARBITRARYGEOMETRY_
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + copyMask * collideField[j + numberOfCells * 18];
                    #else
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    #endif

                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;
                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2 * xlength[0] + 6;
        }
    }
    break;
    }
}
#ifdef _AVX_
void doStreamingAndCollisionAVX(double *collideField, double *streamField,int *flagField,int *xlength, const double tau, int part)
{

    int x ,y ,z, j;
    double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18;
    double density, velocity[3];
    double tmpFeq0, tmpFeq1, tmpFeq2, dotProduct;


    __m256d f0v, f1v, f2v, f3v, f4v, f5v, f6v, f7v, f8v, f9v, f10v, f11v, f12v, f13v, f14v, f15v, f16v, f17v, f18v;
    __m256d densityVector, velocityX, velocityY, velocityZ;
    __m256d tmpFeq0v, tmpFeq1v, tmpFeq2v, dotProductVector;
    const double zero = 0.0;
    const double one = 1.0;
    const double c_s_square = C_S * C_S;
    __m256d cs2v = _mm256_broadcast_sd(&c_s_square);
    __m256d oneVector = _mm256_broadcast_sd(&one);
    __m256d zeroVector = _mm256_broadcast_sd(&zero);
    __m256d tauInvVector = _mm256_div_pd(oneVector, _mm256_broadcast_sd(&tau));
    __m256d weight36 = _mm256_broadcast_sd(&LATTICEWEIGHTS[7]);
    __m256d weight18 = _mm256_broadcast_sd(&LATTICEWEIGHTS[8]);
    __m256d weight3 = _mm256_broadcast_sd(&LATTICEWEIGHTS[9]);



    // Neumann approach
    int numberOfCells = (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);

    switch(part)
    {
    case 1:
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f0v, f1v, f2v, f3v, f4v, f5v, f6v, f7v, f8v, f9v, f10v, f11v, f12v, f13v, f14v, f15v, f16v, f17v, f18v, density, velocity, densityVector, velocityX, velocityY, velocityZ, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, tmpFeq0v, tmpFeq1v, tmpFeq2v, dotProductVector, x, y, j), firstprivate(numberOfCells, tau, tauInvVector, weight3, weight18, weight36, zeroVector, oneVector, cs2v), shared(xlength, streamField, collideField) schedule(static)
        for (z = 2; z <= xlength[2] / 3; z++)
        {
            j = 2 + 2 * (xlength[0] + 2) + (xlength[0] + 2) * (xlength[1] + 2) * z;
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 2; x <= xlength[0] - 5; x += 4)
                {
                    densityVector = velocityX = velocityY = velocityZ = zeroVector;
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                    densityVector = _mm256_add_pd(densityVector, f0v);
                    velocityY = _mm256_sub_pd(velocityY, f0v);
                    velocityZ = _mm256_sub_pd(velocityZ, f0v);

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                    densityVector = _mm256_add_pd(densityVector, f1v);
                    velocityX = _mm256_sub_pd(velocityX, f1v);
                    velocityZ = _mm256_sub_pd(velocityZ, f1v);

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                    densityVector = _mm256_add_pd(densityVector, f2v);
                    velocityZ = _mm256_sub_pd(velocityZ, f2v);

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                    densityVector = _mm256_add_pd(densityVector, f3v);
                    velocityX = _mm256_add_pd(velocityX, f3v);
                    velocityZ = _mm256_sub_pd(velocityZ, f3v);

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                    densityVector = _mm256_add_pd(densityVector, f4v);
                    velocityY = _mm256_add_pd(velocityY, f4v);
                    velocityZ = _mm256_sub_pd(velocityZ, f4v);

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                    densityVector = _mm256_add_pd(densityVector, f5v);
                    velocityX = _mm256_sub_pd(velocityX, f5v);
                    velocityY = _mm256_sub_pd(velocityY, f5v);

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                    densityVector = _mm256_add_pd(densityVector, f6v);
                    velocityY = _mm256_sub_pd(velocityY, f6v);


                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                    densityVector = _mm256_add_pd(densityVector, f7v);
                    velocityX = _mm256_add_pd(velocityX, f7v);
                    velocityY = _mm256_sub_pd(velocityY, f7v);

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                    densityVector = _mm256_add_pd(densityVector, f8v);
                    velocityX = _mm256_sub_pd(velocityX, f8v);

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                    densityVector = _mm256_add_pd(densityVector, f9v);

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                    densityVector = _mm256_add_pd(densityVector, f10v);
                    velocityX = _mm256_add_pd(velocityX, f10v);

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                    densityVector = _mm256_add_pd(densityVector, f11v);
                    velocityX = _mm256_sub_pd(velocityX, f11v);
                    velocityY = _mm256_add_pd(velocityY, f11v);

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                    densityVector = _mm256_add_pd(densityVector, f12v);
                    velocityY = _mm256_add_pd(velocityY, f12v);

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                    densityVector = _mm256_add_pd(densityVector, f13v);
                    velocityX = _mm256_add_pd(velocityX, f13v);
                    velocityY = _mm256_add_pd(velocityY, f13v);

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                    densityVector = _mm256_add_pd(densityVector, f14v);
                    velocityY = _mm256_sub_pd(velocityY, f14v);
                    velocityZ = _mm256_add_pd(velocityZ, f14v);

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                    densityVector = _mm256_add_pd(densityVector, f15v);
                    velocityX = _mm256_sub_pd(velocityX, f15v);
                    velocityZ = _mm256_add_pd(velocityZ, f15v);

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                    densityVector = _mm256_add_pd(densityVector, f16v);
                    velocityZ = _mm256_add_pd(velocityZ, f16v);

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                    densityVector = _mm256_add_pd(densityVector, f17v);
                    velocityX = _mm256_add_pd(velocityX, f17v);
                    velocityZ = _mm256_add_pd(velocityZ, f17v);

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                    densityVector = _mm256_add_pd(densityVector, f18v);
                    velocityY = _mm256_add_pd(velocityY, f18v);
                    velocityZ = _mm256_add_pd(velocityZ, f18v);

                    velocityX = _mm256_div_pd(velocityX, densityVector);
                    velocityY = _mm256_div_pd(velocityY, densityVector);
                    velocityZ = _mm256_div_pd(velocityZ, densityVector);

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                    j += 4;
                }

                /** sequential rest */
                for (; x <= xlength[0] - 1; x++)
                {
                    density = velocity[0] = velocity[1] = velocity[2] = 0;

                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    f8 = collideField[8 * numberOfCells + j + 1];
                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    f10 = collideField[10 * numberOfCells + j - 1];
                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                    velocity[0] = velocity[0] / density;
                    velocity[1] = velocity[1] / density;
                    velocity[2] = velocity[2] / density;

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProduct = - velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                    streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j++;
                }
                j += 4; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += 4 * xlength[0] + 8;
        }
        break;
    case 2:
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f0v, f1v, f2v, f3v, f4v, f5v, f6v, f7v, f8v, f9v, f10v, f11v, f12v, f13v, f14v, f15v, f16v, f17v, f18v, density, velocity, densityVector, velocityX, velocityY, velocityZ, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, tmpFeq0v, tmpFeq1v, tmpFeq2v, dotProductVector, x, y, j), firstprivate(numberOfCells, tau, tauInvVector, weight3, weight18, weight36, zeroVector, oneVector, cs2v), shared(xlength, streamField, collideField) schedule(static)
        for (z = xlength[2] / 3 + 1; z <= 2 * xlength[2] / 3; z++)
        {
            j = 1 + 2 * (xlength[0] + 2) + (xlength[0] + 2) * (xlength[1] + 2) * z;
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 1; x <= xlength[0] - 4; x += 4)
                {
                    densityVector = velocityX = velocityY = velocityZ = zeroVector;
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                    densityVector = _mm256_add_pd(densityVector, f0v);
                    velocityY = _mm256_sub_pd(velocityY, f0v);
                    velocityZ = _mm256_sub_pd(velocityZ, f0v);

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                    densityVector = _mm256_add_pd(densityVector, f1v);
                    velocityX = _mm256_sub_pd(velocityX, f1v);
                    velocityZ = _mm256_sub_pd(velocityZ, f1v);

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                    densityVector = _mm256_add_pd(densityVector, f2v);
                    velocityZ = _mm256_sub_pd(velocityZ, f2v);

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                    densityVector = _mm256_add_pd(densityVector, f3v);
                    velocityX = _mm256_add_pd(velocityX, f3v);
                    velocityZ = _mm256_sub_pd(velocityZ, f3v);

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                    densityVector = _mm256_add_pd(densityVector, f4v);
                    velocityY = _mm256_add_pd(velocityY, f4v);
                    velocityZ = _mm256_sub_pd(velocityZ, f4v);

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                    densityVector = _mm256_add_pd(densityVector, f5v);
                    velocityX = _mm256_sub_pd(velocityX, f5v);
                    velocityY = _mm256_sub_pd(velocityY, f5v);

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                    densityVector = _mm256_add_pd(densityVector, f6v);
                    velocityY = _mm256_sub_pd(velocityY, f6v);


                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                    densityVector = _mm256_add_pd(densityVector, f7v);
                    velocityX = _mm256_add_pd(velocityX, f7v);
                    velocityY = _mm256_sub_pd(velocityY, f7v);

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                    densityVector = _mm256_add_pd(densityVector, f8v);
                    velocityX = _mm256_sub_pd(velocityX, f8v);

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                    densityVector = _mm256_add_pd(densityVector, f9v);

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                    densityVector = _mm256_add_pd(densityVector, f10v);
                    velocityX = _mm256_add_pd(velocityX, f10v);

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                    densityVector = _mm256_add_pd(densityVector, f11v);
                    velocityX = _mm256_sub_pd(velocityX, f11v);
                    velocityY = _mm256_add_pd(velocityY, f11v);

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                    densityVector = _mm256_add_pd(densityVector, f12v);
                    velocityY = _mm256_add_pd(velocityY, f12v);

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                    densityVector = _mm256_add_pd(densityVector, f13v);
                    velocityX = _mm256_add_pd(velocityX, f13v);
                    velocityY = _mm256_add_pd(velocityY, f13v);

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                    densityVector = _mm256_add_pd(densityVector, f14v);
                    velocityY = _mm256_sub_pd(velocityY, f14v);
                    velocityZ = _mm256_add_pd(velocityZ, f14v);

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                    densityVector = _mm256_add_pd(densityVector, f15v);
                    velocityX = _mm256_sub_pd(velocityX, f15v);
                    velocityZ = _mm256_add_pd(velocityZ, f15v);

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                    densityVector = _mm256_add_pd(densityVector, f16v);
                    velocityZ = _mm256_add_pd(velocityZ, f16v);

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                    densityVector = _mm256_add_pd(densityVector, f17v);
                    velocityX = _mm256_add_pd(velocityX, f17v);
                    velocityZ = _mm256_add_pd(velocityZ, f17v);

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                    densityVector = _mm256_add_pd(densityVector, f18v);
                    velocityY = _mm256_add_pd(velocityY, f18v);
                    velocityZ = _mm256_add_pd(velocityZ, f18v);

                    velocityX = _mm256_div_pd(velocityX, densityVector);
                    velocityY = _mm256_div_pd(velocityY, densityVector);
                    velocityZ = _mm256_div_pd(velocityZ, densityVector);

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                    j += 4;
                }

                /** sequential rest */
                for (; x <= xlength[0]; x++)
                {
                    density = velocity[0] = velocity[1] = velocity[2] = 0;

                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    f8 = collideField[8 * numberOfCells + j + 1];
                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    f10 = collideField[10 * numberOfCells + j - 1];
                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                    velocity[0] = velocity[0] / density;
                    velocity[1] = velocity[1] / density;
                    velocity[2] = velocity[2] / density;

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProduct = - velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                    streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j++;
                }
                j += 2; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += 4 * xlength[0] + 8;
        }
        break;
    case 3:
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f0v, f1v, f2v, f3v, f4v, f5v, f6v, f7v, f8v, f9v, f10v, f11v, f12v, f13v, f14v, f15v, f16v, f17v, f18v, density, velocity, densityVector, velocityX, velocityY, velocityZ, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, tmpFeq0v, tmpFeq1v, tmpFeq2v, dotProductVector, x, y, j), firstprivate(numberOfCells, tau, tauInvVector, weight3, weight18, weight36, zeroVector, oneVector, cs2v), shared(xlength, streamField, collideField) schedule(static)
        for (z = 2 * xlength[2] / 3 + 1; z <= xlength[2] - 1; z++)
        {
            j = xlength[0] + 3 + (xlength[0] + 2) * (xlength[1] + 2) * z;
            for (y = 1; y <= xlength[1]; y++)
            {
                for (x = 1; x <= xlength[0] - 4; x += 4)
                {
                    densityVector = velocityX = velocityY = velocityZ = zeroVector;
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                    densityVector = _mm256_add_pd(densityVector, f0v);
                    velocityY = _mm256_sub_pd(velocityY, f0v);
                    velocityZ = _mm256_sub_pd(velocityZ, f0v);

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                    densityVector = _mm256_add_pd(densityVector, f1v);
                    velocityX = _mm256_sub_pd(velocityX, f1v);
                    velocityZ = _mm256_sub_pd(velocityZ, f1v);

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                    densityVector = _mm256_add_pd(densityVector, f2v);
                    velocityZ = _mm256_sub_pd(velocityZ, f2v);

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                    densityVector = _mm256_add_pd(densityVector, f3v);
                    velocityX = _mm256_add_pd(velocityX, f3v);
                    velocityZ = _mm256_sub_pd(velocityZ, f3v);

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                    densityVector = _mm256_add_pd(densityVector, f4v);
                    velocityY = _mm256_add_pd(velocityY, f4v);
                    velocityZ = _mm256_sub_pd(velocityZ, f4v);

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                    densityVector = _mm256_add_pd(densityVector, f5v);
                    velocityX = _mm256_sub_pd(velocityX, f5v);
                    velocityY = _mm256_sub_pd(velocityY, f5v);

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                    densityVector = _mm256_add_pd(densityVector, f6v);
                    velocityY = _mm256_sub_pd(velocityY, f6v);


                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                    densityVector = _mm256_add_pd(densityVector, f7v);
                    velocityX = _mm256_add_pd(velocityX, f7v);
                    velocityY = _mm256_sub_pd(velocityY, f7v);

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                    densityVector = _mm256_add_pd(densityVector, f8v);
                    velocityX = _mm256_sub_pd(velocityX, f8v);

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                    densityVector = _mm256_add_pd(densityVector, f9v);

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                    densityVector = _mm256_add_pd(densityVector, f10v);
                    velocityX = _mm256_add_pd(velocityX, f10v);

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                    densityVector = _mm256_add_pd(densityVector, f11v);
                    velocityX = _mm256_sub_pd(velocityX, f11v);
                    velocityY = _mm256_add_pd(velocityY, f11v);

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                    densityVector = _mm256_add_pd(densityVector, f12v);
                    velocityY = _mm256_add_pd(velocityY, f12v);

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                    densityVector = _mm256_add_pd(densityVector, f13v);
                    velocityX = _mm256_add_pd(velocityX, f13v);
                    velocityY = _mm256_add_pd(velocityY, f13v);

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                    densityVector = _mm256_add_pd(densityVector, f14v);
                    velocityY = _mm256_sub_pd(velocityY, f14v);
                    velocityZ = _mm256_add_pd(velocityZ, f14v);

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                    densityVector = _mm256_add_pd(densityVector, f15v);
                    velocityX = _mm256_sub_pd(velocityX, f15v);
                    velocityZ = _mm256_add_pd(velocityZ, f15v);

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                    densityVector = _mm256_add_pd(densityVector, f16v);
                    velocityZ = _mm256_add_pd(velocityZ, f16v);

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                    densityVector = _mm256_add_pd(densityVector, f17v);
                    velocityX = _mm256_add_pd(velocityX, f17v);
                    velocityZ = _mm256_add_pd(velocityZ, f17v);

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                    densityVector = _mm256_add_pd(densityVector, f18v);
                    velocityY = _mm256_add_pd(velocityY, f18v);
                    velocityZ = _mm256_add_pd(velocityZ, f18v);

                    velocityX = _mm256_div_pd(velocityX, densityVector);
                    velocityY = _mm256_div_pd(velocityY, densityVector);
                    velocityZ = _mm256_div_pd(velocityZ, densityVector);

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                    _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                    _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                    j += 4;
                }

                /** sequential rest */
                for (; x <= xlength[0]; x++)
                {
                    density = velocity[0] = velocity[1] = velocity[2] = 0;

                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
                    f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
                    f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
                    f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
                    f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
                    f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
                    f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
                    f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
                    f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
                    f8 = collideField[8 * numberOfCells + j + 1];
                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[9 * numberOfCells + j];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
                    f10 = collideField[10 * numberOfCells + j - 1];
                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
                    f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
                    f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
                    f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
                    f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
                    f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
                    f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
                    f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
                    f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                    density += f18;
                    velocity[1] += f18;
                    velocity[2] += f18;

                    velocity[0] = velocity[0] / density;
                    velocity[1] = velocity[1] / density;
                    velocity[2] = velocity[2] / density;

                    // compute feq and relaxate; then store to streamField
                    tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                    // i = 0 and i = 18
                    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                    dotProduct = - velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                    streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j++;
                }
                j += 2; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += 2 * xlength[0] + 4;
        }
        break;
    case 4:
        #pragma omp parallel private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f0v, f1v, f2v, f3v, f4v, f5v, f6v, f7v, f8v, f9v, f10v, f11v, f12v, f13v, f14v, f15v, f16v, f17v, f18v, density, velocity, densityVector, velocityX, velocityY, velocityZ, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, tmpFeq0v, tmpFeq1v, tmpFeq2v, dotProductVector, x, y, z, j), firstprivate(numberOfCells, tau, tauInvVector, weight3, weight18, weight36, zeroVector, oneVector, cs2v), shared(xlength, streamField, collideField)
    {
        /** z = 1 */
        #pragma omp for schedule(static)
        for (y = 1; y <= xlength[1]; y++)
        {
            j = 1 + (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2);
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                densityVector = velocityX = velocityY = velocityZ = zeroVector;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                densityVector = _mm256_add_pd(densityVector, f0v);
                velocityY = _mm256_sub_pd(velocityY, f0v);
                velocityZ = _mm256_sub_pd(velocityZ, f0v);

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                densityVector = _mm256_add_pd(densityVector, f1v);
                velocityX = _mm256_sub_pd(velocityX, f1v);
                velocityZ = _mm256_sub_pd(velocityZ, f1v);

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                densityVector = _mm256_add_pd(densityVector, f2v);
                velocityZ = _mm256_sub_pd(velocityZ, f2v);

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                densityVector = _mm256_add_pd(densityVector, f3v);
                velocityX = _mm256_add_pd(velocityX, f3v);
                velocityZ = _mm256_sub_pd(velocityZ, f3v);

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                densityVector = _mm256_add_pd(densityVector, f4v);
                velocityY = _mm256_add_pd(velocityY, f4v);
                velocityZ = _mm256_sub_pd(velocityZ, f4v);

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                densityVector = _mm256_add_pd(densityVector, f5v);
                velocityX = _mm256_sub_pd(velocityX, f5v);
                velocityY = _mm256_sub_pd(velocityY, f5v);

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                densityVector = _mm256_add_pd(densityVector, f6v);
                velocityY = _mm256_sub_pd(velocityY, f6v);


                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                densityVector = _mm256_add_pd(densityVector, f7v);
                velocityX = _mm256_add_pd(velocityX, f7v);
                velocityY = _mm256_sub_pd(velocityY, f7v);

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                densityVector = _mm256_add_pd(densityVector, f8v);
                velocityX = _mm256_sub_pd(velocityX, f8v);

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                densityVector = _mm256_add_pd(densityVector, f9v);

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                densityVector = _mm256_add_pd(densityVector, f10v);
                velocityX = _mm256_add_pd(velocityX, f10v);

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                densityVector = _mm256_add_pd(densityVector, f11v);
                velocityX = _mm256_sub_pd(velocityX, f11v);
                velocityY = _mm256_add_pd(velocityY, f11v);

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                densityVector = _mm256_add_pd(densityVector, f12v);
                velocityY = _mm256_add_pd(velocityY, f12v);

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                densityVector = _mm256_add_pd(densityVector, f13v);
                velocityX = _mm256_add_pd(velocityX, f13v);
                velocityY = _mm256_add_pd(velocityY, f13v);

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                densityVector = _mm256_add_pd(densityVector, f14v);
                velocityY = _mm256_sub_pd(velocityY, f14v);
                velocityZ = _mm256_add_pd(velocityZ, f14v);

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                densityVector = _mm256_add_pd(densityVector, f15v);
                velocityX = _mm256_sub_pd(velocityX, f15v);
                velocityZ = _mm256_add_pd(velocityZ, f15v);

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                densityVector = _mm256_add_pd(densityVector, f16v);
                velocityZ = _mm256_add_pd(velocityZ, f16v);

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                densityVector = _mm256_add_pd(densityVector, f17v);
                velocityX = _mm256_add_pd(velocityX, f17v);
                velocityZ = _mm256_add_pd(velocityZ, f17v);

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                densityVector = _mm256_add_pd(densityVector, f18v);
                velocityY = _mm256_add_pd(velocityY, f18v);
                velocityZ = _mm256_add_pd(velocityZ, f18v);

                velocityX = _mm256_div_pd(velocityX, densityVector);
                velocityY = _mm256_div_pd(velocityY, densityVector);
                velocityZ = _mm256_div_pd(velocityZ, densityVector);

                // compute feq and relaxate; then store to streamField
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                j += 4;
            }

            /** sequential rest */
            for (; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;

                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2; // skip the last two cells of the current row and the first two cell of the next row
        }

        /** 1 < z <= xlength[2] / 3 **/
        #pragma omp for schedule(static)
        for (z = 2; z <= xlength[2] / 3; z++)
        {
            j = z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[0] + 2) + 1;
            // y = 1
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                densityVector = velocityX = velocityY = velocityZ = zeroVector;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                densityVector = _mm256_add_pd(densityVector, f0v);
                velocityY = _mm256_sub_pd(velocityY, f0v);
                velocityZ = _mm256_sub_pd(velocityZ, f0v);

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                densityVector = _mm256_add_pd(densityVector, f1v);
                velocityX = _mm256_sub_pd(velocityX, f1v);
                velocityZ = _mm256_sub_pd(velocityZ, f1v);

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                densityVector = _mm256_add_pd(densityVector, f2v);
                velocityZ = _mm256_sub_pd(velocityZ, f2v);

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                densityVector = _mm256_add_pd(densityVector, f3v);
                velocityX = _mm256_add_pd(velocityX, f3v);
                velocityZ = _mm256_sub_pd(velocityZ, f3v);

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                densityVector = _mm256_add_pd(densityVector, f4v);
                velocityY = _mm256_add_pd(velocityY, f4v);
                velocityZ = _mm256_sub_pd(velocityZ, f4v);

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                densityVector = _mm256_add_pd(densityVector, f5v);
                velocityX = _mm256_sub_pd(velocityX, f5v);
                velocityY = _mm256_sub_pd(velocityY, f5v);

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                densityVector = _mm256_add_pd(densityVector, f6v);
                velocityY = _mm256_sub_pd(velocityY, f6v);


                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                densityVector = _mm256_add_pd(densityVector, f7v);
                velocityX = _mm256_add_pd(velocityX, f7v);
                velocityY = _mm256_sub_pd(velocityY, f7v);

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                densityVector = _mm256_add_pd(densityVector, f8v);
                velocityX = _mm256_sub_pd(velocityX, f8v);

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                densityVector = _mm256_add_pd(densityVector, f9v);

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                densityVector = _mm256_add_pd(densityVector, f10v);
                velocityX = _mm256_add_pd(velocityX, f10v);

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                densityVector = _mm256_add_pd(densityVector, f11v);
                velocityX = _mm256_sub_pd(velocityX, f11v);
                velocityY = _mm256_add_pd(velocityY, f11v);

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                densityVector = _mm256_add_pd(densityVector, f12v);
                velocityY = _mm256_add_pd(velocityY, f12v);

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                densityVector = _mm256_add_pd(densityVector, f13v);
                velocityX = _mm256_add_pd(velocityX, f13v);
                velocityY = _mm256_add_pd(velocityY, f13v);

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                densityVector = _mm256_add_pd(densityVector, f14v);
                velocityY = _mm256_sub_pd(velocityY, f14v);
                velocityZ = _mm256_add_pd(velocityZ, f14v);

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                densityVector = _mm256_add_pd(densityVector, f15v);
                velocityX = _mm256_sub_pd(velocityX, f15v);
                velocityZ = _mm256_add_pd(velocityZ, f15v);

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                densityVector = _mm256_add_pd(densityVector, f16v);
                velocityZ = _mm256_add_pd(velocityZ, f16v);

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                densityVector = _mm256_add_pd(densityVector, f17v);
                velocityX = _mm256_add_pd(velocityX, f17v);
                velocityZ = _mm256_add_pd(velocityZ, f17v);

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                densityVector = _mm256_add_pd(densityVector, f18v);
                velocityY = _mm256_add_pd(velocityY, f18v);
                velocityZ = _mm256_add_pd(velocityZ, f18v);

                velocityX = _mm256_div_pd(velocityX, densityVector);
                velocityY = _mm256_div_pd(velocityY, densityVector);
                velocityZ = _mm256_div_pd(velocityZ, densityVector);

                // compute feq and relaxate; then store to streamField
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                j += 4;
            }

            /** sequential rest */
            for (; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;

                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2;
            // 1 < y < xlength[1]
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                // x = 1
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;

                j += xlength[0] - 1;

                // x = xlength[0]
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += 3;
            }
            // y = xlength[1]
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                densityVector = velocityX = velocityY = velocityZ = zeroVector;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                densityVector = _mm256_add_pd(densityVector, f0v);
                velocityY = _mm256_sub_pd(velocityY, f0v);
                velocityZ = _mm256_sub_pd(velocityZ, f0v);

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                densityVector = _mm256_add_pd(densityVector, f1v);
                velocityX = _mm256_sub_pd(velocityX, f1v);
                velocityZ = _mm256_sub_pd(velocityZ, f1v);

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                densityVector = _mm256_add_pd(densityVector, f2v);
                velocityZ = _mm256_sub_pd(velocityZ, f2v);

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                densityVector = _mm256_add_pd(densityVector, f3v);
                velocityX = _mm256_add_pd(velocityX, f3v);
                velocityZ = _mm256_sub_pd(velocityZ, f3v);

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                densityVector = _mm256_add_pd(densityVector, f4v);
                velocityY = _mm256_add_pd(velocityY, f4v);
                velocityZ = _mm256_sub_pd(velocityZ, f4v);

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                densityVector = _mm256_add_pd(densityVector, f5v);
                velocityX = _mm256_sub_pd(velocityX, f5v);
                velocityY = _mm256_sub_pd(velocityY, f5v);

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                densityVector = _mm256_add_pd(densityVector, f6v);
                velocityY = _mm256_sub_pd(velocityY, f6v);


                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                densityVector = _mm256_add_pd(densityVector, f7v);
                velocityX = _mm256_add_pd(velocityX, f7v);
                velocityY = _mm256_sub_pd(velocityY, f7v);

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                densityVector = _mm256_add_pd(densityVector, f8v);
                velocityX = _mm256_sub_pd(velocityX, f8v);

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                densityVector = _mm256_add_pd(densityVector, f9v);

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                densityVector = _mm256_add_pd(densityVector, f10v);
                velocityX = _mm256_add_pd(velocityX, f10v);

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                densityVector = _mm256_add_pd(densityVector, f11v);
                velocityX = _mm256_sub_pd(velocityX, f11v);
                velocityY = _mm256_add_pd(velocityY, f11v);

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                densityVector = _mm256_add_pd(densityVector, f12v);
                velocityY = _mm256_add_pd(velocityY, f12v);

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                densityVector = _mm256_add_pd(densityVector, f13v);
                velocityX = _mm256_add_pd(velocityX, f13v);
                velocityY = _mm256_add_pd(velocityY, f13v);

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                densityVector = _mm256_add_pd(densityVector, f14v);
                velocityY = _mm256_sub_pd(velocityY, f14v);
                velocityZ = _mm256_add_pd(velocityZ, f14v);

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                densityVector = _mm256_add_pd(densityVector, f15v);
                velocityX = _mm256_sub_pd(velocityX, f15v);
                velocityZ = _mm256_add_pd(velocityZ, f15v);

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                densityVector = _mm256_add_pd(densityVector, f16v);
                velocityZ = _mm256_add_pd(velocityZ, f16v);

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                densityVector = _mm256_add_pd(densityVector, f17v);
                velocityX = _mm256_add_pd(velocityX, f17v);
                velocityZ = _mm256_add_pd(velocityZ, f17v);

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                densityVector = _mm256_add_pd(densityVector, f18v);
                velocityY = _mm256_add_pd(velocityY, f18v);
                velocityZ = _mm256_add_pd(velocityZ, f18v);

                velocityX = _mm256_div_pd(velocityX, densityVector);
                velocityY = _mm256_div_pd(velocityY, densityVector);
                velocityZ = _mm256_div_pd(velocityZ, densityVector);

                // compute feq and relaxate; then store to streamField
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                j += 4;
            }

            /** sequential rest */
            for (; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;

                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2 * xlength[0] + 6;
        }

        /** xlength[2] / 3 < z <= 2 * xlength[2] / 3 */
        #pragma omp for schedule(static)
        for (z = xlength[2] / 3 + 1; z <= 2 * xlength[2] / 3; z++)
        {
            // y = 1
            j = (xlength[0] + 2) * (xlength[1] + 2) * z + (xlength[0] + 2) + 1;
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                densityVector = velocityX = velocityY = velocityZ = zeroVector;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                densityVector = _mm256_add_pd(densityVector, f0v);
                velocityY = _mm256_sub_pd(velocityY, f0v);
                velocityZ = _mm256_sub_pd(velocityZ, f0v);

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                densityVector = _mm256_add_pd(densityVector, f1v);
                velocityX = _mm256_sub_pd(velocityX, f1v);
                velocityZ = _mm256_sub_pd(velocityZ, f1v);

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                densityVector = _mm256_add_pd(densityVector, f2v);
                velocityZ = _mm256_sub_pd(velocityZ, f2v);

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                densityVector = _mm256_add_pd(densityVector, f3v);
                velocityX = _mm256_add_pd(velocityX, f3v);
                velocityZ = _mm256_sub_pd(velocityZ, f3v);

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                densityVector = _mm256_add_pd(densityVector, f4v);
                velocityY = _mm256_add_pd(velocityY, f4v);
                velocityZ = _mm256_sub_pd(velocityZ, f4v);

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                densityVector = _mm256_add_pd(densityVector, f5v);
                velocityX = _mm256_sub_pd(velocityX, f5v);
                velocityY = _mm256_sub_pd(velocityY, f5v);

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                densityVector = _mm256_add_pd(densityVector, f6v);
                velocityY = _mm256_sub_pd(velocityY, f6v);


                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                densityVector = _mm256_add_pd(densityVector, f7v);
                velocityX = _mm256_add_pd(velocityX, f7v);
                velocityY = _mm256_sub_pd(velocityY, f7v);

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                densityVector = _mm256_add_pd(densityVector, f8v);
                velocityX = _mm256_sub_pd(velocityX, f8v);

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                densityVector = _mm256_add_pd(densityVector, f9v);

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                densityVector = _mm256_add_pd(densityVector, f10v);
                velocityX = _mm256_add_pd(velocityX, f10v);

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                densityVector = _mm256_add_pd(densityVector, f11v);
                velocityX = _mm256_sub_pd(velocityX, f11v);
                velocityY = _mm256_add_pd(velocityY, f11v);

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                densityVector = _mm256_add_pd(densityVector, f12v);
                velocityY = _mm256_add_pd(velocityY, f12v);

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                densityVector = _mm256_add_pd(densityVector, f13v);
                velocityX = _mm256_add_pd(velocityX, f13v);
                velocityY = _mm256_add_pd(velocityY, f13v);

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                densityVector = _mm256_add_pd(densityVector, f14v);
                velocityY = _mm256_sub_pd(velocityY, f14v);
                velocityZ = _mm256_add_pd(velocityZ, f14v);

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                densityVector = _mm256_add_pd(densityVector, f15v);
                velocityX = _mm256_sub_pd(velocityX, f15v);
                velocityZ = _mm256_add_pd(velocityZ, f15v);

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                densityVector = _mm256_add_pd(densityVector, f16v);
                velocityZ = _mm256_add_pd(velocityZ, f16v);

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                densityVector = _mm256_add_pd(densityVector, f17v);
                velocityX = _mm256_add_pd(velocityX, f17v);
                velocityZ = _mm256_add_pd(velocityZ, f17v);

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                densityVector = _mm256_add_pd(densityVector, f18v);
                velocityY = _mm256_add_pd(velocityY, f18v);
                velocityZ = _mm256_add_pd(velocityZ, f18v);

                velocityX = _mm256_div_pd(velocityX, densityVector);
                velocityY = _mm256_div_pd(velocityY, densityVector);
                velocityZ = _mm256_div_pd(velocityZ, densityVector);

                // compute feq and relaxate; then store to streamField
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                j += 4;
            }

            /** sequential rest */
            for (; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;

                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }

            j += (xlength[0] + 2) * (xlength[1] - 1) - (xlength[0]);
            // y = xlength[1]
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                densityVector = velocityX = velocityY = velocityZ = zeroVector;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                densityVector = _mm256_add_pd(densityVector, f0v);
                velocityY = _mm256_sub_pd(velocityY, f0v);
                velocityZ = _mm256_sub_pd(velocityZ, f0v);

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                densityVector = _mm256_add_pd(densityVector, f1v);
                velocityX = _mm256_sub_pd(velocityX, f1v);
                velocityZ = _mm256_sub_pd(velocityZ, f1v);

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                densityVector = _mm256_add_pd(densityVector, f2v);
                velocityZ = _mm256_sub_pd(velocityZ, f2v);

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                densityVector = _mm256_add_pd(densityVector, f3v);
                velocityX = _mm256_add_pd(velocityX, f3v);
                velocityZ = _mm256_sub_pd(velocityZ, f3v);

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                densityVector = _mm256_add_pd(densityVector, f4v);
                velocityY = _mm256_add_pd(velocityY, f4v);
                velocityZ = _mm256_sub_pd(velocityZ, f4v);

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                densityVector = _mm256_add_pd(densityVector, f5v);
                velocityX = _mm256_sub_pd(velocityX, f5v);
                velocityY = _mm256_sub_pd(velocityY, f5v);

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                densityVector = _mm256_add_pd(densityVector, f6v);
                velocityY = _mm256_sub_pd(velocityY, f6v);


                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                densityVector = _mm256_add_pd(densityVector, f7v);
                velocityX = _mm256_add_pd(velocityX, f7v);
                velocityY = _mm256_sub_pd(velocityY, f7v);

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                densityVector = _mm256_add_pd(densityVector, f8v);
                velocityX = _mm256_sub_pd(velocityX, f8v);

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                densityVector = _mm256_add_pd(densityVector, f9v);

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                densityVector = _mm256_add_pd(densityVector, f10v);
                velocityX = _mm256_add_pd(velocityX, f10v);

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                densityVector = _mm256_add_pd(densityVector, f11v);
                velocityX = _mm256_sub_pd(velocityX, f11v);
                velocityY = _mm256_add_pd(velocityY, f11v);

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                densityVector = _mm256_add_pd(densityVector, f12v);
                velocityY = _mm256_add_pd(velocityY, f12v);

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                densityVector = _mm256_add_pd(densityVector, f13v);
                velocityX = _mm256_add_pd(velocityX, f13v);
                velocityY = _mm256_add_pd(velocityY, f13v);

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                densityVector = _mm256_add_pd(densityVector, f14v);
                velocityY = _mm256_sub_pd(velocityY, f14v);
                velocityZ = _mm256_add_pd(velocityZ, f14v);

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                densityVector = _mm256_add_pd(densityVector, f15v);
                velocityX = _mm256_sub_pd(velocityX, f15v);
                velocityZ = _mm256_add_pd(velocityZ, f15v);

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                densityVector = _mm256_add_pd(densityVector, f16v);
                velocityZ = _mm256_add_pd(velocityZ, f16v);

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                densityVector = _mm256_add_pd(densityVector, f17v);
                velocityX = _mm256_add_pd(velocityX, f17v);
                velocityZ = _mm256_add_pd(velocityZ, f17v);

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                densityVector = _mm256_add_pd(densityVector, f18v);
                velocityY = _mm256_add_pd(velocityY, f18v);
                velocityZ = _mm256_add_pd(velocityZ, f18v);

                velocityX = _mm256_div_pd(velocityX, densityVector);
                velocityY = _mm256_div_pd(velocityY, densityVector);
                velocityZ = _mm256_div_pd(velocityZ, densityVector);

                // compute feq and relaxate; then store to streamField
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                j += 4;
            }

            /** sequential rest */
            for (; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;

                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2 * xlength[0] + 6;
        }
        /** z = xlength[2] */

        #pragma omp for schedule(static)
        for (y = 1; y <= xlength[1]; y++)
        {
            j = 1 + (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2]) + y * (xlength[0] + 2);
            for (x = 1; x <= xlength[0] - 4; x += 4)
            {
                densityVector = velocityX = velocityY = velocityZ = zeroVector;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0v = _mm256_loadu_pd(collideField + j + (xlength[0] + 2) * (xlength[1] + 3));
                densityVector = _mm256_add_pd(densityVector, f0v);
                velocityY = _mm256_sub_pd(velocityY, f0v);
                velocityZ = _mm256_sub_pd(velocityZ, f0v);

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1v = _mm256_loadu_pd(collideField + numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1);
                densityVector = _mm256_add_pd(densityVector, f1v);
                velocityX = _mm256_sub_pd(velocityX, f1v);
                velocityZ = _mm256_sub_pd(velocityZ, f1v);

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2v = _mm256_loadu_pd(collideField + 2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2));
                densityVector = _mm256_add_pd(densityVector, f2v);
                velocityZ = _mm256_sub_pd(velocityZ, f2v);

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3v = _mm256_loadu_pd(collideField + 3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1);
                densityVector = _mm256_add_pd(densityVector, f3v);
                velocityX = _mm256_add_pd(velocityX, f3v);
                velocityZ = _mm256_sub_pd(velocityZ, f3v);

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4v = _mm256_loadu_pd(collideField + 4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2));
                densityVector = _mm256_add_pd(densityVector, f4v);
                velocityY = _mm256_add_pd(velocityY, f4v);
                velocityZ = _mm256_sub_pd(velocityZ, f4v);

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5v = _mm256_loadu_pd(&collideField[5 * numberOfCells + j + xlength[0] + 3]);
                densityVector = _mm256_add_pd(densityVector, f5v);
                velocityX = _mm256_sub_pd(velocityX, f5v);
                velocityY = _mm256_sub_pd(velocityY, f5v);

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6v = _mm256_loadu_pd(&collideField[6 * numberOfCells + j + xlength[0] + 2]);
                densityVector = _mm256_add_pd(densityVector, f6v);
                velocityY = _mm256_sub_pd(velocityY, f6v);


                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7v = _mm256_loadu_pd(&collideField[7 * numberOfCells + j + xlength[0] + 1]);
                densityVector = _mm256_add_pd(densityVector, f7v);
                velocityX = _mm256_add_pd(velocityX, f7v);
                velocityY = _mm256_sub_pd(velocityY, f7v);

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8v = _mm256_loadu_pd(&collideField[8 * numberOfCells + j + 1]);
                densityVector = _mm256_add_pd(densityVector, f8v);
                velocityX = _mm256_sub_pd(velocityX, f8v);

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9v = _mm256_loadu_pd(&collideField[9 * numberOfCells + j]);
                densityVector = _mm256_add_pd(densityVector, f9v);

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10v = _mm256_loadu_pd(&collideField[10 * numberOfCells + j - 1]);
                densityVector = _mm256_add_pd(densityVector, f10v);
                velocityX = _mm256_add_pd(velocityX, f10v);

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11v = _mm256_loadu_pd(&collideField[11 * numberOfCells + j - xlength[0] - 1]);
                densityVector = _mm256_add_pd(densityVector, f11v);
                velocityX = _mm256_sub_pd(velocityX, f11v);
                velocityY = _mm256_add_pd(velocityY, f11v);

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12v = _mm256_loadu_pd(&collideField[12 * numberOfCells + j - xlength[0] - 2]);
                densityVector = _mm256_add_pd(densityVector, f12v);
                velocityY = _mm256_add_pd(velocityY, f12v);

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13v = _mm256_loadu_pd(&collideField[13 * numberOfCells + j - xlength[0] - 3]);
                densityVector = _mm256_add_pd(densityVector, f13v);
                velocityX = _mm256_add_pd(velocityX, f13v);
                velocityY = _mm256_add_pd(velocityY, f13v);

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14v = _mm256_loadu_pd(&collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)]);
                densityVector = _mm256_add_pd(densityVector, f14v);
                velocityY = _mm256_sub_pd(velocityY, f14v);
                velocityZ = _mm256_add_pd(velocityZ, f14v);

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15v = _mm256_loadu_pd(&collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1]);
                densityVector = _mm256_add_pd(densityVector, f15v);
                velocityX = _mm256_sub_pd(velocityX, f15v);
                velocityZ = _mm256_add_pd(velocityZ, f15v);

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16v = _mm256_loadu_pd(&collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)]);
                densityVector = _mm256_add_pd(densityVector, f16v);
                velocityZ = _mm256_add_pd(velocityZ, f16v);

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17v = _mm256_loadu_pd(&collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1]);
                densityVector = _mm256_add_pd(densityVector, f17v);
                velocityX = _mm256_add_pd(velocityX, f17v);
                velocityZ = _mm256_add_pd(velocityZ, f17v);

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18v = _mm256_loadu_pd(&collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)]);
                densityVector = _mm256_add_pd(densityVector, f18v);
                velocityY = _mm256_add_pd(velocityY, f18v);
                velocityZ = _mm256_add_pd(velocityZ, f18v);

                velocityX = _mm256_div_pd(velocityX, densityVector);
                velocityY = _mm256_div_pd(velocityY, densityVector);
                velocityZ = _mm256_div_pd(velocityZ, densityVector);

                // compute feq and relaxate; then store to streamField
                tmpFeq0v = _mm256_sub_pd(oneVector, _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ)), _mm256_add_pd(cs2v, cs2v)));

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f0v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 18 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f18v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f1v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 17 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f17v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 2 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f2v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 16 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f16v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 3 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f3v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 15 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f15v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 4 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f4v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 14 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f14v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 5 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f5v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 13 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f13v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 6 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f6v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 12 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f12v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProductVector = _mm256_sub_pd(velocityX, velocityY);
                tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 7 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f7v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 11 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f11v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
                tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
                _mm256_storeu_pd(streamField + 8 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f8v), _mm256_mul_pd(tauInvVector, tmpFeq1v)));
                _mm256_storeu_pd(streamField + 10 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f10v), _mm256_mul_pd(tauInvVector, tmpFeq2v)));

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                _mm256_storeu_pd(streamField + 9 * numberOfCells + j, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(oneVector, tauInvVector), f9v), _mm256_mul_pd(weight3, _mm256_mul_pd( tauInvVector, _mm256_mul_pd( densityVector, tmpFeq0v)))));
                j += 4;
            }

            /** sequential rest */
            for (; x <= xlength[0]; x++)
            {
                density = velocity[0] = velocity[1] = velocity[2] = 0;

                // load direction i = 0;
                // c_0 = (0, -1, -1)
                f0 = collideField[j + (xlength[0] + 2) * (xlength[1] + 3)];
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
                f1 = collideField[numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
                f2 = collideField[2 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2)];
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
                f3 = collideField[3 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
                f4 = collideField[4 * numberOfCells + j + (xlength[0] + 2) * (xlength[1] + 2) - (xlength[0] + 2)];
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
                f5 = collideField[5 * numberOfCells + j + xlength[0] + 3];
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
                f6 = collideField[6 * numberOfCells + j + xlength[0] + 2];
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
                f7 = collideField[7 * numberOfCells + j + xlength[0] + 1];
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
                f8 = collideField[8 * numberOfCells + j + 1];
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[9 * numberOfCells + j];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
                f10 = collideField[10 * numberOfCells + j - 1];
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
                f11 = collideField[11 * numberOfCells + j - xlength[0] - 1];
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
                f12 = collideField[12 * numberOfCells + j - xlength[0] - 2];
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
                f13 = collideField[13 * numberOfCells + j - xlength[0] - 3];
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
                f14 = collideField[14 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 1)];
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
                f15 = collideField[15 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) + 1];
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
                f16 = collideField[16 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2)];
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
                f17 = collideField[17 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 2) - 1];
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
                f18 = collideField[18 * numberOfCells + j - (xlength[0] + 2) * (xlength[1] + 3)];
                density += f18;
                velocity[1] += f18;
                velocity[2] += f18;

                velocity[0] = velocity[0] / density;
                velocity[1] = velocity[1] / density;
                velocity[2] = velocity[2] / density;

                // compute feq and relaxate; then store to streamField
                tmpFeq0 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2] ) / (2 * C_S * C_S);

                // i = 0 and i = 18
                // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
                dotProduct = - velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j] = (1 - 1./tau) * f0 + 1./tau * tmpFeq1;
                streamField[18 * numberOfCells + j] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[numberOfCells + j] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[17 * numberOfCells + j] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[2 * numberOfCells + j] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[16 * numberOfCells + j] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[3 * numberOfCells + j] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[15 * numberOfCells + j] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[4 * numberOfCells + j] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[14 * numberOfCells + j] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[5 * numberOfCells + j] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[13 * numberOfCells + j] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[6 * numberOfCells + j] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[12 * numberOfCells + j] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[7 * numberOfCells + j] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[11 * numberOfCells + j] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[8 * numberOfCells + j] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[10 * numberOfCells + j] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[9 * numberOfCells + j] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j++;
            }
            j += 2; // skip the last two cells of the current row and the first two cell of the next row
        }

    }
    break;
    }

}
#endif
