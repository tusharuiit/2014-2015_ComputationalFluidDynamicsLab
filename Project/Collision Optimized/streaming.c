#include "streaming.h"
#include "LBDefinitions.h"
#include "helper.h"
#include <omp.h>
#include <unistd.h>

void doStreamingAndCollision(double *collideField, double *streamField,int *flagField,int *xlength, const double tau, int part)
{
    int x ,y ,z, j;
    double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18;
    double density, velocity[3];
    double tmpFeq0, tmpFeq1, tmpFeq2, dotProduct;
#ifdef _ARBITRARYGEOMETRIES_
    int streamMask, copyMask;
#endif // _ARBITRARYGEOMETRY_

    // Philip Neumann approach

    // During communication the innermost parts of the domain can already be processed. After sending in direction of the x-axis
    // we process a third of the domain (part == 1). After that we send in direction of the z-axis and then continue with the second
    // third of the domain (part == 2). Finally we send in direction of the y-axis and process the last third of the domain (part == 3).
    // Once we recieved all data from neighbours in y-direction we can stream from the local boundary cells as well (part == 4).

    switch(part)
    {
    case 1:
#ifdef _ARBITRARYGEOMETRIES_
        #pragma omp parallel for private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared(xlength, flagField, streamField, collideField) schedule(static)
#else
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared(xlength, streamField, collideField) schedule(static)
#endif
        for (z = 2; z <= xlength[2] / 3; z++)
        {
            j = PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + 2 * (xlength[0] + 2) + 2);
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 2; x <= xlength[0] - 1; x++)
                {
#ifdef _ARBITRARYGEOMETRIES_
                    streamMask = !flagField[j / PARAMQ];
                    copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                    density = velocity[0] = velocity[1] = velocity[2] = 0;
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                    f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                    f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                    f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                    f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                    f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                    f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                    f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                    f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                    f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[j + 9];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                    f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                    f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                    f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                    f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                    f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                    f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                    f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                    f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                    f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                    streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j += PARAMQ;
                }
                j += 4 * PARAMQ ; // skip the last cell of the current row and the first cell of the next row
            }
            j += PARAMQ * (4 * xlength[0] + 8);
        }
        break;
    case 2:
#ifdef _ARBITRARYGEOMETRIES_
        #pragma omp parallel for private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared(xlength, flagField, streamField, collideField) schedule(static)
#else
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared(xlength, streamField, collideField) schedule(static)
#endif
        for (z = xlength[2] / 3 + 1; z <= 2 * xlength[2] / 3; z++)
        {
            j = PARAMQ * (1 + 2 * (xlength[0] + 2) + (xlength[0] + 2) * (xlength[1] + 2) * z);
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 1; x <= xlength[0]; x++)
                {
#ifdef _ARBITRARYGEOMETRIES_
                    streamMask = !flagField[j / PARAMQ];
                    copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                    density = velocity[0] = velocity[1] = velocity[2] = 0;
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                    f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                    f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                    f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                    f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                    f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                    f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                    f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                    f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                    f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[j + 9];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                    f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                    f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                    f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                    f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                    f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                    f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                    f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                    f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                    f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                    streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j += PARAMQ;
                }
                j += 2 * PARAMQ; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += PARAMQ * (4 * xlength[0] + 8);
        }
        break;
    case 3:
#ifdef _ARBITRARYGEOMETRIES_
        #pragma omp parallel for private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared(xlength, flagField, streamField, collideField) schedule(static)
#else
        #pragma omp parallel for private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared( xlength, streamField, collideField) schedule(static)
#endif
        for (z = 2 * xlength[2] / 3 + 1; z <= xlength[2]; z++)
        {
            j = PARAMQ * (1 + 2 * (xlength[0] + 2) + (xlength[0] + 2) * (xlength[1] + 2) * z);
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                for (x = 1; x <= xlength[0]; x++)
                {
#ifdef _ARBITRARYGEOMETRIES_
                    streamMask = !flagField[j / PARAMQ];
                    copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                    density = velocity[0] = velocity[1] = velocity[2] = 0;
                    // load direction i = 0;
                    // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                    f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                    density += f0;
                    velocity[1] -= f0;
                    velocity[2] -= f0;

                    // load direction i = 1;
                    // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                    f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                    density += f1;
                    velocity[0] -= f1;
                    velocity[2] -= f1;

                    // load direction i = 2;
                    // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                    f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                    density += f2;
                    velocity[2] -= f2;

                    // load direction i = 3;
                    // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                    f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                    density += f3;
                    velocity[0] += f3;
                    velocity[2] -= f3;

                    // load direction i = 4;
                    // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                    f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                    f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                    density += f4;
                    velocity[1] += f4;
                    velocity[2] -= f4;

                    // load direction i = 5;
                    // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                    f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                    density += f5;
                    velocity[0] -= f5;
                    velocity[1] -= f5;

                    // load direction i = 6;
                    // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                    f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                    density += f6;
                    velocity[1] -= f6;

                    // load direction i = 7;
                    // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                    f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                    density += f7;
                    velocity[0] += f7;
                    velocity[1] -= f7;

                    // load direction i = 8;
                    // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                    f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                    density += f8;
                    velocity[0] -= f8;

                    // load direction i = 9;
                    // c_9 = (0, 0, 0)
                    f9 = collideField[j + 9];
                    density += f9;

                    // load direction i = 10;
                    // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                    f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                    density += f10;
                    velocity[0] += f10;

                    // load direction i = 11;
                    // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                    f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                    density += f11;
                    velocity[0] -= f11;
                    velocity[1] += f11;

                    // load direction i = 12;
                    // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                    f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                    density += f12;
                    velocity[1] += f12;

                    // load direction i = 13;
                    // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                    f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                    f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                    density += f13;
                    velocity[0] += f13;
                    velocity[1] += f13;

                    // load direction i = 14;
                    // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                    f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                    density += f14;
                    velocity[1] -= f14;
                    velocity[2] += f14;

                    // load direction i = 15;
                    // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                    f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                    density += f15;
                    velocity[0] -= f15;
                    velocity[2] += f15;

                    // load direction i = 16;
                    // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                    f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                    density += f16;
                    velocity[2] += f16;

                    // load direction i = 17;
                    // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                    f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                    density += f17;
                    velocity[0] += f17;
                    velocity[2] += f17;

                    // load direction i = 18;
                    // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                    f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                    f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                    streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                    // i = 1 and i = 17
                    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                    dotProduct = - velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                    streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                    // i = 2 and i = 16
                    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                    dotProduct = - velocity[2];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                    streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                    // i = 3 and i = 15
                    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                    dotProduct = velocity[0] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                    streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                    // i = 4 and i = 14
                    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                    dotProduct = velocity[1] - velocity[2];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                    streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                    // i = 5 and i = 13
                    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                    dotProduct = - velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                    streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                    // i = 6 and i = 12
                    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                    dotProduct = - velocity[1];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                    streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                    // i = 7 and i = 11
                    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                    dotProduct = velocity[0] - velocity[1];
                    tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                    streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                    // i = 8 and i = 10
                    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                    dotProduct = - velocity[0];
                    tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                    streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                    streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                    // i = 9
                    // c_9 = (0, 0, 0), w_9 = 1/3
                    streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                    j += PARAMQ;
                }
                j += 2 * PARAMQ; // skip the last two cells of the current row and the first two cell of the next row
            }
            j += PARAMQ * (4 * xlength[0] + 8);
        }
        break;
    case 4:
#ifdef _ARBITRARYGEOMETRIES_
        #pragma omp parallel private(copyMask, streamMask, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared(xlength, flagField, streamField, collideField)
#else
        #pragma omp parallel private(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, density, velocity, tmpFeq0, tmpFeq1, tmpFeq2, dotProduct, x, y, z, j) shared( xlength, streamField, collideField)
#endif
    {
        /** z = 1 */
        #pragma omp for schedule(static)
        for (y = 1; y <= xlength[1]; y++)
        {
            j = PARAMQ * (1 + (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2));
            for (x = 1; x <= xlength[0]; x++)
            {
#ifdef _ARBITRARYGEOMETRIES_
                streamMask = !flagField[j / PARAMQ];
                copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[j + 9];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += PARAMQ;
            }
            j += 2 * PARAMQ; // skip the last two cells of the current row and the first two cell of the next row
        }

        /** 1 < z <= xlength[2] / 3 **/
        #pragma omp for schedule(static)
        for (z = 2; z <= xlength[2] / 3; z++)
        {
            j = PARAMQ * (z * (xlength[0] + 2) * (xlength[1] + 2) + (xlength[0] + 2) + 1);
            // y = 1
            for (x = 1; x <= xlength[0]; x++)
            {
#ifdef _ARBITRARYGEOMETRIES_
                streamMask = !flagField[j / PARAMQ];
                copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[j + 9];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += PARAMQ;
            }
            j += 2 * PARAMQ;
            // 1 < y < xlength[1]
            for (y = 2; y <= xlength[1] - 1; y++)
            {
                // x = 1
#ifdef _ARBITRARYGEOMETRIES_
                streamMask = !flagField[j / PARAMQ];
                copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
               density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[j + 9];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += PARAMQ * (xlength[0] - 1);

                // x = xlength[0]
#ifdef _ARBITRARYGEOMETRIES_
                streamMask = !flagField[j / PARAMQ];
                copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[j + 9];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;

                j += 3 * PARAMQ;
            }
            // y = xlength[1]
            for (x = 1; x <= xlength[0]; x++)
            {
#ifdef _ARBITRARYGEOMETRIES_
                streamMask = !flagField[j / PARAMQ];
                copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[j + 9];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += PARAMQ;
            }
            j += PARAMQ * (2 * xlength[0] + 6);
        }

        /** xlength[2] / 3 < z <= xlength[2] */
        #pragma omp for schedule(static)
        for (z = xlength[2] / 3 + 1; z <= xlength[2]; z++)
        {
            // y = 1
            j = PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) * z + (xlength[0] + 2) + 1);
            for (x = 1; x <= xlength[0]; x++)
            {
#ifdef _ARBITRARYGEOMETRIES_
                streamMask = !flagField[j / PARAMQ];
                copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[j + 9];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += PARAMQ;
            }

            j += PARAMQ * ((xlength[0] + 2) * (xlength[1] - 1) - xlength[0]);
            // y = xlength[1]
            for (x = 1; x <= xlength[0]; x++)
            {
#ifdef _ARBITRARYGEOMETRIES_
                streamMask = !flagField[j / PARAMQ];
                copyMask = !streamMask;
#endif // _ARBITRARYGEOMETRIES_
                density = velocity[0] = velocity[1] = velocity[2] = 0;
                // load direction i = 0;
                // c_0 = (0, -1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)] * streamMask + collideField[j] * copyMask;
#else
                f0 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 3)];
#endif // _ARBITRARYGEOMETRIES_
                density += f0;
                velocity[1] -= f0;
                velocity[2] -= f0;

                // load direction i = 1;
                // c_1 = (-1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1] * streamMask + collideField[j + 1] * copyMask;
#else
                f1 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 1];
#endif // _ARBITRARYGEOMETRIES_
                density += f1;
                velocity[0] -= f1;
                velocity[2] -= f1;

                // load direction i = 2;
                // c_2 = (0, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2] * streamMask + collideField[2 + j] * copyMask;
#else
                f2 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2)) + 2];
#endif // _ARBITRARYGEOMETRIES_
                density += f2;
                velocity[2] -= f2;

                // load direction i = 3;
                // c_3 = (1, 0, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3] * streamMask + collideField[3 + j] * copyMask;
#else
                f3 = collideField[j + PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 3];
#endif // _ARBITRARYGEOMETRIES_
                density += f3;
                velocity[0] += f3;
                velocity[2] -= f3;

                // load direction i = 4;
                // c_4 = (0, 1, -1)
#ifdef _ARBITRARYGEOMETRIES_
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4] * streamMask + collideField[4 + j] * copyMask;
#else
                f4 = collideField[j + PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 4];
#endif // _ARBITRARYGEOMETRIES_
                density += f4;
                velocity[1] += f4;
                velocity[2] -= f4;

                // load direction i = 5;
                // c_5 = (-1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5] * streamMask + collideField[5 + j] * copyMask;
#else
                f5 = collideField[j + PARAMQ * (xlength[0] + 3) + 5];
#endif // _ARBITRARYGEOMETRIES_
                density += f5;
                velocity[0] -= f5;
                velocity[1] -= f5;

                // load direction i = 6;
                // c_6 = (0, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6] * streamMask + collideField[6 + j] * copyMask;
#else
                f6 = collideField[j + PARAMQ * (xlength[0] + 2) + 6];
#endif // _ARBITRARYGEOMETRIES_
                density += f6;
                velocity[1] -= f6;

                // load direction i = 7;
                // c_7 = (1, -1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7] * streamMask + collideField[7 + j] * copyMask;
#else
                f7 = collideField[j + PARAMQ * (xlength[0] + 1) + 7];
#endif // _ARBITRARYGEOMETRIES_
                density += f7;
                velocity[0] += f7;
                velocity[1] -= f7;

                // load direction i = 8;
                // c_8 = (-1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f8 = collideField[j + PARAMQ + 8] * streamMask + collideField[8 + j] * copyMask;
#else
                f8 = collideField[j + PARAMQ + 8];
#endif // _ARBITRARYGEOMETRIES_
                density += f8;
                velocity[0] -= f8;

                // load direction i = 9;
                // c_9 = (0, 0, 0)
                f9 = collideField[j + 9];
                density += f9;

                // load direction i = 10;
                // c_10 = (1, 0, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f10 = collideField[j - PARAMQ + 10] * streamMask + collideField[10 + j] * copyMask;
#else
                f10 = collideField[j - PARAMQ + 10];
#endif // _ARBITRARYGEOMETRIES_
                density += f10;
                velocity[0] += f10;

                // load direction i = 11;
                // c_11 = (-1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11] * streamMask + collideField[11 + j] * copyMask;
#else
                f11 = collideField[j -  PARAMQ * (xlength[0] + 1) + 11];
#endif // _ARBITRARYGEOMETRIES_
                density += f11;
                velocity[0] -= f11;
                velocity[1] += f11;

                // load direction i = 12;
                // c_12 = (0, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12] * streamMask + collideField[12 + j] * copyMask;
#else
                f12 = collideField[j - PARAMQ * (xlength[0] + 2) + 12];
#endif // _ARBITRARYGEOMETRIES_
                density += f12;
                velocity[1] += f12;

                // load direction i = 13;
                // c_13 = (1, 1, 0)
#ifdef _ARBITRARYGEOMETRIES_
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13] * streamMask + collideField[13 + j] * copyMask;
#else
                f13 = collideField[j - PARAMQ * (xlength[0] + 3) + 13];
#endif // _ARBITRARYGEOMETRIES_
                density += f13;
                velocity[0] += f13;
                velocity[1] += f13;

                // load direction i = 14;
                // c_14 = (0, -1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14] * streamMask + collideField[14 + j] * copyMask;
#else
                f14 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 1) + 14];
#endif // _ARBITRARYGEOMETRIES_
                density += f14;
                velocity[1] -= f14;
                velocity[2] += f14;

                // load direction i = 15;
                // c_15 = (-1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15] * streamMask + collideField[15 + j] * copyMask;
#else
                f15 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) - 1) + 15];
#endif // _ARBITRARYGEOMETRIES_
                density += f15;
                velocity[0] -= f15;
                velocity[2] += f15;

                // load direction i = 16;
                // c_16 = (0, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16] * streamMask + collideField[16 + j] * copyMask;
#else
                f16 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 2) + 16];
#endif // _ARBITRARYGEOMETRIES_
                density += f16;
                velocity[2] += f16;

                // load direction i = 17;
                // c_17 = (1, 0, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17] * streamMask + collideField[17 + j] * copyMask;
#else
                f17 = collideField[j - PARAMQ * ((xlength[0] + 2) * (xlength[1] + 2) + 1) + 17];
#endif // _ARBITRARYGEOMETRIES_
                density += f17;
                velocity[0] += f17;
                velocity[2] += f17;

                // load direction i = 18;
                // c_18 = (0, 1, 1)
#ifdef _ARBITRARYGEOMETRIES_
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18] * streamMask + collideField[18 + j] * copyMask;
#else
                f18 = collideField[j - PARAMQ * (xlength[0] + 2) * (xlength[1] + 3) + 18];
#endif // _ARBITRARYGEOMETRIES_
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
                streamField[j + 18] = (1 - 1./tau) * f18 + 1./tau * tmpFeq2;

                // i = 1 and i = 17
                // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
                dotProduct = - velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 1] = (1 - 1./tau) * f1 + 1./tau * tmpFeq1;
                streamField[j + 17] = (1 - 1./tau) * f17 + 1./tau * tmpFeq2;

                // i = 2 and i = 16
                // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
                dotProduct = - velocity[2];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 2] = (1 - 1./tau) * f2 + 1./tau * tmpFeq1;
                streamField[j + 16] = (1 - 1./tau) * f16 + 1./tau * tmpFeq2;

                // i = 3 and i = 15
                // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
                dotProduct = velocity[0] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 3] = (1 - 1./tau) * f3 + 1./tau * tmpFeq1;
                streamField[j + 15] = (1 - 1./tau) * f15 + 1./tau * tmpFeq2;

                // i = 4 and i = 14
                // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
                dotProduct = velocity[1] - velocity[2];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 4] = (1 - 1./tau) * f4 + 1./tau * tmpFeq1;
                streamField[j + 14] = (1 - 1./tau) * f14 + 1./tau * tmpFeq2;

                // i = 5 and i = 13
                // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
                dotProduct = - velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 5] = (1 - 1./tau) * f5 + 1./tau * tmpFeq1;
                streamField[j + 13] = (1 - 1./tau) * f13 + 1./tau * tmpFeq2;

                // i = 6 and i = 12
                // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
                dotProduct = - velocity[1];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 6] = (1 - 1./tau) * f6+ 1./tau * tmpFeq1;
                streamField[j + 12] = (1 - 1./tau) * f12 + 1./tau * tmpFeq2;

                // i = 7 and i = 11
                // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
                dotProduct = velocity[0] - velocity[1];
                tmpFeq1 = 1./36 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./36 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 7] = (1 - 1./tau) * f7 + 1./tau * tmpFeq1;
                streamField[j + 11] = (1 - 1./tau) * f11 + 1./tau * tmpFeq2;

                // i = 8 and i = 10
                // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
                dotProduct = - velocity[0];
                tmpFeq1 = 1./18 * density * (tmpFeq0 + dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                tmpFeq2 = 1./18 * density * (tmpFeq0 - dotProduct/(C_S * C_S) + dotProduct * dotProduct / (2 * C_S * C_S * C_S * C_S));
                streamField[j + 8] = (1 - 1./tau) * f8 + 1./tau * tmpFeq1;
                streamField[j + 10] = (1 - 1./tau) * f10 + 1./tau * tmpFeq2;

                // i = 9
                // c_9 = (0, 0, 0), w_9 = 1/3
                streamField[j + 9] = (1 - 1. / tau) * f9 + 1. / (3 * tau) * density * tmpFeq0;


                j += PARAMQ;
            }
            j += PARAMQ * (2 * xlength[0] + 6);
        }
    }
    break;
    }
}
