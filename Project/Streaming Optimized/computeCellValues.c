#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include <stdlib.h>

#include "helper.h"


void computeDensity(const double *const currentCell, double *density, const int tot_cells)
{
    *density = 0.0;
    int i;
    for (i = 0; i < PARAMQ; i++)
        *density += *(currentCell + i * tot_cells);
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity, const int tot_cells)
{
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;
    for (int i = 0; i < PARAMQ; i++)
    {
        velocity[0] += *(currentCell + i * tot_cells) * LATTICEVELOCITIES[i][0];

        velocity[1] += *(currentCell + i * tot_cells) * LATTICEVELOCITIES[i][1];

        velocity[2] += *(currentCell + i * tot_cells) * LATTICEVELOCITIES[i][2];
    }

    velocity[0] = velocity[0] / (*density);
    velocity[1] = velocity[1] / (*density);
    velocity[2] = velocity[2] / (*density);

}

void computeFeq(const double * const density, const double * const velocity, double *feq)
{
    double velocityMagnitudeSqr, dotProduct;
    double c_s_square = C_S * C_S;
    velocityMagnitudeSqr = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
    for (int i = 0; i < PARAMQ; i++)
    {
        dotProduct = velocity[0] * LATTICEVELOCITIES[i][0] + velocity[1] * LATTICEVELOCITIES[i][1] + velocity[2] * LATTICEVELOCITIES[i][2];
        feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + dotProduct / c_s_square + dotProduct * dotProduct / (2  * c_s_square * c_s_square) - velocityMagnitudeSqr / (2 * c_s_square));
    }
}
//
//void computeFeqAVXv2(double density[4], double velocity[3][4], double feq[PARAMQ][4])
//{
//    const double c_s_square = C_S * C_S;
//    __m256d cs2v = _mm256_broadcast_sd(&c_s_square);
//    __m256d latticeWeightsVector, latticeVelocityX, latticeVelocityY, latticeVelocityZ, dotProductVector, tmpFeq1v, tmpFeq2v;
//    const double zero = 0.0;
//    const double one = 1.0;
//    __m256d oneVector = _mm256_broadcast_sd(&one);
//    __m256d zeroVector = _mm256_broadcast_sd(&zero);
//    __m256d weight36 = _mm256_broadcast_sd(&LATTICEWEIGHTS[7]);
//    __m256d weight18 = _mm256_broadcast_sd(&LATTICEWEIGHTS[8]);
//    __m256d weight3 = _mm256_broadcast_sd(&LATTICEWEIGHTS[9]);
//
//    __m256d densityVector = _mm256_load_pd(density);
//    __m256d velocityX = _mm256_load_pd(velocity[0]);
//    __m256d velocityY = _mm256_load_pd(velocity[1]);
//    __m256d velocityZ = _mm256_load_pd(velocity[2]);
//    __m256d tmpFeq0v = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, velocityX), _mm256_mul_pd(velocityY, velocityY)), _mm256_mul_pd(velocityZ, velocityZ));
//    tmpFeq0v = _mm256_sub_pd( oneVector, _mm256_div_pd(tmpFeq0v, _mm256_add_pd(cs2v, cs2v)));
//
//
//
//    // i = 0 and i = 18
//    // c_0 = (0, -1, -1), c_18 = (0, 1, 1), w_0 = w_18 = 1/36
//    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityY, velocityZ));
//    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[0], tmpFeq1v );
//    _mm256_store_pd( feq[18], tmpFeq2v);
//
//    // i = 1 and i = 17
//    // c_1 = (-1, 0, -1), c_17 = (1, 0, 1), w_1 = w_17 = 1/36
//    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityZ));
//    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[1], tmpFeq1v );
//    _mm256_store_pd( feq[17], tmpFeq2v);
//
//    // i = 2 and i = 16
//    // c_2 = (0, 0, -1), c_16 = (0, 0, 1), w_2 = w_16 = 1/18
//    dotProductVector = _mm256_sub_pd(zeroVector, velocityZ);
//    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[2], tmpFeq1v );
//    _mm256_store_pd( feq[16], tmpFeq2v);
//
//    // i = 3 and i = 15
//    // c_3 = (1, 0, -1), c_16 = (-1, 0, 1), w_3 = w_15 = 1/36
//    dotProductVector = _mm256_sub_pd(velocityX, velocityZ);
//    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[3], tmpFeq1v );
//    _mm256_store_pd( feq[15], tmpFeq2v);
//
//    // i = 4 and i = 14
//    // c_4 = (0, 1, -1), c_14 = (0, -1, 1), w_4 = w_14 = 1/36
//    dotProductVector = _mm256_sub_pd(velocityY, velocityZ);
//    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[4], tmpFeq1v );
//    _mm256_store_pd( feq[14], tmpFeq2v);
//
//    // i = 5 and i = 13
//    // c_5 = (-1, -1, 0), c_13 = (1, 1, 0), w_5 = w_13 = 1/36
//    dotProductVector = _mm256_sub_pd(zeroVector, _mm256_add_pd(velocityX, velocityY));
//    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[5], tmpFeq1v );
//    _mm256_store_pd( feq[13], tmpFeq2v);
//
//    // i = 6 and i = 12
//    // c_6 = (0, -1, 0), c_14 = (0, 1, 0), w_6 = w_12 = 1/18
//    dotProductVector = _mm256_sub_pd(zeroVector, velocityY);
//    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[6], tmpFeq1v );
//    _mm256_store_pd( feq[12], tmpFeq2v);
//
//    // i = 7 and i = 11
//    // c_7 = (1, -1, 0), c_11 = (-1, 1, 0), w_7 = w_11 = 1/36
//    dotProductVector = _mm256_sub_pd(velocityX, velocityY);
//    tmpFeq1v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight36, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[7], tmpFeq1v );
//    _mm256_store_pd( feq[11], tmpFeq2v);
//
//    // i = 8 and i = 10
//    // c_8 = (-1, 0, 0), c_10 = (1, 0, 0), w_8 = w_10 = 1/18
//    dotProductVector = _mm256_sub_pd(zeroVector, velocityX);
//    tmpFeq1v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_add_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//    tmpFeq2v = _mm256_mul_pd(weight18, _mm256_mul_pd( densityVector, _mm256_add_pd(_mm256_sub_pd(tmpFeq0v, _mm256_div_pd(dotProductVector, cs2v)), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(cs2v, cs2v), _mm256_mul_pd(cs2v, cs2v))))));
//
//    _mm256_store_pd( feq[8], tmpFeq1v );
//    _mm256_store_pd( feq[10], tmpFeq2v);
//
//    // i = 9
//    // c_9 = (0, 0, 0), w_9 = 1/3
//    _mm256_store_pd( feq[9], _mm256_mul_pd(weight3, _mm256_mul_pd( densityVector, tmpFeq0v)));
//
//}
//void computeFeqAVX(const double * const density, const double * const velocity, double *feq)
//{
//    double velocityMagnitudeSqr, dotProduct;
//    double const one = 1.0;
//    __m256d dotProductVector, velocityX, velocityY, velocityZ, latWeights, result;
//    __m128i temp;
//    __m256d vl0, vl1, vl2;
//    int i;
//    double c_s_square = C_S * C_S;
//    velocityX = _mm256_broadcast_sd(velocity);
//    velocityY = _mm256_broadcast_sd(velocity + 1);
//    velocityZ = _mm256_broadcast_sd(velocity + 2);
//
//    velocityMagnitudeSqr = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
//    __m256d c_s_squareVector = _mm256_broadcast_sd(&c_s_square);
//    __m256d velocityMagnitudeSqrVector = _mm256_broadcast_sd(&velocityMagnitudeSqr);
//    __m256d densityVector = _mm256_broadcast_sd(density);
//    __m256d oneVector = _mm256_broadcast_sd(&one);
//    for (i = 0; i < PARAMQ - 3; i += 4)
//    {
//        latWeights = _mm256_loadu_pd(&LATTICEWEIGHTS[i]);
//        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[0][i]);
//        vl0 = _mm256_cvtepi32_pd(temp);
//        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[1][i]);
//        vl1 = _mm256_cvtepi32_pd(temp);
//        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[2][i]);
//        vl2 = _mm256_cvtepi32_pd(temp);
//        dotProductVector = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, vl0), _mm256_mul_pd(velocityY, vl1)), _mm256_mul_pd(velocityZ, vl2));
//        result = _mm256_mul_pd(latWeights, _mm256_mul_pd(densityVector, _mm256_sub_pd(_mm256_add_pd(oneVector, _mm256_add_pd(_mm256_div_pd(dotProductVector, c_s_squareVector), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(c_s_squareVector, c_s_squareVector),_mm256_mul_pd(c_s_squareVector, c_s_squareVector))))), _mm256_div_pd(velocityMagnitudeSqrVector, _mm256_add_pd(c_s_squareVector, c_s_squareVector)))));
//        _mm256_store_pd(&feq[i], result);
//    }
//    for (; i < PARAMQ; i++)
//    {
//        dotProduct = velocity[0] * LATTICEVELOCITIES[i][0] + velocity[1] * LATTICEVELOCITIES[i][1] + velocity[2] * LATTICEVELOCITIES[i][2];
//        feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + dotProduct / c_s_square + dotProduct * dotProduct / (2  * c_s_square * c_s_square) - velocityMagnitudeSqr / (2 * c_s_square));
//    }
//}
