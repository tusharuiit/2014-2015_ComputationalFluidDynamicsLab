#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <immintrin.h>


 void computeDensity(const double *const currentCell, double *density)
{
    *density = 0.0;
    int i;
    for (i = 0; i < PARAMQ; i++)
        *density += currentCell[i];
}

void computeDensityAVX(const double * const currentCell, double *density)
{
    double const zero = 0.0;
    __m256d vsum = _mm256_broadcast_sd(&zero);
    int i;
    for (i = 0; i < PARAMQ - 3; i += 4)
    {
        __m256d v = _mm256_loadu_pd(&currentCell[i]);
        vsum = _mm256_add_pd(vsum, v);
    }
    vsum = _mm256_hadd_pd(vsum, vsum);
    *density = ((double*)&vsum)[0] + ((double*)&vsum)[2];
    for (; i < PARAMQ; i++)
        *density += currentCell[i];


}
void computeDensitySSE(const double * const currentCell, double *density)
{
    __m128d vsum = _mm_set1_pd(0.0);
    int i;
    for (i = 0; i < PARAMQ - 1; i += 2)
    {
        __m128d v = _mm_loadu_pd(&currentCell[i]);
        vsum = _mm_add_pd(vsum, v);
    }
    vsum = _mm_hadd_pd(vsum, vsum);
    _mm_storeh_pd(density, vsum);
    if (i < PARAMQ)
    {
        *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity)
{
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;
    for (int i = 0; i < PARAMQ; i++)
    {
        velocity[0] += currentCell[i] * LATTICEVELOCITIES[i][0];
    }

    for (int i = 0; i < PARAMQ; i++)
    {
        velocity[1] += currentCell[i] * LATTICEVELOCITIES[i][1];
    }

    for (int i = 0; i < PARAMQ; i++)
    {
        velocity[2] += currentCell[i] * LATTICEVELOCITIES[i][2];
    }

    velocity[0] = velocity[0] / (*density);
    velocity[1] = velocity[1] / (*density);
    velocity[2] = velocity[2] / (*density);

}
void computeVelocityAVX(const double * const currentCell, const double * const density, double *velocity)
{
    __m256d v0, v1, v2;
    int i;
    double const zero = 0.0;
    v0 = v1 = v2 = _mm256_broadcast_sd(&zero);
    for (i = 0; i < PARAMQ - 3; i += 4)
    {
        __m256d vc, vl0, vl1, vl2;
        __m128i vtemp;
        vc = _mm256_loadu_pd(&currentCell[i]);
        vc = _mm256_loadu_pd(&currentCell[i]);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[0][i]);
        vl0 = _mm256_cvtepi32_pd(vtemp);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[1][i]);
        vl1 = _mm256_cvtepi32_pd(vtemp);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[2][i]);
        vl2 = _mm256_cvtepi32_pd(vtemp);
        v0 = _mm256_add_pd(v0, _mm256_mul_pd(vc, vl0));
        v1 = _mm256_add_pd(v1, _mm256_mul_pd(vc, vl1));
        v2 = _mm256_add_pd(v2, _mm256_mul_pd(vc, vl2));
    }
    v0 = _mm256_hadd_pd(v0, v0);
    v1 = _mm256_hadd_pd(v1, v1);
    v2 = _mm256_hadd_pd(v2, v2);
    velocity[0] = ((double*)&v0)[0] + ((double*)&v0)[2];
    velocity[1] = ((double*)&v1)[0] + ((double*)&v1)[2];
    velocity[2] = ((double*)&v2)[0] + ((double*)&v2)[2];
    for (; i < PARAMQ; i++)
    {
        velocity[0] += currentCell[i] * LATTICEVELOCITIES2[0][i];
        velocity[1] += currentCell[i] * LATTICEVELOCITIES2[1][i];
        velocity[2] += currentCell[i] * LATTICEVELOCITIES2[2][i];
    }
    velocity[0] = velocity[0] / (*density);
    velocity[1] = velocity[1] / (*density);
    velocity[2] = velocity[2] / (*density);

}
void computeVelocitySSE(const double * const currentCell, const double * const density, double *velocity)
{
    __m128d v0, v1, v2;
    int i;
    v0 = v1 = v2 = _mm_setzero_pd();
    for (i = 0; i < PARAMQ - 1; i += 2)
    {
        __m128d vc, vl0, vl1, vl2;
        __m128i vtemp;

        vc = _mm_loadu_pd(&currentCell[i]);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[0][i]);
        vl0 = _mm_cvtepi32_pd(vtemp);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[1][i]);
        vl1 = _mm_cvtepi32_pd(vtemp);
        vtemp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[2][i]);
        vl2 = _mm_cvtepi32_pd(vtemp);
        v0 = _mm_add_pd(v0, _mm_mul_pd(vc, vl0));
        v1 = _mm_add_pd(v1, _mm_mul_pd(vc, vl1));
        v2 = _mm_add_pd(v2, _mm_mul_pd(vc, vl2));
    }
    v0 = _mm_hadd_pd(v0, v0);
    v1 = _mm_hadd_pd(v1, v1);
    v2 = _mm_hadd_pd(v2, v2);
    _mm_store_sd (&velocity[0], v0);
    _mm_store_sd (&velocity[1], v1);
    _mm_store_sd (&velocity[2], v2);
    if (i < PARAMQ)
    {
        velocity[0] += currentCell[i] * LATTICEVELOCITIES2[0][i];
        velocity[1] += currentCell[i] * LATTICEVELOCITIES2[1][i];
        velocity[2] += currentCell[i] * LATTICEVELOCITIES2[2][i];
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
void computeFeqAVX(const double * const density, const double * const velocity, double *feq)
{
    double velocityMagnitudeSqr, dotProduct;
    double const one = 1.0;
    __m256d dotProductVector, velocityX, velocityY, velocityZ, latWeights, result;
    __m128i temp;
    __m256d vl0, vl1, vl2;
    int i;
    double c_s_square = C_S * C_S;
    velocityX = _mm256_broadcast_sd(velocity);
    velocityY = _mm256_broadcast_sd(velocity + 1);
    velocityZ = _mm256_broadcast_sd(velocity + 2);

    velocityMagnitudeSqr = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
    __m256d c_s_squareVector = _mm256_broadcast_sd(&c_s_square);
    __m256d velocityMagnitudeSqrVector = _mm256_broadcast_sd(&velocityMagnitudeSqr);
    __m256d densityVector = _mm256_broadcast_sd(density);
    __m256d oneVector = _mm256_broadcast_sd(&one);
    for (i = 0; i < PARAMQ - 3; i += 4)
    {
        latWeights = _mm256_loadu_pd(&LATTICEWEIGHTS[i]);
        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[0][i]);
        vl0 = _mm256_cvtepi32_pd(temp);
        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[1][i]);
        vl1 = _mm256_cvtepi32_pd(temp);
        temp = _mm_loadu_si128((__m128i *)&LATTICEVELOCITIES2[2][i]);
        vl2 = _mm256_cvtepi32_pd(temp);
        dotProductVector = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(velocityX, vl0), _mm256_mul_pd(velocityY, vl1)), _mm256_mul_pd(velocityZ, vl2));
        result = _mm256_mul_pd(latWeights, _mm256_mul_pd(densityVector, _mm256_sub_pd(_mm256_add_pd(oneVector, _mm256_add_pd(_mm256_div_pd(dotProductVector, c_s_squareVector), _mm256_div_pd(_mm256_mul_pd(dotProductVector, dotProductVector), _mm256_add_pd(_mm256_mul_pd(c_s_squareVector, c_s_squareVector),_mm256_mul_pd(c_s_squareVector, c_s_squareVector))))), _mm256_div_pd(velocityMagnitudeSqrVector, _mm256_add_pd(c_s_squareVector, c_s_squareVector)))));
        _mm256_storeu_pd(&feq[i], result);
    }
    for (; i < PARAMQ; i++)
    {
        dotProduct = velocity[0] * LATTICEVELOCITIES[i][0] + velocity[1] * LATTICEVELOCITIES[i][1] + velocity[2] * LATTICEVELOCITIES[i][2];
        feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + dotProduct / c_s_square + dotProduct * dotProduct / (2  * c_s_square * c_s_square) - velocityMagnitudeSqr / (2 * c_s_square));
    }
}
