#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
 void computeDensity(const double *const currentCell, double *density)
{
    *density = 0.0;
    int i;
    for (i = 0; i < PARAMQ; i++)
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

