#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#define PARAMQ 19
#define FLUID 0
#define NO_SLIP 1
#define MOVING_WALL 2
#define FREE_SLIP 3
#define INFLOW 4
#define OUTFLOW 5
#define PRESSURE_IN 6
#define PARALLEL_BOUNDARY 7
#include <stdint.h>

//#define _AVX_                     // Enables the AVX2 Streaming-and-Collision function which is slower
#define _ARBITRARYGEOMETRY_       // Enables support for scenarios of worksheet 3
//#define _DEBUG_
//#define _NOPROGRESS_              // Disables the output of progress information
#define _VTK_                       // Enables the output of VTK files
// #define _MPI_VTK_                   // Instead of VTK files create a single binary file
static const int VELOCITIESRIGHTOUT[5] = { 3, 7, 13, 17, 10 };
static const int VELOCITIESLEFTOUT[5] = { 1, 5, 11, 15, 8 };
static const int VELOCITIESBOTTOMOUT[5] = { 0, 5, 7, 14, 6 };
static const int VELOCITIESTOPOUT[5] = { 4, 11, 13, 18, 12 };
static const int VELOCITIESFRONTOUT[5] = { 14, 15, 17, 18, 16 };
static const int VELOCITIESBACKOUT[5] = { 0, 1, 3, 4, 2 };

static const int32_t LATTICEVELOCITIES[PARAMQ][3] = {{0, -1, -1},
                                          {-1, 0, -1},
                                          {0, 0, -1},
                                          {1, 0, -1},
                                          {0, 1, -1},
                                          {-1, -1, 0},
                                          {0, -1, 0},
                                          {1, -1, 0},
                                          {-1, 0, 0},
                                          {0, 0, 0},
                                          {1, 0, 0},
                                          {-1, 1, 0},
                                          {0, 1, 0},
                                          {1, 1, 0},
                                          {0, -1, 1},
                                          {-1, 0, 1},
                                          {0, 0, 1},
                                          {1, 0, 1},
                                          {0, 1, 1}
                                         };
static const double LATTICEWEIGHTS[PARAMQ] =  { 1./36,
                                         1./36,
                                         2./36,
                                         1./36,
                                         1./36,
                                         1./36,
                                         2./36,
                                         1./36,
                                         2./36,
                                         12./36,
                                         2./36,
                                         1./36,
                                         2./36,
                                         1./36,
                                         1./36,
                                         1./36,
                                         2./36,
                                         1./36,
                                         1./36
                                       };

static const int32_t LATTICEVELOCITIES2[3][PARAMQ] = {
    { 0, -1, 0, 1, 0, -1, 0, 1, -1, 0, 1, -1, 0, 1, 0, -1, 0, 1, 0 },
    { -1, 0, 0, 0, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, 0, 0, 0, 1 },
    { -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 }
};

// 1. / sqrt(3)
static const double C_S = (1.0 / 1.73205080756887729);

#endif

