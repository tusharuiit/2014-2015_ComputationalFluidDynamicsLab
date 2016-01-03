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

//#define DEBUG   // check correctness by comparison with pregenerated test files
#define _ARBITRARYGEOMETRIES_   // Allows for scenarios other than driven cavity
//#define _NOPROGRESS_    // Set this flag to suppress the display of progress
#define _VTK_           // enables VTK output
//#define _VTK_BINARY_    // enables the use of MPI_IO to output a binary file of all collected outputs
                          // which can later be converted to csv output
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

