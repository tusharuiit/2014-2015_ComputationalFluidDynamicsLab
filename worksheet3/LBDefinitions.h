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
#include <stdint.h>

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

