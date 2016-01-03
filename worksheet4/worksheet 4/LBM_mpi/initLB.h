#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *xlength,                       /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *bddParams,                  /* an array of all parameters needed for boundary treatment, parameter name depends on the scenario */
    int *iProc,                         /* number of processes in X direction */
    int *jProc,                         /* number of processes in Y direction */
    int *kProc,                         /* number of processes in Z direction */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    char *problem,                      /* name of the scenario to simulate*/
    char *pgmInput,                      /* pgm File for the scenario used in the tilted plate case */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[]                        /* argv[1] shall contain the path to the config file */
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(double *collideField, double *streamField,int *flagField, int *xlength, char *problem, char *pgmInput, int rank, int iProc, int jProc, int kProc);

#endif

