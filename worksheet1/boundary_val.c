#include "helper.h"
#include "boundary_val.h"

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
)
{
    /** Set No-slip boundary condition*/
    int i;

    for(i = 1; i <= jmax; ++i)
    {
    	U[0][i] = 0.0;
    	U[imax][i] = 0.0;

    	V[0][i] = -V[1][i];
    	V[imax + 1][i] = -V[imax][i];
    }

    for(i = 1; i <= imax; ++i)
    {
    	V[i][0] = 0.0;
    	V[i][jmax] = 0.0;

    	U[i][0] = -U[i][1];
    	U[i][jmax + 1] = 2.0 - U[i][jmax];
    }

}
