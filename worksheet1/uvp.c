#include <math.h>
#include "helper.h"
#include "uvp.h"

/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
    /** Initialization */
    int i, j;
    double val_1, val_2, val_3; // store 3 values inside min( , , )
    double temp_val1, temp_val2; // store temporary variables
    double umax = 0.0, vmax = 0.0; // store the maximal values of u and v in the internal domain

    temp_val1 = 1 / dx;
    temp_val2 = 1 / dy;

    /** Calculate the tree entries */
    val_1 = Re / 2 / ( temp_val1 * temp_val1 + temp_val2 * temp_val2 );

    for(i = 1; i <= imax; ++i)
    	for (j = 1; j <= jmax; ++j)
    	{
    	    if( fabs(U[i][j]) > umax ) umax = fabs(U[i][j]);
    	    if( fabs(V[i][j]) > vmax ) vmax = fabs(V[i][j]);
    	}

    val_2 = dx / umax;
    val_3 = dy / vmax;

    (*dt) = (val_1 > val_2) ? val_2 : val_1;
    (*dt) = ( (*dt) > val_3) ? val_3 : (*dt);
    (*dt) = tau * (*dt);

}


/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
)
{
    int i, j;
    double discrLaplU = 0.0;
    double discrLaplV = 0.0;
    double dx_u2 = 0.0;
    double dy_uv = 0.0;
    double dy_v2 = 0.0;
    double dx_uv = 0.0;

    for (i = 0; i < imax; i++)
    {
        for (j = 1; j <= jmax; j++)
        {
            if (i == 0)
            {
                F[i][j] = U[i][j];
                F[imax][j] = U[imax][j];
            }
            else
            {
                discrLaplU = 1./(dx*dx) * (U[i+1][j] - 2*U[i][j] + U[i-1][j]) +
                             1./(dy*dy) * (U[i][j+1] - 2*U[i][j] + U[i][j-1]);
                dx_u2 = 1./(4*dx) * (pow(U[i][j] + U[i+1][j],2) - pow(U[i-1][j] + U[i][j],2))
                        + alpha/(4.0*dx) * (fabs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j]) - fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]));
                dy_uv = 1./(4*dy) * ((U[i][j+1]+U[i][j])*(V[i][j]+V[i+1][j]) - (U[i][j]+U[i][j-1])*(V[i][j-1]+V[i+1][j-1]))
                        + alpha/(4.0*dy) * (fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1]) - fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]));
                F[i][j] = U[i][j] + dt*(1./Re * discrLaplU - dx_u2 - dy_uv + GX);
            }
        }
    }
    
    for (i = 1; i <= imax; i++)
    {
        for (j = 0; j < jmax; j++)
        {
            if (j == 0)
            {
                G[i][0] = V[i][0];
                G[i][jmax] = V[i][jmax];
            }
            else
            {
                discrLaplV = 1./(dx*dx) * (V[i+1][j] - 2*V[i][j] + V[i-1][j]) +
                             1./(dy * dy) * (V[i][j+1] - 2*V[i][j] + V[i][j-1]);
                dy_v2 = 1./(4*dy) * (pow(V[i][j] + V[i][j+1],2) - pow(V[i][j] + V[i][j-1],2)) +
                        alpha/(4.0*dy) * (fabs(V[i][j] + V[i][j+1])*(V[i][j] - V[i][j+1]) - fabs(V[i][j] + V[i][j-1])*(V[i][j] - V[i][j-1]));
                dx_uv = 1./(4*dx) * ((U[i][j]+U[i][j+1])*(V[i+1][j]+V[i][j]) -
                                     (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])) +
                        alpha/(4.0*dx) * (fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j]) - fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]));
                G[i][j] = V[i][j] + dt*(1./Re * discrLaplV - dx_uv - dy_v2 + GY);
            }
        }
    }
}


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
	int i, j;

	for(i = 1; i <= imax; ++i)
		for(j = 1; j <= jmax; ++j)
		{
			RS[i][j] = 1./dt * ( (F[i][j] - F[i - 1][j]) / dx + (G[i][j] - G[i][j - 1]) / dy );
		}

}

/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
    int i, j;

    for(i = 1; i <= imax; ++i)
    	for(j = 1; j <= jmax; ++j)
    	{
    		/* boundary values will be overwritten anyway */
		U[i][j] = F[i][j] - dt / dx * (P[i + 1][j] - P[i][j]);
		V[i][j] = G[i][j] - dt / dy * (P[i][j + 1] - P[i][j]);

    	}


}
