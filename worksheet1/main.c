#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include <time.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args)
{

    /** Program variables */
    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dT, dx, dy, alpha, omg, tau, eps, dt_value; // dT: fixed dt
    double t, dt = 0.0, it, res = 10.0, **U, **V, **P, **F, **G, **RS;
    double output_t = 0.0;
    int imax, jmax, itermax, n, output_count = 0;


    /** statistics variables */
    int sumIter, maxIter, minIter;
    double minTimeStep, maxTimestep;
    time_t start, stop;

    /** Command line check */
    if ( argn != 3 )
    {
        fprintf(stderr, "Usage: %s <input_file> <output_directory>\n", args[0]);
        return -1;
    }

    char *file_in = args[1];
    char file_out[500];
    sprintf(file_out, "%s/%s", args[2], args[1]);

    /** read problem parameters */
    read_parameters( file_in, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength,
                     &ylength, &dT, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
                     &itermax, &eps, &dt_value );

    /** set t and n */
    t = 0.0;
    n = 0;

    /** initialize u v p */
    U = matrix(0, imax, 0, jmax + 1);
    V = matrix(0, imax + 1, 0, jmax);
    P = matrix(0, imax + 1, 0, jmax + 1);
    RS = matrix(0, imax + 1, 0, jmax + 1);
    F = matrix(0, imax, 0, jmax);
    G = matrix(0, imax, 0, jmax);
    init_uvp( UI, VI, PI, imax, jmax, U, V, P );

    /** initialize statistics variables */
    minIter = itermax;
    minTimeStep = t_end;
    sumIter = 0;
    maxIter = 0;
    maxTimestep = 0;

    // save initial conditions in vtk
    write_vtkFile(file_out,0,xlength,ylength,imax,jmax,dx,dy,U,V,P);
    output_count++;

    time(&start);
    while (t < t_end)
    {

        if( tau < 0.0 )
        {
            dt = dT;
        }
        else
        {
            calculate_dt( Re, tau, &dt, dx, dy, imax, jmax, U, V );
        }
        if (dt < minTimeStep)
            minTimeStep = dt;
        if (dt > maxTimestep)
            maxTimestep = dt;

        boundaryvalues( imax, jmax, U, V );
        calculate_fg( Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G );
        calculate_rs( dt, dx, dy, imax, jmax, F, G, RS);

        it = 0;
        res = 10.0;
        while (it < itermax && res > eps)
        {
            sor( omg, dx, dy, imax, jmax, P, RS, &res );
            it++;
        }
        if (res > eps)
            fprintf(stderr, "WARNING: SOR did not converge in timestep %d, residual: %f\n", n + 1, res);
        if (it > maxIter)
            maxIter = it;
        if (it < minIter)
            minIter = it;
        sumIter += it;

        calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);

        // generate vtk file
        t = t + dt;
        if(t >= output_t)
        {
            write_vtkFile(file_out, output_count, xlength, ylength, imax, jmax, dx, dy, U, V, P);
            output_t = output_t + dt_value;
            output_count++;
        }
        n++;
    }
    time(&stop);

    printf("Execution time: %.1f seconds. \n", difftime(stop, start));
    printf("Min. timestep: %.3f     Max. timestep: %.3f    Avg. timestep: %.3f\n", minTimeStep, maxTimestep, t_end/n);
    printf("Min. # SOR iterations: %d     Max. # SOR iterations: %d     Avg: %.3f\n", minIter, maxIter, sumIter/(double) n);

    // free dynamic memory
    free_matrix( U, 0, imax , 0, jmax + 1 );
    free_matrix( V, 0, imax + 1, 0, jmax );
    free_matrix( P, 0, imax + 1, 0, jmax + 1 );
    free_matrix( RS, 0, imax + 1, 0, jmax + 1 );
    free_matrix( F, 0, imax, 0, jmax );
    free_matrix( G, 0, imax, 0, jmax );

    return 0;
}
