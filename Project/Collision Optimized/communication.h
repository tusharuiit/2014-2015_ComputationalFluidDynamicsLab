#ifndef _COMMUNICATION_H_
#define _COMMUNICATION_H_
#include <mpi.h>

void communicateBoundaryValues (  MPI_Datatype * MPISendTypes, MPI_Datatype * MPIRecvTypes, double * collideField, int * neighbours, MPI_Request * MPISendReq, MPI_Request * MPIRecvReq, int axis);
void computeNeighbours( int iProc, int jProc, int kProc, int * neighbours );
void computePosition( int iProc, int jProc, int kProc, int * iCoord, int * jCoord, int * kCoord );
void computeLocalDomainSize( int * xlength, int * local_xlength, int iProc, int jProc, int kProc );
void checkRequestCompletion( MPI_Request * MPISendReq, MPI_Request * MPIRecvReq, int axis, int * neighbours);
void getSendRecvCount( const int * const local_xlength, const int * const neighbours, const int * const flagField, int * sendCount, int * recvCount );
void getSendRecvIndices( const int * const local_xlength, const int * const neighbours, const int * const flagField, int ** sendIndices, int ** recvIndices );
void initialiseMPITypes( const int * const local_xlength, const int * const neighbours, const int * const flagField, MPI_Datatype * MPISendTypes, MPI_Datatype * MPIRecvTypes );
#endif
