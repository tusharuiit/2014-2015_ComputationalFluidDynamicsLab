#ifndef _STREAMING_H_
#define _STREAMING_H_

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(double *collideField, double *streamField,int *flagField,int *xlength);
void doStreamingSSE(double *collideField, double *streamField,int *flagField,int *xlength, double * densities, double ** velocities);
void doStreamingAVX(double *collideField, double *streamField,int *flagField,int *xlength);
void doStreamingAndCollision(double *collideField, double *streamField,int *flagField,int *xlength, const double tau, int part);
void doStreamingAndCollisionAVX(double *collideField, double *streamField,int *flagField,int *xlength, const double tau, int part);
void doStreamingAVXv2(double *collideField, double *streamField,int *flagField,int *xlength, double * densities, double ** velocities);

#endif

