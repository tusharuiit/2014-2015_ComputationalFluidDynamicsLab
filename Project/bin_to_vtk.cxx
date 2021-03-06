#include <iostream>
#include <fstream>
#include <list>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <sstream>
#include <iterator>

using namespace std;

void write_vtkHeader( FILE *fp, int * local_xlength )
{
    if( fp == NULL )
    {
        char szBuff[80];
        sprintf( szBuff, "Null pointer in write_vtkHeader" );
        return;
    }

    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output (written by Ye Tao) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i %i \n", local_xlength[0], local_xlength[1], local_xlength[2]);
    fprintf(fp,"POINTS %i float\n", (local_xlength[0]) * (local_xlength[1]) * (local_xlength[2]));
    fprintf(fp,"\n");
}


int main(int argc, char **argv){

	int num_atttr = 7;

	if(argc != 7){
		cerr << "Usage: <./exec> <input_dir/filename> <output_dir/> <number output files> <xCells> <yCells> <zCells>" << endl;
		return -1;
	}

	int timesteps = atoi(argv[3]);
	int local_xlength[3];
	local_xlength[0] = atoi(argv[3]);
	local_xlength[1] = atoi(argv[4]);
	local_xlength[2] = atoi(argv[6]);

	int tot_cells = local_xlength[0] * local_xlength[1] * local_xlength[2];
	string prefix(argv[2]); 
	prefix = prefix + "drivenCavity.";

	FILE *fp = fopen(argv[1], "rb");
	double buffer[7];
	
	// read file
	/*	FILE *fp;
	fp = fopen(argv[1], "rb");	
	double buffer[7];
	FILE *ofp = fopen("a.csv", "w");
	int i =0;
        while(i < 27000){
                i++;
                fread(buffer, sizeof(double), 7, fp);
                fprintf(ofp, "%f %f %f %f\n", buffer[0], buffer[1], buffer[2], buffer[6]);
		}*/

	
	list< list<double> > dataset;

	// read file
	for(int i = 0; i < timesteps * tot_cells; i++ ){
		list<double> ls;
		fread(buffer, sizeof(double), 7, fp);

		for(int j = 0; j < 7; j++) ls.push_back( buffer[j] );

		dataset.push_back(ls);
	}
	fclose(fp);
	
	list< list<double> >::const_iterator it1 = dataset.begin();
	for(int i = 0; i < timesteps; i++){

	  stringstream ss;
	  ss << i;
	  string filename = prefix + ss.str() + ".csv";
	  
	  FILE *ofp = fopen(filename.c_str(), "w");
	  int j =0;
	  
	  while(j < tot_cells){
	        j++;
		list<double>::const_iterator it2 = (*it1).begin();
                fprintf(ofp, "%f %f %f %f %f %f %f\n", *next(it2,0), *next(it2,1), *next(it2,2), *next(it2,3), *next(it2,4), *next(it2,5), *next(it2,6));
		it1++;
	  }
	  fclose(ofp);
	}
	
	  //list< list<double> >::const_iterator data_it = dataset.begin();
	/*
	// write file
	for(int i = 0; i < timesteps; i++){
	  stringstream ss;
	  ss << i;
	  string filename = prefix + ss.str() + ".vtk";
	  
	  fp = fopen(filename.c_str(), "w");
	  write_vtkHeader( fp, local_xlength );
	  
	  for(int k = 0; k < tot_cells; k++){
	    list<double>::const_iterator ls_it = (*next(data_it, k)).begin();
	    fprintf(fp, "%f %f %f\n", *next(ls_it, 0), *next(ls_it, 1),  *next(ls_it,2));
	  }

	  fprintf(fp,"POINT_DATA %i \n", (local_xlength[0]) * (local_xlength[1]) * (local_xlength[2]));
cout<< "OK" << endl;
	  fprintf(fp,"\n");
	  fprintf(fp, "VECTORS velocity float\n");

	  for(int k = 0; k < tot_cells; k++){
	    list<double>::const_iterator ls_it = (*next(data_it, k)).begin();
	    fprintf(fp, "%f %f %f\n", *next(ls_it, 3), *next(ls_it, 4),  *next(ls_it,5));
	  }

	  fprintf(fp,"\n");

	  fprintf(fp, "SCALARS density float 1 \n");
	  fprintf(fp, "LOOKUP_TABLE default \n");
	  for(int k = 0; k < tot_cells; k++){
	    list<double>::const_iterator ls_it = (*next(data_it, k)).begin();
	    fprintf(fp, "%f\n", *next(ls_it, 6));
	  }
	  
	  fclose(fp);
	  
	  data_it = next(data_it, tot_cells);

	}
	*/
	return 0;
}
