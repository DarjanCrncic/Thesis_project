#ifndef FUNC_H
#define FUNC_H


void quicksort(float number[],int first,int last);

int getNeighbors ( float eps, int pointIndex, float distances[], int dataSizeTest);

float calcCoreDist ( int pointIndx, int neighCount, int Nmin, float distances[]);

void sortDistances(float distances[], int size);

int getUnproces (int dataSizeTest);

int update ( int neighCount, int seedNum, int pointIndx,  int change, float distances[]);

void sortSeeds( int seedNum);

void optics ( float eps, int Nmin, int dataSizeTest );



float gradientDet(float x_r, float y_r, float z_r, float w);

float inflectionIndex(float x_r, float y_r, float z_r, float w);

float vectAbs(float rx, float ry, float w);

int gradientClustering( int Nmin, float t, float w, float largeClustPerc, float mergePerc, int minMax, int dataSizeTest);

void getClusterIndices(int clustNum, int dataSizeTest);



int opticsMerging( int clustNum,  int iter, int newClustIndex, int dataSizeTest, int crossPoint);

void featureExtraction( int clustNum, float timeStampEnd, float timeStampStart, int Nmin, int dataSizeTest );
	
#endif
