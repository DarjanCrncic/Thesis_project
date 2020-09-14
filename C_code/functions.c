
#include <math.h>
#include <arm_math.h>
#include <arm_const_structs.h>
#include "EventRecorder.h" 

#define dataSize 2000
#define crosPointNumber 1000

#define maxdist 1024
#define maxClust 500
#define simMatrixSize 100
#define mergingConst 0.7f

extern float data_freq[dataSize];
extern float data_time[dataSize];

// matrix init.			
extern short int neighIndices[dataSize];
extern short int seeds[dataSize];
extern short int orderedList[dataSize];
extern float coreDistList[dataSize];
extern signed char procesList[dataSize];
extern float reachDistList[dataSize];
extern float orderedReachDistList[dataSize];

extern short int clusterIndices[dataSize];
extern short int prevVectorIndices[crosPointNumber];
extern short int currVectorIndices[crosPointNumber];
extern short int prevCluster[simMatrixSize];
extern short int currCluster[simMatrixSize];
extern short int newClustTest[simMatrixSize];
	
extern int setClustersStart[maxClust];
extern int setClustersEnd[maxClust];
	
extern float timeMin[simMatrixSize];
extern float timeMax[simMatrixSize];
extern float yMinEnd[simMatrixSize];
extern float yMaxEnd[simMatrixSize];
extern float yMinStart[simMatrixSize];
extern float yMaxStart[simMatrixSize];
extern float prevyMaxEnd[simMatrixSize];
extern float prevyMinEnd[simMatrixSize];
// sorting function - quicksort algorithm

void quicksort(float number[],int first,int last){
   
	int i, j, pivot, mem;
	float temp;
		
   if(first<last){
      pivot=first;
      i=first;
      j=last;

      while(i<j){
         while(number[i]<=number[pivot]&&i<last)
            i++;
         while(number[j]>number[pivot])
            j--;
         if(i<j){
            temp=number[i];
					  mem=neighIndices[i];
					 
            number[i]=number[j];
						neighIndices[i]=neighIndices[j];
					 
            number[j]=temp;
						neighIndices[j]=mem;
         }
      }

      temp=number[pivot];
			mem=neighIndices[pivot];
			
      number[pivot]=number[j];
			neighIndices[pivot]=neighIndices[j];
			
      number[j]=temp;
			neighIndices[j]=mem;
			
      quicksort(number,first,j-1);
      quicksort(number,j+1,last);

   }
}


// function for sorting seeds - insertion sort algorithm

void sortSeeds( int seedNum){
  //EventStartD(0);
	int i, mem, hole;
		
	for(i=1; i<seedNum; i++)
	{
		float value = reachDistList[seeds[i]];
		mem = seeds[i];
		
		hole = i;
		while( hole>0 && reachDistList[seeds[hole-1]]>value)
		{
			seeds[hole]=seeds[hole-1];
			hole--;
		}
		seeds[hole] =  mem;
		
	}
	//EventStopD(0);
}

// function that returns points neighbors indices and their total number

int getNeighbors (float eps, int pointIndex, float distances[], int dataSizeTest){
	
	//EventStartC(0);
	
	
	float result;
	int k = 0;
	int i;
	
	for(i=0;i<dataSizeTest;i++){
		
			if (i==pointIndex)
				continue;
			
			arm_sqrt_f32((data_freq[pointIndex]-data_freq[i])*(data_freq[pointIndex]-data_freq[i])+(data_time[pointIndex]-data_time[i])*(data_time[pointIndex]-data_time[i]), &result);
			if (result<=eps){
				neighIndices[k]=i;
				distances[k]=result;
				k++;				
			}					
	}

	//EventStopC(0);
	
return k;
}

// function that calculates points core distance

float calcCoreDist ( int pointIndex, int neighCount, int Nmin, float distances[]){
	
	//EventStartC(1);
	
	float coreDist;
	
	if (coreDistList[pointIndex]!=-1)
		coreDist = coreDistList[pointIndex];
	
	else {		
		if (neighCount >= Nmin-1){			
			//EventStartB(0);
			quicksort(distances, 0, neighCount-1);
			//EventStopB(0);
			coreDist = distances[Nmin-2];	
		}else
		coreDist=-1;
	}
	//EventStopC(1);
	return coreDist;
}

// function that returns the next unprocessed point

int getUnproces (int dataSizeTest){
	
	int indUnproc = -1;
	int i;
	
	for (i=0;i<=dataSizeTest;i++){
		if (procesList[i]==0){
			indUnproc=i;
			break;
		}
	}
	return indUnproc;
	
}


// function for updating seeds list and the seeds reachability distance

int update ( int neighCount, int seedNum, int pointIndex , int change, float distances[]){
	
	//EventStartC(2);
	
	int i, ind_neigh;
	float  new_reach;
	change=0;
	
	for(i=0;i<neighCount;i++){
		
		ind_neigh=neighIndices[i];
			
		if(procesList[ind_neigh]==0){
			
			if(coreDistList[pointIndex] > distances[i])
				new_reach=coreDistList[pointIndex];
			else			
				new_reach = distances[i];
						
			if(reachDistList[ind_neigh]==-1){
				
				reachDistList[ind_neigh] = new_reach;
				seeds[seedNum] = ind_neigh;
				seedNum = seedNum+1;
				change=1;
				
			}else{
				
			if(new_reach<reachDistList[ind_neigh]){
				reachDistList[ind_neigh] = new_reach;
				change=1;
			}
			}
		}
	}
	
	//EventStopC(2);
return seedNum;
}


// optics function
 
void optics ( float eps, int Nmin, int dataSizeTest ){
	
int neighCount=0, seedNum=0, unprocesCount=dataSizeTest, pointIndex, orderedCount=0, neighbor, change;
int i, neighCount_mark, seedNumPrev;
float distances[dataSize]={0};	

while(unprocesCount){
	
	pointIndex = getUnproces(dataSizeTest);
	
	neighCount = getNeighbors (  eps, pointIndex, distances, dataSizeTest);
	
	procesList[pointIndex] = 1;
	
	orderedList[orderedCount]=pointIndex;
	orderedCount+=1;
	unprocesCount-=1;
	
	coreDistList[pointIndex]=calcCoreDist ( pointIndex, neighCount, Nmin,  distances);
	
	if(coreDistList[pointIndex]!=-1){
			
			for(i=0;i<dataSizeTest;i++)
				seeds[i]=-1;
			
			seedNum=0;
			seedNum = update ( neighCount, seedNum, pointIndex, change, distances);
	
			while(seeds[0]!=-1){
				
						if (seedNumPrev!=seedNum || change)
							sortSeeds(  seedNum );
						
						neighbor = seeds[0];
						procesList[neighbor]=1;
						orderedList[orderedCount] = neighbor;
						orderedCount = orderedCount + 1;
            unprocesCount = unprocesCount - 1;
			
						for(i=1;i<dataSizeTest;i++)
							seeds[i-1]=seeds[i];
						seeds[dataSizeTest-1]=-1;
						seedNum-=1;
			
						neighCount_mark=getNeighbors( eps, neighbor, distances, dataSizeTest);
			
						coreDistList[neighbor]=calcCoreDist( neighbor, neighCount_mark, Nmin, distances);
						
						seedNumPrev=seedNum;
					
						if(coreDistList[neighbor]!=-1)
							seedNum = update ( neighCount_mark, seedNum, neighbor, change, distances);
	
					}
	}	
}	
}

// function for returning cluster index for every point from original dataset

void getClusterIndices( int clustNum, int dataSizeTest){
	
	int i, j, index;
	
	for(i=0;i<dataSizeTest;i++){
			index = orderedList[i];
			for(j=0;j<clustNum;j++){				
					if(i>=setClustersStart[j] && i<=setClustersEnd[j]){
						clusterIndices[index]=j+1;
						break;
					}
			}
	}
}

// function for merging clustering results from two windows

int opticsMerging( int clustNum,  int iter, int newClustIndex, int dataSizeTest, int crossPoint){

	int i,j;
	uint32_t index;
	float value, result=0, result2=0;
	int simMatrix[simMatrixSize][simMatrixSize]={0};
	float pom[simMatrixSize]={0}, pom2[simMatrixSize]={0};
	
	if(iter>0 && clustNum>0 ){
		
		for(i=0;i<crossPoint;i++){
			currVectorIndices[i]=clusterIndices[i];
			if( prevVectorIndices[i]!=0 && currVectorIndices[i]!=0){
            simMatrix[prevVectorIndices[i]][currVectorIndices[i]]+=1;    
      }
		}
		
		for(i=0;i<clustNum;i++){
		   
			 result = 0;
			 result2 = 0;
			 for(j=1;j<simMatrixSize;j++){
				 pom[j]=simMatrix[j][i+1];
			 }
			 arm_max_f32(pom,simMatrixSize, &value, &index);
			 
			 for(j=1;j<simMatrixSize;j++){
				 pom2[j]=simMatrix[index][j];
			 }
			 
			 for(j=0;j<simMatrixSize-1;j++){
				 result=result+pom[j];
				 result2=result2+pom2[j];
		   }
			 
			 if(value >= mergingConst*result && value >= mergingConst*result2 && value!=0){
				 currCluster[i]=prevCluster[index-1];
			 }else{
				 currCluster[i]=newClustIndex;
				 newClustIndex+=1;
				 newClustTest[i]=1;
			 }
			 
		}
			 			
	}else{
		for(i=0;i<clustNum;i++){
			newClustTest[i]=1;
			currCluster[i]=i+1;			
		}
		newClustIndex=clustNum+1;
		for(i=0;i<crossPoint;i++)   {
			currVectorIndices[i]=clusterIndices[dataSizeTest-crossPoint+i];
	  }
	}
	
	return newClustIndex;
}	

// function for extracting features of clusters

void featureExtraction( int clustNum, float timeStampEnd, float timeStampStart, int Nmin, int dataSizeTest){
	
	int i, index, counter=0, start=0, j;
	//float maxTpoint[simMatrixSize]={0}, minTpoint[simMatrixSize]={0};
	
	for(i=0;i<dataSizeTest;i++){
		index = clusterIndices[i];
			
		if(index != 0){
			
			if(timeMin[index]==0){
				timeMin[index]=data_time[i];
				//minTpoint[index]=i;
			}
			timeMax[index]=data_time[i];
			//maxTpoint[index]=i;
			
			if(data_freq[i]>yMaxEnd[index])
				yMaxEnd[index]=data_freq[i];
			if(data_freq[i]<yMinEnd[index] || yMinEnd[index]==0)
				yMinEnd[index]=data_freq[i];
						
		}
	}
	
	for(i=1;i<=clustNum;i++){
		
		// cluster starts in this window
		if(newClustTest[i-1]==1){
			
			yMinStart[i]=0;
			yMaxStart[i]=0;
			counter=0;
			start=0;
			
			while(counter<Nmin && start<dataSizeTest){
				index=clusterIndices[start];
				
				if(index==i){
					
					counter+=1;
					
					if(data_freq[start]<yMinStart[i] || yMinStart[i]==0)
						yMinStart[i]=data_freq[start];
					if(data_freq[start]>yMaxStart[i])
						yMaxStart[i]=data_freq[start];
				}
				start+=1;
			}
		}
		
		// cluster ends in this window
		if(timeMax[i]-timeStampStart < 0.81f*(timeStampEnd-timeStampStart)){
			
			yMinEnd[i]=0;
			yMaxEnd[i]=0;
			start=dataSizeTest-1;
			counter=0;
			
			while(counter<Nmin && start>=0){
				index=clusterIndices[start];
				
				if(index==i){
					
					counter+=1;
					
					if(data_freq[start]<yMinEnd[i] || yMinEnd[i]==0)
						yMinEnd[i]=data_freq[start];
					if(data_freq[start]>yMaxEnd[i])
						yMaxEnd[i]=data_freq[start];
				}
				start-=1;
			}

		}else{
			timeMax[i]=timeStampEnd;
		}
		
		// cluster continues to an older cluster
		if(newClustTest[i-1]==0){
			
			timeMin[i]=timeStampStart;
			
			if(timeMax[i]<timeStampStart)
				timeMax[i]=timeStampStart;
			
			index=currCluster[i-1];
			j=1;
			
			while(prevCluster[j-1]!=0){				
				if(prevCluster[j-1]==index){				
					yMinStart[i]=prevyMinEnd[j];
					yMaxStart[i]=prevyMaxEnd[j];	
					break;
				}
				j+=1;
			}
		}
		
	}
}
	
///////////////////////////////////////////////////////////////////////////////////////////////////
// gradient algorithm

// function, returns size of a vector

float vectAbs(float rx, float ry, float wi){
	float result;
	arm_sqrt_f32((ry-rx)*(ry-rx)+wi*wi, &result);
	return result;
}

// function for determening inflection index of a point

float inflectionIndex(float x_r, float y_r, float z_r, float wi){
	
	float prevVector, nextVector;
	
	prevVector = vectAbs(x_r, y_r, wi);
	nextVector = vectAbs(y_r, z_r, wi);

	return (-wi*wi + (x_r-y_r)*(z_r-y_r))/(prevVector*nextVector); 
}

// function for determening the gradiant of a point

float gradientDet(float x_r, float y_r, float z_r, float wi){
	
	return wi*(y_r-x_r) - wi*(z_r-y_r);
}

// gradient clustering algorithm

int gradientClustering( int Nmin, float t, float wi, float largeClustPerc, float mergePerc, int minMax, int dataSizeTest) {

	int clustNum=0, i, j, k, lastEnd, currClusterStart=0, currClusterEnd=0, startPointNum=1, tempClustStart, tempClustEnd, clustSize, wrong;
	int startPoints[maxClust]={0};
	
	t=cos(t*PI/180);
	lastEnd=dataSizeTest-1;
	
	for(i=0;i<maxClust;i++)
		startPoints[i]=-1;
	
	startPoints[0]=0;
	
	for(i=1;i<dataSizeTest-2;i++){
		
		if(inflectionIndex(orderedReachDistList[i-1],orderedReachDistList[i], orderedReachDistList[i+1],wi)>t){
			
			if(gradientDet(orderedReachDistList[i-1],orderedReachDistList[i], orderedReachDistList[i+1],wi)>0){
				
				if(currClusterEnd-currClusterStart+1>=Nmin){
					setClustersStart[clustNum]=currClusterStart;
					setClustersEnd[clustNum]=currClusterEnd;
					clustNum+=1;				
				}
				currClusterEnd=0;
				currClusterStart=0;
			
				if(startPoints[0]!=-1){
					if(orderedReachDistList[startPoints[startPointNum-1]] <= orderedReachDistList[i]){
						startPoints[startPointNum-1]=-1;
						startPointNum-=1;
					}
				}
				
				if(startPoints[0]!=-1){
					while(orderedReachDistList[startPoints[startPointNum-1]]<orderedReachDistList[i]){
						tempClustStart=startPoints[startPointNum-1];
						tempClustEnd=lastEnd;
						
						if(tempClustEnd-tempClustStart+1>=Nmin){
							setClustersEnd[clustNum]=tempClustEnd;
							setClustersStart[clustNum]=tempClustStart;
							clustNum+=1;
						}
						startPoints[startPointNum-1]=-1;
						startPointNum-=1;
					}
					
					tempClustStart=startPoints[startPointNum-1];
					tempClustEnd=lastEnd;
					
					if(tempClustEnd-tempClustStart+1>=Nmin){
							setClustersEnd[clustNum]=tempClustEnd;
							setClustersStart[clustNum]=tempClustStart;
							clustNum+=1;
					}
					
				}
				
				if(orderedReachDistList[i+1]<orderedReachDistList[i]){
					startPoints[startPointNum]=i;
					startPointNum+=1;
				}
			}
			else {
				
				if(orderedReachDistList[i+1]>orderedReachDistList[i]){
					lastEnd=i;
					if(startPointNum==0){
						currClusterStart=0;
						currClusterEnd=lastEnd;
					}	
					else{
						currClusterStart=startPoints[startPointNum-1];
						currClusterEnd=lastEnd;
					}
				
				}
			}	
	}
	
}	

while(startPoints[0]!=-1){
	
	currClusterStart=startPoints[startPointNum];
	currClusterEnd=dataSizeTest-1;
	
	if(orderedReachDistList[startPointNum-1]>orderedReachDistList[dataSizeTest-1] && (currClusterEnd-currClusterStart+1)>=Nmin){
		
		setClustersEnd[clustNum]=tempClustEnd;
		setClustersStart[clustNum]=tempClustStart;
		clustNum+=1;
				
	}
	startPoints[startPointNum-1]=-1;
	startPointNum-=1;	

}

// add case when there are 0 clusters found
if (clustNum==0){
	return 0;
}


if(clustNum>1){
	
// remove large clusters, larger than largeClusterPerc, removing points of undetermened reachability from clusters
	
	i=0;
	while(i<clustNum){
		clustSize=setClustersEnd[i]-setClustersStart[i];
	
		if(clustSize>=largeClustPerc*dataSizeTest){			
			for(j=i;j<clustNum-1;j++){
				setClustersStart[j]=setClustersStart[j+1];
				setClustersEnd[j]=setClustersEnd[j+1];	
			}
			setClustersStart[clustNum-1]=0;
			setClustersEnd[clustNum-1]=0;
			clustNum-=1;
			continue;
		}
		if(orderedReachDistList[setClustersStart[i]]==maxdist){
			setClustersStart[i]=setClustersStart[i]+1;
			continue;
		}
		if(orderedReachDistList[setClustersEnd[i]]==maxdist){
			setClustersEnd[i]=setClustersEnd[i]-1;
			continue;
		}
		i+=1;				
	}

// removing clusters that spread over two clusters divided by noise point
	
	i=0;
	while(i<clustNum){
		wrong=0;
		
		for(j=setClustersStart[i];j<=setClustersEnd[i];j++){
			if(orderedReachDistList[j]==maxdist){
				wrong=1;
				break;
			}
		}
		
		if(wrong==1){
			for(k=i;k<clustNum-1;k++){
				setClustersStart[k]=setClustersStart[k+1];
				setClustersEnd[k]=setClustersEnd[k+1];
			} 
				setClustersStart[clustNum-1]=0;
				setClustersEnd[clustNum-1]=0;
				clustNum-=1;
				continue;
		}
		i=i+1;
	}
	

// merge simmilar cluesters, depending on mergeperc

	// sort clusters by size
	for(i=1;i<clustNum;i++){
		clustSize=setClustersEnd[i]-setClustersStart[i];
		tempClustEnd=setClustersEnd[i];
		tempClustStart=setClustersStart[i];
		k=i;
		
		while(k>0 && setClustersEnd[k-1]-setClustersStart[k-1]<clustSize){
			setClustersStart[k]=setClustersStart[k-1];
			setClustersEnd[k]=setClustersEnd[k-1];
			k=k-1;			
		}
		
		setClustersEnd[k]=tempClustEnd;
		setClustersStart[k]=tempClustStart;		
	}
	
	i=0;
	while(i<clustNum-1){
		
			j=i+1;
			while(j<clustNum){
				if(setClustersEnd[i]>=setClustersEnd[j] && setClustersStart[i]<=setClustersStart[j]){ //cluster j is a part of cluster i
					
					clustSize=setClustersEnd[i]-setClustersStart[i];
					if((setClustersEnd[j]-setClustersStart[j])>=clustSize*mergePerc){										//if size of cluster i is big enough to merge
						for(k=j;k<clustNum-1;k++){
							setClustersStart[k]=setClustersStart[k+1];
							setClustersEnd[k]=setClustersEnd[k+1];
						} 
						setClustersStart[clustNum-1]=0;
						setClustersEnd[clustNum-1]=0;
						clustNum-=1;
						continue;
					}
				}
				j+=1;
		  }
			i+=1;
	}


// keeping only the largest clusters, could contain two or more smaller ones
	if(minMax==2){
		
		i=0;
		while(i<clustNum-1){
			j=i+1;			
			while(j<=clustNum-1){
							
				if(setClustersEnd[j]<=setClustersEnd[i] && setClustersStart[j]>=setClustersStart[i]){
					for(k=j;k<clustNum-1;k++){
									setClustersStart[k]=setClustersStart[k+1];
									setClustersEnd[k]=setClustersEnd[k+1];
					} 
					setClustersStart[clustNum-1]=0;
					setClustersEnd[clustNum-1]=0;
					clustNum-=1;
					continue;					
				}
				j=j+1;
			}
			i=i+1;
		}
	}
}

for(i=0;i<clustNum;i++){
	setClustersEnd[i]=setClustersEnd[i]-1;
}

return clustNum;
	
} 
