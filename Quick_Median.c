// N16064103, Chun Yii Liu, HW3 Quick Median
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define SHOW 1
#define DEBUG 0
int quickmedian(int *x, int left, int right);

int main()
{
	int *x, *y, s;
	int i, j, N;
	
	srand( time(NULL) );
	
	for(N=10;N<=10;N*=2)
	{
		x = (int *) malloc( N * sizeof(int) );
		y = (int *) malloc( N * sizeof(int) );
		
		for(i=0;i<N;++i)
		{
			y[i] = x[i] = rand() % N;
		}
		
		#if SHOW
		for(i=0;i<N;++i)
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		#endif
		
		printf("-----------------\n");
		quickmedian(y,0,N);
		#if SHOW
		for(i=0;i<N;++i)
		{
			printf("x[%d]=%d\n",i,y[i]);
		}
		#endif
		
		if((N%2) == 0)	// The average of the middle two
		{
			printf("Median = %f\n",(y[(N-2)/2]+y[N/2])/2.0);
		}
        else			// The middle value
        {
        	printf("Median = %d\n",y[(N-1)/2]);
		}
	}
}

int quickmedian(int *x, int left, int right)
{
	int i, j, k, N = right;
	
	int pivot, t;
	
	if(left < right-1)
	{
		pivot = x[left];
    	i = left+1;
    	j = right-1;
    	
    	while(1)
    	{
    		while(i < right && pivot >= x[i]) i++; // search rightward to find 1st  pivot < x[i]
      		while(j >  left && pivot <  x[j]) j--; // search leftward to find 1st  pivot >= x[j]
      		
			#if DEBUG 
			printf("%d %d %d\n", i,j,pivot);
			#endif
			
			if(i>=j) break;	// Need no swap
			t = x[i];		// if not -> swap x[i], x[j]
      		x[i] = x[j];
      		x[j] = t;
      		
      		#if DEBUG
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			system("pause");
			#endif
		}
		x[left] = x[j];		// Swap
        x[j] = pivot;
        
        #if DEBUG
        printf("i=%d,j=%d\n",i,j);
		for(k=left;k<right;++k)
		{
			printf("x[%d]=%d\n",k,x[k]);
		}
		system("pause");
        #endif
        
        if((N%2) == 0)	// If N is an even number
		{
			if(j <= (N-2)/2)		// Need only right side for median
				quickmedian(x, j+1, right);
			else if(j > (N-2)/2)	// Need only left side
				quickmedian(x, left, j);
		}
        else			// If N is an odd number
        {
        	if(j == (N-1)/2)		// Already found
        		return 1;
        	else if(j < (N-1)/2)	// Need only right side for median
        		quickmedian(x, j+1, right);
        	else if(j > (N-1)/2)	// Need only left side
        		quickmedian(x, left, j);
		}
	}
	else 
    {
    	return 1;
	}
}
