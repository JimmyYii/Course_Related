// Chun Yii Liu, N16064103, HW2
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main()
{
	clock_t t1, t2;				// variables for computing clocks
	double **A, *x, *b, *c, T1;	// x = 1234, a location in memory, x[0] the value at memory 1234, 1235, 1236, ... 1041
	int i, j, k, L, M, N=20000;
	
	for(N=20000;N<=20000;N*=2)
	{
		A = (double **) malloc( N * sizeof(double*) ); // Take N double* memory,
													   // and put the initial pointer at A, A[0], A[1],...
													   // But the values of A[0], A[1], A[2],... A[N-1] are not given
													   
		A[0] = (double *) malloc( N*N*sizeof(double)); // Take N*N double* memory,													   
													   // and put the initial pointer at A[0]
													   // A[0][0], A[0][1],..., A[0][N*N-1] are not given yet
		
		#pragma omp parallel for
		for(i=1;i<N;++i) A[i] = A[0] + i*N;
		// A[1]=A[0]+N --> A[2]=A[0]+2N --> A[3]=A[0]+3N : Which is able to be parallelized
		
		x = (double *) malloc( N * sizeof(double) );
		b = (double *) malloc( N * sizeof(double) );
		c = (double *) malloc( N * sizeof(double) );
		
		//------------------------------------- Create random elements of matrix and vector
		
		M = N/4;
		//srand() below need a overall parallel here 
		#pragma omp parallel num_threads(4) private(i,j,k,L)	// If a parallel is here
		{														// parallels below is not needed and useless
			k = omp_get_thread_num();
			//printf("The seed at thread %d is : %d\n",k,time(NULL)>>k);	// >> : divided by 2
			srand(time(NULL)>>k);		// Setting initial value in each thread 
			//#pragma omp parallel for	// Then go on to take random number
			for(i=k*M;i<(k+1)*M;++i)
			{
				//L = omp_get_thread_num();
				//printf("thread %d, %d\n",k,L);
				for(j=0;j<N;++j)
				{
					A[i][j] = rand();
				}
				x[i] = rand();
			}
		}
		
		// ------------------------------------------------ Do Ax = b
		
		// For correct answer
		double t;
		t1 = clock();
		for(i=0;i<N;++i) 
		{
			t = 0.0;
			for(j=0;j<N;++j)
			{
				t += A[i][j]*x[j];
			}
			b[i] = t;	// Ax = b
		}
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Time of Ax = b (w/o parallelization): %f\n",T1);
		
		// --------------------------------------------------------- After parallelization
		t1 = clock();
		M = N/8;	// One thread do the work of 1/8
					// Fastening effect: 1/8 > 1/4 > 1/2 after testing
		#pragma omp parallel num_threads(8) private(i,j,k,t)
		{
			k = omp_get_thread_num();	// Prepare to use separate threads
			for(i=k*M;i<(k+1)*M;++i)
			{
				t = 0.0;
				L = omp_get_thread_num();
				//#pragma omp parallel for 
				for(j=0;j<N;++j)
				{
					//L = omp_get_thread_num();
					//printf("%d %d\n",k,L);
					t += A[i][j]*x[j];
				}
				c[i] = t;	// Ax = c
			}
		}
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Time of Ax = b (with parallelization): %f\n",T1);
		
		// Check correctness with answer without parallelization
		for(i=0;i<N;++i)
		{
			//printf("%f\n",fabs(b[i]-c[i]));
			if((b[i]-c[i]) != 0)
			{
				printf("Error in parallelization\n");
				break;
			}
		}
		
		free(b);
		free(c);
		free(x);
		free(A[0]);
		free(A);	
	} 
	return 0;
}
