// N16064103, ¼B«T¼Ý Jimmy Liu
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
#define DEBUG 0
#define SHOW 1		// To see the answers: 0->1

int main()
{
	int i, N;
	double *x_r, *x_i; // y = fft(x);
	clock_t t1, t2;
	
	// Input N 
	printf("N = ");
	scanf("%d", &N);
	printf("N = %d\n", N);
	
	x_r = (double *) malloc(N * sizeof(double));
	x_i = (double *) malloc(N * sizeof(double));
	
	x_r[0] = 0.0;		// Only need this for x_r
	for(i=0;i<N;++i)
	{
		x_i[i] = 0.0;
	}
	
	t1 = clock();
	bit_reverse(x_r, x_i, N);
	butterfly_235(x_r, x_i, N);
	t2 = clock();
	printf("Time = %f s\n",1.0*(t2-t1)/CLOCKS_PER_SEC);
	system("pause");
	
	#if SHOW
	for(i=0;i<N;++i)
	{
		printf("y[%d] : %f + %f i\n", i, x_r[i], x_i[i]);
	}
	#endif
	
	return 0;
}

int bit_reverse(double *x_re, double *x_im, int N)
{
	int i, j, N2, N_b, Nb_c, *bt, p, q, k;
	
	bt = (int *) malloc(N * sizeof(int));
	N2 = N;
	N_b  = 0;
	while(N2>1)
	{
		if(N2 % 2 == 0) { N2 /= 2; bt[N_b] = 2;}
		else if(N2 % 3 == 0) { N2 /= 3; bt[N_b] = 3;}
		else if(N2 % 5 == 0) { N2 /= 5; bt[N_b] = 5;}
		N_b++;
	}
	# if DEBUG
	printf("No. of bit = %d\n", N_b);
	for(i=0;i<N_b;++i)
	{
		printf("%d ", bt[i]);
	}
	#endif
	
	Nb_c = N_b;			// Set a constant
	q = N/bt[N_b-1];	// 4 (N = 12)
	for(p=1;p<N-1;++p)
    {
		x_re[p] = q;
        k = N/bt[Nb_c-1];	// Always need to be 4
        N_b = Nb_c;			// final 2 or 3 or 5
        while(q >= (bt[N_b-1]-1)*k)			// q >=k ²Ä (log_2 k + 1)¦ì¬O1,  
        {
            q = q - (bt[N_b-1]-1)*k;		// 1->0
            if(k % bt[N_b-1] != 0) {N_b--;}	// e.g. 2*2*3*3: next k still /3
            k = k/bt[N_b-1];				// To check next bit ->
        }
        q = q+k;
    }
    x_re[N-1] = q;
    
    return 0;
}

int butterfly_235(double *x_re, double *x_im, int N)
{
	int i, k, m, p, q, r, s, t, N_bt, N2, *bt;
	double t1, t2, t3, t4, t5, t6, t7, t8, a, b; // For clearer w_3
	double w_N_re, w_N_im, w_N_re2, w_N_im2, w_N_re3, w_N_im3, w_N_re4, w_N_im4;
	double w_re, w_im, w_re2, w_im2, w_re3, w_im3, w_re4, w_im4, *w5;
	
	bt = (int *) malloc(N * sizeof(int));
	N2 = N;
	N_bt  = 0;
	while(N2>1)
	{
		if(N2 % 2 == 0) { N2 /= 2; bt[N_bt] = 2;}
		else if(N2 % 3 == 0) { N2 /= 3; bt[N_bt] = 3;}
		else if(N2 % 5 == 0) { N2 /= 5; bt[N_bt] = 5;}
		N_bt++;
	}

	a = -0.5; b = -sqrt(3)/2;
	m = 1;
	for (i=N_bt-1;i>=0;--i)		// *The order of butterfly should be opposite to Bit-reverse
	{
		if (bt[i] == 2)
		{
			// N = 2^p
			w_re = 1.0;
			w_im = 0.0;
			w_N_re =  cos(M_PI/m);
			w_N_im = -sin(M_PI/m);
			for(k=0;k<m;++k)	// m = 1: run once, m=2: run twice ...
			{
				for(p=k;p<N;p+=2*m)
				{
					q = p + m;
					t1 = x_re[q]; 
					x_re[q] = w_re*x_re[q] - w_im*x_im[q];
					x_im[q] = w_re*x_im[q] + w_im*t1;
					
					t1 = x_re[p];
					x_re[p] = x_re[p] + x_re[q];	// x0 + x4, x2 + x6 ...
					x_re[q] = t1      - x_re[q]; 	// x0 - x4, x2 - x6 ...
					t1 = x_im[p];
					x_im[p] = x_im[p] + x_im[q];
					x_im[q] = t1      - x_im[q];
				}
				t1   = w_re;
				w_re = w_N_re*w_re - w_N_im*w_im;
				w_im = w_N_re*w_im + w_N_im*t1;
			}
		}
		else if (bt[i] == 3)
		{
			// N = 3^q
			w_re = 1.0;
			w_im = 0.0;
			w_re2 = 1.0;
			w_im2 = 0.0;
			w_N_re =  cos(2.0*M_PI/3/m);
			w_N_im = -sin(2.0*M_PI/3/m);
			w_N_re2 =  cos(4.0*M_PI/3/m);
			w_N_im2 = -sin(4.0*M_PI/3/m);
			for(k=0;k<m;++k)	// m = 1: run once, m=3: run trice ...
			{
				for(p=k;p<N;p+=3*m)
				{
					q = p + m;
					r = p + 2*m;
					t1 = x_re[q];
					x_re[q] = w_re*x_re[q] - w_im*x_im[q];
					x_im[q] = w_re*x_im[q] + w_im*t1;
					t2 = x_re[r];
					x_re[r] = w_re2*x_re[r] - w_im2*x_im[r];
					x_im[r] = w_re2*x_im[r] + w_im2*t2;
					
					t1 = x_re[p];
					t2 = x_re[q];
					t3 = x_re[r];
					x_re[p] = x_re[p] + x_re[q] + x_re[r];	//
					x_re[q] = t1 + (x_re[q]*a - x_im[q]*b) + (x_re[r]*a + x_im[r]*b);
					x_re[r] = t1 + (t2*a + x_im[q]*b)      + (x_re[r]*a - x_im[r]*b);
					t1 = x_im[p];
					t4 = x_im[q];
					x_im[p] = x_im[p] + x_im[q] + x_im[r];
					x_im[q] = t1 + (x_im[q]*a + t2*b) + (x_im[r]*a - t3*b);
					x_im[r] = t1 + (t4*a - t2*b)      + (x_im[r]*a + t3*b);
				}
				t1    = w_re;
				w_re = w_N_re*w_re - w_N_im*w_im;
				w_im = w_N_re*w_im + w_N_im*t1;
				t1    = w_re2;
				w_re2 = w_N_re2*w_re2 - w_N_im2*w_im2;
				w_im2 = w_N_re2*w_im2 + w_N_im2*t1;
			}
		}
		else if (bt[i] == 5)
		{
			// N = 5^r
			w_re = 1.0;
			w_im = 0.0;
			w_re2 = 1.0;
			w_im2 = 0.0;
			w_re3 = 1.0;
			w_im3 = 0.0;
			w_re4 = 1.0;
			w_im4 = 0.0;
			w5 = (double *) malloc(8 * sizeof(double));	// w5^1 ~ w_5^4 
			for (k=0;k<4;++k)
			{
				w5[2*k] = cos(2.0*(k+1)*M_PI/5);
				w5[2*k+1] = -sin(2.0*(k+1)*M_PI/5);
			}
			w_N_re =  cos(2.0*M_PI/5/m);
			w_N_im = -sin(2.0*M_PI/5/m);
			w_N_re2 =  cos(4.0*M_PI/5/m);
			w_N_im2 = -sin(4.0*M_PI/5/m);
			w_N_re3 =  cos(6.0*M_PI/5/m);
			w_N_im3 = -sin(6.0*M_PI/5/m);
			w_N_re4 =  cos(8.0*M_PI/5/m);
			w_N_im4 = -sin(8.0*M_PI/5/m);
			for(k=0;k<m;++k)	// m = 1: run once, m=5: run 5 times ...
			{
				for(p=k;p<N;p+=5*m)
				{
					q = p + m;
					r = p + 2*m;
					s = p + 3*m;
					t = p + 4*m;
					t1 = x_re[q];
					x_re[q] = w_re*x_re[q] - w_im*x_im[q];
					x_im[q] = w_re*x_im[q] + w_im*t1;
					t1 = x_re[r];
					x_re[r] = w_re2*x_re[r] - w_im2*x_im[r];
					x_im[r] = w_re2*x_im[r] + w_im2*t1;
					t1 = x_re[s];
					x_re[s] = w_re3*x_re[s] - w_im3*x_im[s];
					x_im[s] = w_re3*x_im[s] + w_im3*t1;
					t1 = x_re[t];
					x_re[t] = w_re4*x_re[t] - w_im4*x_im[t];
					x_im[t] = w_re4*x_im[t] + w_im4*t1;
					
					t1 = x_re[p];
					t2 = x_re[q];
					t3 = x_re[r];
					t4 = x_re[s];
					t5 = x_re[t];
					x_re[p] = x_re[p] + x_re[q] + x_re[r] + x_re[s] + x_re[t];
					x_re[q] = t1 + (w5[0]*x_re[q] - w5[1]*x_im[q]) + (w5[2]*x_re[r] - w5[3]*x_im[r]) + (w5[4]*x_re[s] - w5[5]*x_im[s]) + (w5[6]*x_re[t] - w5[7]*x_im[t]);
					x_re[r] = t1 + (w5[2]*t2 - w5[3]*x_im[q]) + (w5[6]*x_re[r] - w5[7]*x_im[r]) + (w5[0]*x_re[s] - w5[1]*x_im[s]) + (w5[4]*x_re[t] - w5[5]*x_im[t]);
					x_re[s] = t1 + (w5[4]*t2 - w5[5]*x_im[q]) + (w5[0]*t3 - w5[1]*x_im[r]) + (w5[6]*x_re[s] - w5[7]*x_im[s]) + (w5[2]*x_re[t] - w5[3]*x_im[t]);
					x_re[t] = t1 + (w5[6]*t2 - w5[7]*x_im[q]) + (w5[4]*t3 - w5[5]*x_im[r]) + (w5[2]*t4 - w5[3]*x_im[s]) + (w5[0]*x_re[t] - w5[1]*x_im[t]);
					t1 = x_im[p];
					t6 = x_im[q];
					t7 = x_im[r];
					t8 = x_im[s];
					x_im[p] = x_im[p] + x_im[q] + x_im[r] + x_im[s] + x_im[t];
					x_im[q] = t1 + (w5[0]*x_im[q] + w5[1]*t2) + (w5[2]*x_im[r] + w5[3]*t3) + (w5[4]*x_im[s] + w5[5]*t4) + (w5[6]*x_im[t] + w5[7]*t5);
					x_im[r] = t1 + (w5[2]*t6 + w5[3]*t2) + (w5[6]*x_im[r] + w5[7]*t3) + (w5[0]*x_im[s] + w5[1]*t4) + (w5[4]*x_im[t] + w5[5]*t5);
					x_im[s] = t1 + (w5[4]*t6 + w5[5]*t2) + (w5[0]*t7 + w5[1]*t3) + (w5[6]*x_im[s] + w5[7]*t4) + (w5[2]*x_im[t] + w5[3]*t5);
					x_im[t] = t1 + (w5[6]*t6 + w5[7]*t2) + (w5[4]*t7 + w5[5]*t3) + (w5[2]*t8 + w5[3]*t4) + (w5[0]*x_im[t] + w5[1]*t5);
				}
				t1    = w_re;
				w_re = w_N_re*w_re - w_N_im*w_im;
				w_im = w_N_re*w_im + w_N_im*t1;
				t1    = w_re2;
				w_re2 = w_N_re2*w_re2 - w_N_im2*w_im2;
				w_im2 = w_N_re2*w_im2 + w_N_im2*t1;
				t1    = w_re3;
				w_re3 = w_N_re3*w_re3 - w_N_im3*w_im3;
				w_im3 = w_N_re3*w_im3 + w_N_im3*t1;
				t1    = w_re4;
				w_re4 = w_N_re4*w_re4 - w_N_im4*w_im4;
				w_im4 = w_N_re4*w_im4 + w_N_im4*t1;
			}
		}
		m = m*bt[i];
	}
	
	return;
}
