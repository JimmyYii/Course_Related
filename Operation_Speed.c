// N16064103 Chun-Yii Liu, HW1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>	// To use sin()
#include <time.h>

int main()
{
	clock_t t1, t2;						// t1 and t2 is in ms
	int i, j=1, k=2, N=1000000000;		// N is a billion
	double T0, T1, a=1.234, b=2.345;
	printf("For computational time of 1000000000 (a billion) '+', '*' and 'sin', \n");
	
	// Do the empty loop first
	t1 = clock();
	for(i=0;i<1000000000;++i)
	{
	}
	t2 = clock();
	T0 = (t2-t1)/(double) CLOCKS_PER_SEC;	// CLOCKS_PER_SEC = 1000, T0 is in sec
	printf("Empty loop: %f\n\n",T0);
	
	// For "+" operation ---------------------------------------------------------------
	t1 = clock();
	for(i=0;i<1000000000;++i)
	{
		j = j+1;		// j=j+1: make sure everytime j changed
		//a = a+1;		// Try float number
		//a = a+1.234;	// more complicated
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("With one addition: %f\n",T1);
	printf("Real (+) time: %f\n\n",T1-T0);
	//printf("j = %d\n",j);	// For valid optimization
	
	// For "*" operation --------------------------
	t1 = clock();
	for(i=0;i<1000000000;++i)
	{
		k = k*2;		// k = k*2: make sure everytime k changed
		//b = b*2;		// Try float number
		//b = b*2.345;	// more complicated
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("With one multiply: %f\n",T1);
	printf("Real (*) time: %f\n\n",T1-T0);
	//printf("k = %d\n",k);	// For valid optimization
	
	// For "sin" operation ----------------------
	t1 = clock();
	for(i=0;i<1000000000;++i)
	{
		a = sin(a);		// time = 0 for highest optimization, if a = sin(a) highest optimization is valid
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("With one sin: %f\n",T1);
	printf("Real (sin) time: %f\n",T1-T0);
	//printf("a = %f\n",a);	// For valid optimization
}
