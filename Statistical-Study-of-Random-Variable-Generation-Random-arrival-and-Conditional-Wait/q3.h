#include <stdio.h>
#include<iostream>
#include<stdint.h>
#include<ctime>
#include<math.h>
#include<fstream>
using namespace std;
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti = N + 1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
	mt[0] = s & 0xffffffffUL;
	for (mti = 1; mti<N; mti++) {
		mt[mti] =
			(1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		mt[mti] &= 0xffffffffUL;

	}
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
	int i, j, k;
	init_genrand(19650218UL);
	i = 1; j = 0;
	k = (N>key_length ? N : key_length);
	for (; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
			+ init_key[j] + j; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++; j++;
		if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
		if (j >= key_length) j = 0;
	}
	for (k = N - 1; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
			- i; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++;
		if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
	}

	mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
	unsigned long y;
	static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N) { /* generate N words at one time */
		int kk;

		if (mti == N + 1)   /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* a default initial seed is used */

		for (kk = 0; kk<N - M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk<N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}
/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
	unsigned long a = genrand_int32() >> 5, b = genrand_int32() >> 6;
	return(a*67108864.0 + b)*(1.0 / 9007199254740992.0);
}

int main()
{
	double x,y;
	int* hist[2];
	int bin;
	double Num = 100000;
	double binSize = 0.1;
	double maxVal = 10.0;
	int nBins = maxVal/binSize;
	double t = 1/12;
    double s1 = 0.0, s2 = 0.0;
    double l = 0.333333, a = 0.333333;
    double b = 0.6666667;
	// Initialize bins to zero
	for (int j = 0; j < 2; j++)
	{
	    hist[j] = new int[nBins];
		for (int i = 0; i < nBins; i++)
		{
			hist[j][i] = 0.0;
		}
	}
	init_genrand((unsigned long)time(NULL));

	for (int i = 0; i < Num; i++)
	{
            x = -l * log(1 - genrand_res53());
            if (x >= t)
			{
            x = x - t;
            s1++;
			}
			bin = floor(x / binSize);
            if(bin<nBins)
            hist[0][bin]++;
            y = (pow((1.0 - genrand_res53()), (-a)) - 1)*b;
			while (y >= t)
			{
            y = y - t;
            s2++;
			}
			bin = floor(y / binSize);
            if(bin<nBins)
            hist[1][bin]++;
	}
    FILE *fp1, *fp2;
	fp1 = fopen("exp_bus.txt","w");
	fp2 = fopen("par_bus.txt","w");
	fprintf(fp1,"%f %g \n", 0.0, 1.0);
	fprintf(fp2,"%f %g \n", 0.0, 1.0);
	double count = 0.0000;
	for (int i = 0; i < nBins; i++)
		{

			count = count + hist[0][i]/s1;
			if((1-count) >= 1/s1)
            fprintf(fp1,"%f %g \n",(i+1)*binSize ,1-count);
		}
	count = 0.0000;
		for (int i = 0; i < nBins; i++)
		{
            count = count + hist[1][i]/s2;
			fprintf(fp2,"%f %g \n",(i+1)*binSize ,1-count);
		}
	fclose(fp1);
	fclose(fp2);
	return 0;
}
