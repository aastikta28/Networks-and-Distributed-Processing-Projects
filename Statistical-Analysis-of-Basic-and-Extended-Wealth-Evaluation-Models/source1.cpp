#include <stdio.h>
#include<iostream>
#include<stdint.h>
#include<ctime>
#include<math.h>
#include<fstream>
# define NUM (1<<24)
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
	double x=1.0;
	int* hist[1];
	double* val[1];
	int bin;
	double Num = 100000;
	double n;
	double binSize = 100;
	double maxVal = 10000.0;
	int nBins = maxVal/binSize;
    val[0] = new double[(int)Num];
	// Initialize bins to zero
	    hist[0] = new int[nBins];
		for (int i = 0; i < nBins; i++)
			hist[0][i] = 0.0;

	// Generating Random number and subsequent exponential and pareto variables
	init_genrand((unsigned long)time(NULL));
double max_wealth=0, avg_wealth=0;
double md;
	//printf("n is %f \n", n);
for(int j = 0; j < Num; j++)
{
    n = floor(10+ 30*genrand_res53());
    //printf("n is %f \n", n);
    x = 1.0;
    for (int i = 0; i < n; i++)
	{
		x = x*(1.0 + (-0.3 + genrand_res53()));
    }
    val[0][j] = x;

    if(x > max_wealth)
        max_wealth = x;
    avg_wealth=avg_wealth+x;
        bin = floor(x/binSize);
        if(bin<nBins)
        hist[0][bin]++;
}
for(int k=0;k<Num;k++)
{
    for(int a=k;a<Num;a++)
    {
        if(val[0][k] < val[0][a])
        {
            md = val[0][k];
            val[0][k] = val[0][a];
            val[0][a] = md;
        }
    }
    printf("%f \n",val[0][k]);
}
int mid=Num/100;
double upper_one=val[0][mid];
mid=99*Num/100;
double below_one=val[0][mid];
/*for(int a1=0;a1<(Num/100);a1++)
    below_one=below_one+val[0][a1];
for(int a2=0;a2<(Num*99/100);a2++)
    upper_one=upper_one+val[0][a2];*/
    double delta=upper_one/below_one;
avg_wealth=avg_wealth/Num;
double gamma = max_wealth/avg_wealth;
	/*FILE *fp1;
	fp1 = fopen("wealth_new.txt","w");
	double count = 0.0;
	for (int i = 0; i < nBins; i++)
		{
            //if((hist[0][i]/Num) >= 1/Num)
            fprintf(fp1,"%f %g \n",(i+1)*binSize ,hist[0][i]/Num);
		}
	fclose(fp1);*/

    printf("Ratio(gamma) of %f and %f is %f \n",max_wealth,avg_wealth, gamma);
    printf("Ratio(delta) of %f and %f is %f \n",upper_one,below_one, delta);
	//printf("Ratio is %f \n", gamma);
	return 0;
}

