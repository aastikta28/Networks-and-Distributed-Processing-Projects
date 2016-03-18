
#include <stdio.h>
#include<iostream>
#include<stdint.h>
#include<ctime>
#include<math.h>
#include<fstream>
#include<vector>
#include<algorithm>
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
	double on,off, sum=0.0;
	int* hist[1];
	int bin;
	double Num = 100000;
	double binSize = 1000;
	double maxVal = 10000000.0;
	int nBins = maxVal/binSize;
    vector<double> Y;
    vector<double> Z;
	//double t = 0.08333333;
  //  int s1 = 0, s2 = 0;
    double a = 3.0, b = 1.0, l=1.0;
	// Initialize bins to zero
    hist[0] = new int[nBins];
	for (int i = 0; i < nBins; i++)
		hist[0][i] = 0.0;
    double T = 100.0;
    double h = 0.16667;
    double gmin = 100000;
    int index1 = T, index2 = T-h;
   // cout << index1 << " " << index2 << endl;
	init_genrand((unsigned long)time(NULL));

	  while(sum <= T*Num)
        {
            Y.push_back(sum);
            on = (pow((1.0 - genrand_res53()), -1/a) - 1) * Num;
            off = -1 * Num * log(1.0 - genrand_res53());
            sum += (on+off);
        }
        sort(Y.begin(), Y.end());
       for(int i=0; i<Y.size(); i++)
        {
            bin = floor(Y[i] / binSize);
         //   cout << bin << endl;
            if(bin<nBins)
            hist[0][bin]++;
        }

	FILE *fp1, *fp2;
	fp1 = fopen("superposition3.txt","w");
	fprintf(fp1,"%f %g \n", 0.0, 1.0);
	double count = 0.0000;
	for (int i = 0; i < nBins; i++)
		{
		    cout<< hist[0][i] << endl;
		    count = count + hist[0][i]/(double)Y.size();//(double)s1);
			if((1-count) >= 1/(double)Y.size())
            fprintf(fp1,"%f %g \n",(i+1)*binSize ,1-count);
		}
    fclose(fp1);
   // Y.clear();
   /* int var;
    for(var=0; var<Num; var++)
    {
        sum = 0;
        while(sum <= T)
        {
            Y.push_back(sum);
            on = (pow((1.0 - genrand_res53()), -1/a) - 1);
            off = -1 * log(1.0 - genrand_res53());
            sum += (on+off);
        }

        for(var=0; var<Y.size(); var++)
        {
            double val = Y[var+1] - Y[var];
            if(val<gmin)
                gmin=val;
        }

        cout << gmin << endl;

        bin = floor(gmin / binSize);
        if(bin<nBins)
        hist[0][bin]++;

        Y.erase (Y.begin(),Y.end());
    } */

  /*  for(var=0; var<Z.size(); var++)
        {
            bin = floor(Z[var] / binSize);
            if(bin<nBins)
            hist[0][bin]++;
        } */

  /*   FILE* fp2 = fopen("superposition2.txt","w");
	fprintf(fp2,"%f %g \n", 0.0, 1.0);
    double count = 0.0000;
	for (var = 0; var < nBins; var++)
		{
		//    cout<< hist[0][var] << endl;
		    count = count + hist[0][var]/(double)Z.size();//(double)s1);
			if((1-count) >= 1/(double)Z.size())
            fprintf(fp2,"%f %g \n",(var+1)*binSize ,1-count);
		}
    fclose(fp2); */


	return 0;
}


