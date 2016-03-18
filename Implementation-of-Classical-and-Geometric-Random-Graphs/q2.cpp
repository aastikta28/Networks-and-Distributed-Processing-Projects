
#include <stdio.h>
#include<iostream>
#include<stdint.h>
#include<ctime>
#include<math.h>
#include<fstream>
#include<stdlib.h>
#include<list>
#include<queue>
#include<vector>
#include<algorithm>
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


class Graph
{
    private:
        int V;
        list<int> *adj;
        void DFSUtil(int v, bool visited[]);
    public:
        Graph(int V)
        {
            this->V = V;
            adj = new list<int>[V];
        }
        ~Graph()
        {
            delete [] adj;
        }
        void addEdge(int v, int w);
        bool isConnected();
        Graph getTranspose();
        void clearGraph();
};

/*
 *  A recursive function to print DFS starting from v
 */
void Graph::DFSUtil(int v, bool visited[])
{
    visited[v] = true;
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            DFSUtil(*i, visited);
}

/*
 * Function that returns reverse (or transpose) of this graph
 */
Graph Graph::getTranspose()
{
    Graph g(V);
    for (int v = 0; v < V; v++)
    {
        list<int>::iterator i;
        for(i = adj[v].begin(); i != adj[v].end(); ++i)
        {
            g.adj[*i].push_back(v);
        }
    }
    return g;
}
/*
 * Add Edge to connect v and w
 */
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    //adj[w].push_back(v);
}

/*
 * Check if Graph is Connected
 */
bool Graph::isConnected()
{
    bool visited[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;
    DFSUtil(0, visited);
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
             return false;
    Graph gr = getTranspose();
    for(int i = 0; i < V; i++)
        visited[i] = false;
    gr.DFSUtil(0, visited);
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
              return false;
    return true;
}

/*
 * Main Contains Menu
 */
int main()
{
    // Generating Random number and subsequent exponential and pareto variables
	init_genrand((unsigned long)time(NULL));

    int Num = 10000;

        int n = 100;
        double c = log(n);
        double p = c/n;
        int Davg;
        double r = 0.1;
      //  int tot = floor(n*(n-1)/2);
      //  int act =  floor(n*(n-1)*p/2);

    int cnt = 0;
    for(int i=0; i<Num; i++)
    {

        Graph g1(n);
        Davg = (int) (n-1) * 22 / 7 * r * r;
        for(int k=0; k<n; k++)
        {
        for(int j=0; j<Davg; j++)
        {
            int rnd1 = (int)n*genrand_res53();
            g1.addEdge(k, rnd1);
        }
        }
        if (g1.isConnected())
        {
            cout<<"The Graph is Connected"<<endl;
            cnt++;
        }
        else
            cout<<"The Graph is not Connected"<<endl;



    }
    cout << "Count of connected graphs is: " << cnt << endl;

       return 0;
}
