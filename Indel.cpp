#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>

#define N	2000
#define S	10000
#define SampleSize 6
//#define RandomConversionRate

using namespace std;

char genomes[2][N][S][2];
double u = 0.1;
double EPS = -0.001;
int FixedN, LostN;
int freq[S];
double real();	// A random real number between 0 and 1 is returned.
int TotalMutN;
int f[2*N], P1[2*N+1];

int main(){
	int time, n, i, j, s, k;
	char ch;
	cout<<"Input a different character to generate different random numbers: "; 
	cin>>ch;
	for (i=0; i<ch; i++)rand();
//	cout<<"RAND_MAX="<<RAND_MAX<<endl;
#ifdef RandomConversionRate
	cout<<"RandomConversionRate	";
#endif
	cout<<"N = "<<N<<"	u = "<<u<<"	eps = "<<EPS<<endl;
	int MaxTime = N*100;
	for (time=0; time<MaxTime; time++){
		int mutsite=0, mutated=0;
		int t0 = time%2, t1 = !t0;
		for (n=0; n<N*2; n++){
			float r=real();
			if (r<u){					// Mutation.
				while( freq[mutsite] && mutsite<S ) mutsite++;
				if (mutsite>=S) {
					cout<<"Total site number is not big enough or mutation rate is too high!"<<endl;
					cin>>n; 
					return 1 ;
				}
				genomes[t0][n/2][mutsite][n%2]=1;
				freq[mutsite]++;
				TotalMutN++;
				mutated=1;
			}
			for (s=0; s<S; s++){		// Gene conversion in heterozygotes.
				char * gen = genomes[t0][n/2][s];
				if (gen[n%2] && !gen[!(n%2)] ){
					double eps = EPS;
#ifdef RandomConversionRate
					eps += (real()-0.5)*0.01;
#endif
					if(eps>0 && real()<eps) gen[!(n%2)]=1;
					if(eps<0 && real()<-eps) gen[(n%2)]=0;
				}
			}
		}
		for (s=0; s<S; s++){			// Reproduction.
			int fr = freq[s];
			freq[s]=0;
			for (n=0; n<2*N; n++){
				char g=0;
				if (fr) {
					int nn = rand() % (2*N); 
					g = genomes[t0][nn/2][s][nn%2];
				}
				genomes[t1][n/2][s][n%2] = g;
				if (g) freq[s]++;
			}
			if (mutated && s==mutsite) P1[freq[mutsite]]++;
			if (freq[s]==2*N){				// Remove fixed genes in the population.
				freq[s]=0;
				for (n=0; n<2*N; n++) genomes[t1][n/2][s][n%2] = 0;
				FixedN++;
			}
			if (freq[s]) f[freq[s]]++;
		}
		if (time==MaxTime-1){
			time=time;
		}
		if (time%N==0){
//			cout<<"time = "<<time<<"..."<<endl;
			if (time>3*N) {
				int sample[SampleSize];
				for (i=0; i<SampleSize; i++){
					int SampleNeeded = 1;
					while (SampleNeeded){
						sample[i] = rand()%(2*N);
						SampleNeeded=0;
						for (j=0; j<i; j++) 
							if (sample[j]==sample[i]) SampleNeeded=1;
					}
				}
				int K[S];
				for (s=0; s<S; s++){
					K[s]=0;
					for (i=0; i<SampleSize; i++)
						K[s] += genomes[t0][sample[i]/2][s][sample[i]%2];
				}
				int Nk[SampleSize];
				for (k=0; k<SampleSize; k++){
					Nk[k]=0;
					for (s=0; s<S; s++){
						if (K[s]==k+1) Nk[k]++;
					}
					cout<<Nk[k]<<'\t';
				}
				cout<<endl;
			}
		}
	}

	return 0;
}



double real(){
	return double(rand())/RAND_MAX;
}
