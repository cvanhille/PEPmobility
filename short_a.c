#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*ALLOW REPETITION OF PARTICLES IN ORDER OF MOVEMENT*/
/*NO ATTRACTION*/

void main(){

float A[4]={0.001,0.002,0.005,0.008};
float F[4]={0.07,0.2,0.5,0.8};

int a_i;
int f_i;

for (a_i=0; a_i<4; a_i++) {
for (f_i=0; f_i<4; f_i++) {

/*VARIABLES AND PARAMETERS*/
float alpha=A[a_i], fi=F[f_i];						/*Tumble probability and particle density (respectively)*/
int N=2000;								/*Number of slots*/
int M=fi*N;								/*Number of particles*/
int i;
srand(time(NULL));
int Tmax=1e6;								/*Maximum number of time steps*/
int measures=1e3;							/*Measure variable*/
int step;								/*Time step variable*/
int array[N];								/*Array of slots, with the occupied ones assigned to 1*/
int arraycluster[N];							/*Vector with the clusters on the above array*/
int C[N/2];								/*Correlations array*/
int Nclusters=0;
int positions[M];							/*Array of particles' positions*/	
int directions[M];							/*Array of particles' directions (1=+x; 0=-x;)*/
int cluster[M];								/*Vector of clusters the particles are in (0 for no cluster)*/
int Nsize[M];								/*Vector with the number of clusters of the same size (overall)*/
int Ndir[3];								/*Vector with the number of clusters in each direction : +, -, 0*/
for (i=0; i<M; i++) {Nsize[i]=0;}
int dist;
for (dist=0; dist<N/2; dist++) {C[dist]=0;}


/*FILE PARAMETERS*/
FILE *sc;
char buffersc[96];
sprintf(buffersc, "sizedistr_a%.3f_f%.3f_t%d_N%d.dat", alpha, fi, Tmax, N); /*Saves the number of clusters of the same size overall*/
FILE *dc;
char bufferdc[96];
sprintf(bufferdc, "dirdistr_a%.3f_f%.3f_t%d_N%d.dat", alpha, fi, Tmax, N); /*Saves the number of clusters of the same size overall*/
FILE *corr;
char buffercorr[96];
sprintf(buffercorr, "corr_a%.3f_f%.3f_t%d_N%d.dat", alpha, fi, Tmax, N); /*Saves the total number of clusters*/

/*REPETITION LOOP FOR GIVEN PARAMETERS*/
int m_i;
for (m_i=0; m_i<measures; m_i++) {

/*INITIAL STATE OF THE SYSTEM - RANDOM*/
for (i=0; i<M; i++) {positions[i]=0; cluster[i]=0;}	
for (i=0; i<M; i++) {						/*Particles are M, ordered from 1 to M*/
	int r1=rand()%N+1;					/*Random position (1 to N, not 0 to N-1, !!)*/
	int r2=rand()%2;					/*Random direction (1=+x; 0=-x;)*/
	int j, used;
	for (j=0; j<M; j++) {					/*Check if that position is occupied*/
		if (positions[j]==r1) {used=1; break;} 
		else {used=0;}
	}
	if (used==0) {positions[i]=r1; directions[i]=r2;}	/*Save the position if free*/
	else {i=i-1;}						/*Restart for same particle if occupied*/
}

/*DYNAMICS LOOP - TIME STEPS AND CHANGE OF POSITION AND DIRECTION: EVALUATE THE DYNAMICS*/
for (step=0; step<Tmax; step++) {					/*Loop for time steps - OPEN*/

		/*CHOOSE A RANDOM ORDER OF MOVEMENT*/
		int orderofmovement[M];					/*Choice of an order for moving the particles*/
		for (i=0; i<M; i++) {
			int par=rand()%M+1;				/*Random particle to move at turn i (1-M, not 0-m-1!!)*/
			orderofmovement[i]=par;				/*Save particle par in turn i (even if repetition happens)*/
		}

		/*CHANGE EACH PARTICLE STATE IN THE ABOVE ORDER*/
		int k;
		int c=1;
		for (k=0; k<M; k++) {					/*Movement of each particle*/

			int ptm=orderofmovement[k]-1;			/*Particle to move*/

			/*Allow for a change of direction*/
			int prob=rand();
			if (prob<=alpha*RAND_MAX) {
				int tmbl=rand();
				if (tmbl<=0.5*RAND_MAX) {
					if (directions[ptm]==0) {directions[ptm]=1;} 
					else {directions[ptm]=0;}
				}
			}

			/*Find the expected new position*/
			int npos;					/*npos is the new expected position of the particle ptm*/
			if (positions[ptm]==1) {
				if (directions[ptm]==0) {npos=N;} 
				else {npos=positions[ptm]+1;}
			}
			else if (positions[ptm]==N) {
				if (directions[ptm]==1) {npos=1;} 
				else {npos=positions[ptm]-1;}
			}
			else {
				if (directions[ptm]==0) {npos=positions[ptm]-1;} 
				else {npos=positions[ptm]+1;}
			}

			/*Find if the movement is allowed*/
			int j, used, ratt;
			for (j=0; j<M; j++) {				/*Check if that position is occupied*/
				if (positions[j]==npos) {used=1; break;} 
				else {used=0;}
			}
			if (used==0) {					/*Movement is allowed - position is free*/
				positions[ptm]=npos;			/*Change position*/
			}

		}

		/*Evaluate clusters*/
		for (i=0; i<N; i++) {array[i]=0; arraycluster[i]=0;}	/*Empty the array and the cluster sizes (reset)*/
		for (i=0; i<M; i++) {
			array[positions[i]-1]=i+1;
		}

		int d, p, pp, nn, j=0;
		for (i=0; i<N; i++) {if (array[i]==0) {d=i; break;}}	/*Displacement to avoid separating a cluster due to end of count*/
		for (i=0; i<N; i++) {
			p=i+d;						/*Position*/
			nn=p+1;						/*Next one*/
			pp=p-1;						/*Previous one*/
			if (p>=N) {p=p-N;}
			if (nn>=N) {nn=nn-N;}
			if (pp>=N) {pp=pp-N;}
			if (pp<0) {pp=pp+N;}
			if (array[p]!=0) {				/*Save cluster it is in for each position*/
				if (array[pp]!=0 || array[nn]!=0) {	/*Only if there is another particle before or after (cond. 4 cluster)*/
					arraycluster[p]=j+1;
				}
			}
			else if (array[p]==0 && arraycluster[pp]!=0) {j++;}	/*If not change cluster*/
		}
		Nclusters=j;
		if (arraycluster[p]!=0) {Nclusters++;}

		/*ASSIGN CLUSTER NUMBER TO PARTICLES (0 IF FREE)*/
		for (i=0; i<N; i++) {
			for (j=0; j<M; j++) {
				if (positions[j]-1==i) {cluster[j]=arraycluster[i];}
			}
		}

		/*CHANGE POSITION OF CLUSTERS*/
		int dir, size;
		float clusterprob;
		for (i=0; i<Nclusters; i++) {
			dir=0, size=0;
			for (p=0; p<M; p++) {				/*Find the overall direction of movement for the cluster*/
				if (cluster[p]==i+1) {
					size++;
					if (directions[p]==0) {dir=dir-1;}
					else {dir=dir+1;}
				}
			}

			if ((step+1)==Tmax) {Nsize[size-1]++;}	/*Save number of clusters of same size at the END*/
			
			if ((step+1)==Tmax) {
						if (dir<0) {Ndir[0]++;}
						else if (dir==0) {Ndir[1]++;}
						else if (dir>0) {Ndir[2]++;}
			}						/*Save number of clusters in the same direction at the END*/

			clusterprob=abs(dir)/(float)size;
			int rclust=rand();
			int cltop=0, clbot=N, condmovcl=1;
			if (rclust<=clusterprob*RAND_MAX) {		/*Probability for whole cluster met*/
			for (p=0; p<M; p++) {
				if (dir>0) {if (cluster[p]==i+1 && positions[p]>cltop) {cltop=positions[p];}}
				if (dir<0) {if (cluster[p]==i+1 && positions[p]<clbot) {clbot=positions[p];}}
			}
			int ncltop=cltop+1;
			if (ncltop>N) {ncltop=ncltop-N;}
			int pclbot=clbot-1;
			if (pclbot<=0) {pclbot=pclbot+N;}
			for (p=0; p<M; p++) {
				if (dir>0) {if (positions[p]==ncltop) {condmovcl=0;}}
				if (dir<0) {if (positions[p]==pclbot) {condmovcl=0;}}
			}
			if (dir>0 && condmovcl!=0) {			/*Move whole cluster up one position*/
				for (p=0; p<M; p++) {
					if (cluster[p]==i+1) {
						positions[p]=positions[p]+1;
						if (positions[p]>N) {positions[p]=positions[p]-N;}
					}
				}
			}
			else if (dir<0 && condmovcl!=0) {		/*Move whole cluster down one position*/
				for (p=0; p<M; p++) {
					if (cluster[p]==i+1) {
						positions[p]=positions[p]-1;
						if (positions[p]<=0) {positions[p]=positions[p]+N;}
					}
				}
			}
			}
		}

		char message[96];
		sprintf(message,"%.2f %% elapsed of repetition %d out of %d",(step+1)/(double)Tmax*100,m_i+1,measures);
		if (step==0) {printf("%s",message);}
		else {
			int mm;
			for (mm=0; mm<96; mm++) {printf("\b");}
			printf("%s",message);
		}
			
}					/*End of the time loop - CLOSED*/

		/*SPATIAL CORRELATIONS*/
		int tag;
		int dist;
		int rest_i;
		int rest_j;
		for (tag=0; tag<N; tag++) {
			if (array[tag]==0) {continue;}
			for (dist=1; dist<N/2; dist++) {
				rest_i=tag+dist;
				rest_j=tag-dist;
				if (rest_i>=N) {rest_i=rest_i-N;}
				if (rest_j<=0) {rest_j=rest_j+N;}
				if (array[rest_i]!=0) {C[dist-1]++;}
				if (array[rest_j]!=0) {C[dist-1]++;}
			}
		}

}					/*END OF REPETITION LOOP*/

	/*SAVE CLUSTER SIZE DISTRIBUTION (NSIZE)*/
	sc=fopen(buffersc, "wb");
	for (i=0; i<M; i++) {if (Nsize[i]!=0) {fprintf(sc, "%d	%d\n", i+1, Nsize[i]);}}
	fclose(sc);

	/*SAVE CLUSTER DIRECTION DISTRIBUTION (NDIR)*/
	dc=fopen(bufferdc, "wb");
	for (i=0; i<3; i++) {if (Ndir[i]!=0) {fprintf(dc, "%d	%d\n", i-1, Ndir[i]);}}
	fclose(dc);

	/*SAVE CORRELATIONS*/
	corr=fopen(buffercorr, "wb");
	for (dist=0; dist<N/2; dist++) {fprintf(corr, "%d	%d\n", i-1, C[dist]);}
	fclose(corr);
}

}		/*Close phi loop*/
}		/*Close alpha loop*/
