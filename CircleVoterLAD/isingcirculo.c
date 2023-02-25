/*************************************************************************
*                     Votante Latoski-Arenzon-Dantas 2D                  *
*                             V1.07 18/06/2021                           *
*************************************************************************/

/***************************************************************
 *                    OBRIGATORY DEFINITIONS
 **************************************************************/
// -D L="LATTICE_SIZE"
// -D DETA="ETA_INCREMENT"

/***************************************************************
 *                     OPTIONAL DEFINITIONS
 **************************************************************/
// -DNBINARY [non binary case]
// -DRESET  [full reset case]
// -DGRESET [gamma reset case]
// -DINTRANS [intrans case]

// -DSPEEDTEST [sweep and measure speed test]
// -DDEBUG [debug program]
// -DVISUAL [live gif of the evolution]
// -DSNAPSHOTS -I ~/VotanteLAD/liblat2eps/ -llat2eps [snapshots of the system]

/***************************************************************
 *                            INCLUDES                      
 **************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#ifdef SNAPSHOTS
  #include <lat2eps.h>
#endif
#include "mc.h"

/****************************************************************
 *                       PARAMETERS DEFINITIONS                      
 ***************************************************************/

#define N          (L*L)  //Lattice volume
#define MCS         1000000 //Max evolution time
#define THRESHOLD   1. //Certainty's treshold
#define MEASURES    40
#define ALPHA       1.  //Transiten probability 1
#define BETA        1.  //Transiten probability 2
#define GAMMA       1.05  //Gamma reset scale

/****************************************************************
 *                            SETTINGS 
 ***************************************************************/

#ifdef RESET
  #define RESET       1 
#endif
#ifdef GRESET
  #define RESET       2
#endif

#define LOGSCALE    1 // 0 --> measures logaritmically spaced, 1 --> measures in logscale.
#define SIMPLIFIED  1 // 1 --> simplify the algorithm to alpha=beta=1 to avoid calculations

/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/

void initialize(void); 
void openfiles(void); 
void sweep(void); 
void visualize(int,unsigned long); 
void states(void); 
void measures1(void); 
void measures2(void);
#ifdef SNAPSHOTS
  void snap(void);  
#endif
void hoshen_kopelman(void);
void connections(int,int);
int percolates2d(int);
bool exists(const char*);
bool probcheck(double);

/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

FILE *fp1,*fp2;
int *spin,**neigh,*memory,*measures,*zealot,*right,*left,*up, *down, sum, sumz, activesum;
int *siz, *label, *his, *qt, cl1, numc, mx1, mx2;
int probperc0,probperc1;
unsigned long seed;
double *certainty;

/***************************************************************
 *                          MAIN PROGRAM  
 **************************************************************/
int main(void){

  #if(SEED==0)
    seed = time(0);
  #else
    seed = SEED;
  #endif
  
  #if(SNAPSHOTS==0 && VISUAL==0)
    openfiles(); 
    
  #else
    #if(DEBUG==1)
      seed = 1111111111;
    #endif
  #endif

  int k=0;
  initialize();

  for (int j=0;j<=MCS+1;j++)  {
    #if(VISUAL==1)
      visualize(j,seed);
      sweep();
    #else
      if( ( qt[0]==0 ) | ( qt[1]==0 ) ){
        states();
        hoshen_kopelman();
        fprintf(fp1,"%d %.8f %.8f %.8f %.8f %.8f %d %.8f %d\n",j,(double)sum/N,(double)sumz/N,(double)activesum/N,(double)numc/N,(double)mx1/N,probperc0,(double)mx2/N,probperc1);
        while(measures[k]!=0){
          fprintf(fp1,"%d %.8f %.8f %.8f %.8f %.8f %d %.8f %d\n",measures[k],(double)sum/N,(double)sumz/N,(double)activesum/N,(double)numc/N,(double)mx1/N,probperc0,(double)mx2/N,probperc1);
          k++;
        }           
        break;
      }
      if (measures[k]==j) {  
        #if(SNAPSHOTS==1)
          snap();   
          k++;   
        #else  
          #if(SPEEDTEST==1)
            clock_t t=clock();
          #endif
          states();
          hoshen_kopelman();
          fprintf(fp1,"%d %.8f %.8f %.8f %.8f %.8f %d %.8f %d\n",j,(double)sum/N,(double)sumz/N,(double)activesum/N,(double)numc/N,(double)mx1/N,probperc0,(double)mx2/N,probperc1);
          k++;
          #if(SPEEDTEST==1)
            t=clock()-t;
            double time_taken = ((double)t)/CLOCKS_PER_SEC;
            printf("one measure: %fs\n",time_taken);
          #endif
        #endif
      }
      #if(SPEEDTEST==1)
        clock_t t=clock();
        sweep();
        t=clock()-t;
        double time_taken = ((double)t)/CLOCKS_PER_SEC;
        printf("onesweep: %fs\n",time_taken);
      #else
        sweep();
      #endif
    #endif
  }

  #if(SNAPSHOTS==0)
  fclose(fp1);
  #endif

}
/***************************************************************
 *                        INICIALIZAÇÃO  
 **************************************************************/
void initialize(void) {
 
  start_randomic(seed);

  his = malloc(N*sizeof(int));
  spin = malloc(N*sizeof(int));
  neigh = (int**)malloc(N*sizeof(int*));
  memory = malloc(N*sizeof(int));
  measures = malloc(MCS*sizeof(int));
  zealot = malloc(N*sizeof(int));
  right = malloc(N*sizeof(int));
  left = malloc(N*sizeof(int));
  up = malloc(N*sizeof(int));
  down = malloc(N*sizeof(int));
  certainty = malloc(N*sizeof(double));

    #if(NBINARY==0)
      qt = malloc(2*sizeof(int));
    #else
      qt = malloc(N*sizeof(int));
    #endif

  for(int i=0; i<N; i++){
    neigh[i] = (int*)malloc(4*sizeof(int));
    #if(NBINARY==1)
      qt[n] = 0;
    #endif
  }
  int k;
  #if(WALL==0)
    for(int n=0; n<N; n++) {
      certainty[n] = 0;
      zealot[n] = 0;
      memory[n] = 0;
      double a = pow(abs(n%L - L/2),2);
      double b = pow(abs(n/L - L/2),2);
      if( ( a + b ) <= N/16 ) {
        k=1;
      }
      else k=0;
      spin[n] = k*2 - 1; 
      qt[k]++;
    } 
  #else
    for(int n=0; n<N; n++) {
      certainty[n] = 0;
      zealot[n] = 1;
      memory[n] = 0;
      if( n < 2*L) {
        k=0;
      }
      else if( n >= N/2+8*L) {
        k=0;
      }
      else k=1;
      spin[n] = k*2 - 1; 
      qt[k]++;
    } 
  #endif
  for (int i = 0; i < N; i++) {    
    neigh[i][0] = (i+1)%L + (i/L)*L; //right
    neigh[i][1] = (i-1+L)%L + (i/L)*L; //left
    neigh[i][2] = (i-L+N)%N; //up
    neigh[i][3] = (i+L)%N; //down
    right[i] = neigh[i][0];
    left[i] = neigh[i][1];
    up[i] = neigh[i][2];
    down[i] = neigh[i][3];
  }

  #if(LOGSCALE==1)
    measures1(); 
  #else
    measures2();
  #endif

}

/****************************************************************
 *               MCS routine
 ***************************************************************/
void sweep(void) {

  for (int n=0; n<N; n++) {
    int site = FRANDOM*N;
    int E1=0;
    int E2=0;
    for(int i=0; i<4; i++){
      int neighbour = neigh[site][i];
      if(spin[neighbour]==spin[site])E1++;
      else E2++; 
    }
    if(E2>=E1)spin[site]=-spin[site];
  }
}

/****************************************************************
 *               Check states numbers
 ***************************************************************/
void states(void) {
  #if(NBINARY==0)  
    sum=N;
  #else
    sum=0;
  #endif
  sumz=0;
  activesum=0;
  for (int i=0; i<N; i++) {
    #if(NBINARY==0)
      if (memory[i]!=0) sum--;
    #else
      if (memory[i]>=0) sum++;
    #endif
    
    if (zealot[i]!=0) sumz++;
    if (spin[right[i]]!=spin[i]) activesum++;
    if (spin[down[i]]!=spin[i]) activesum++;
  }
}
 /**************************************************************
 *                       Measures Vector     // Medidas Logarítmicas
 *************************************************************/
void measures1(void){
  int n=0;
  int m=0;
  for(int i=0; i<=MCS; i++) {
    if(i%((int)pow(10,n))==0) {
      measures[m]=i;
      m++;
    }
    if(i/(pow(10,n))==10)n++;
  }
  return;
}

 /**************************************************************
 *                       Measures Vector    // Número de medidas
 *************************************************************/
void measures2(void){
  int m=0;
  measures[0]=0;
  double temp = pow((double)MCS,1.0/(MEASURES-1));
  for(int i=0;i<MEASURES;i++){
    measures[i] = (int)pow(temp,(double) i);
    if(measures[i]<=m)measures[i]=m+1;
    m=measures[i];
  }
  return;
}

/**************************************************************
 *                      Teste
 *************************************************************/

bool probcheck(double _ALPHA) {
  double r = FRANDOM;
  if(r<=_ALPHA)return true;
  else return false;
}

/**************************************************************
 *                       Vizualização                   
 *************************************************************/
void visualize(int _j,unsigned long _seed) {
  int l;
  printf("pl '-' matrix w image t 'time = %d seed = %ld'\n",_j,_seed);
  for(l = N-1; l >= 0; l--) {
    #if(NBINARY==0)
      if(zealot[l]==1)printf("%d ", spin[l]+1);
      else printf("%d ", spin[l]);
    #else
      if(zealot[l]==1)printf("%d ", -spin[l]);
      else printf("%d ", spin[l]);
    #endif
    if( l%L == 0 ) printf("\n");
  }
  printf("e\n");
}

#ifdef SNAPSHOTS
/**************************************************************
 *                       Snapshots                   
 *************************************************************/
  void snap(void) {
    int l;
    int identifier = 0;
    char teste[100];

    lat2eps_init(L,L);
    lat2eps_set_color(0,0x00000); //black
    lat2eps_set_color(1,0xFFFFFF); //white
    lat2eps_set_color(2,0x333434); // black gray
    lat2eps_set_color(3,0xC2C2C2); // white gray
    lat2eps_set_color(6,0xFF0000); // red
    lat2eps_set_color(7,0x0000FF); // blue
    lat2eps_set_color(4,0xff4545); // gray red
    lat2eps_set_color(5,0x81C2EF); // gray blue

    for(l=0; l<N; l++) {
      if(spin[l]==1) {
        if(certainty[l]>=1)lat2eps_set_site(l%L,l/L,4);
        else lat2eps_set_site(l%L,l/L,6);
      }
      else {
        if(certainty[l]>=1)lat2eps_set_site(l%L,l/L,5);
        else lat2eps_set_site(l%L,l/L,5);
      } 
    }

    snprintf(teste,sizeof teste,"sd%ld[%d].eps",seed,identifier);
    while(exists(teste)==true) {
      identifier++;
      snprintf(teste,sizeof teste,"sd%ld[%d].eps",seed,identifier);
    }
    lat2eps_gen_eps(teste,0,0,L,L,1,3);
    lat2eps_release();
  }
#endif


/**************************************************************
 *                    Cluster measures                   
 *************************************************************/
void hoshen_kopelman(void) {
  
  int i,j,temp1,temp2;
  int bigst1=0,bigst2=0;
  mx1=0;
  mx2=0;
  cl1=0;
  numc=0;
  probperc0=0;
  probperc1=0;
  his[N] = 0;
  
  label = malloc(N*sizeof(int));
  siz = malloc(N*sizeof(int));

  for (i=0; i<N; ++i) {
    label[i] = i;
    siz[i] = 0;
    his[i] = 0;
  }

  for (i=0;i < N; ++i) {
    if (spin[i]==spin[right[i]]) connections(i,right[i]);
    if (spin[i]==spin[down[i]]) connections(i,down[i]);
  }

  for (i=0; i<N; ++i) {
    j = i;
    while (label[j] != j)
    {
      j = label[j];
    }
    ++siz[label[j]];                        
    label[i] = label[j];                    
  }

  temp1 = 0;                                 
  temp2 = 0;                                 
  numc = 0;                                  

  for (i=0; i<N; ++i) {
    if (siz[i]>0) {
      ++his[siz[i]];
      if (siz[i]>=temp1) {
        temp2 = temp1;
        temp1 = siz[i];
        bigst2 = bigst1;
        bigst1 = i;
      }
      else if (siz[i]>temp2){
        temp2 = siz[i];
        bigst2 = i;
      } 
      ++(numc); /* Count total number of clusters */
    }
  }

  probperc0 = percolates2d(bigst1);
  if(temp2>0)probperc1 = percolates2d(bigst2);
  
  mx1 = temp1;
  mx2 = temp2;
  
  return;
}

/*****************************************************************************************
 *                              Percolation routines                                     *
 *                           Last modified: 31/07/2020                                   *
 *                                                                                       *
 * Returns: 0 (not spanning/percolating), 1 (spanning) or 2 (percolating)                *
 ****************************************************************************************/
int percolates2d(int site) {

  int ok1,ok2;
  int i,j;
  
  for (i=0; i<L; ++i) {
    ok1=0;
    for (j=i; j<N; j+=L) {
      if (label[j]==site) {
        ok1=1;
        break;         
      }
    }
    if (ok1==0) break;
  }
  for (i=0; i<N; i+=L) {
    ok2=0;
    for (j=i; j<i+L; ++j) {
      if (label[j]==site) {
        ok2=1;
        break;
      }         
    }
    if (ok2==0) break;
  }
  return ok1+ok2;
}

/*******************************************************************************
*                   Instantaneous Percolation Measurements                     *
*                     Last modified:  13/06/2005                               *
*                                                                              *
* Remember that a cluster is identified by the last site in the chain of       *
* pointers.                                                                    *
*******************************************************************************/
void connections(int i,int j) {

  unsigned long i1,j1;

  i1 = label[i];                            /* check where i points to          */
  while (label[i1] != i1)i1 = label[i1];    /* while it doesn't point to itself */

  j1 = label[j];                            /* check where j points to          */
  while (label[j1] != j1)j1 = label[j1];    /* while it doesn't point to itself */

  if (label[i1] > label[j1])label[i1] = label[j1];
  else label[j1] = label[i1];

  return;

}

/*****************************************************************************
 *                            	   Comparison (greater -> smaller)            *
 ****************************************************************************/
int comp (const void *x, const void *y) {
	return - (int) (*(int *)x - *(int*)y);
}

/**************************************************************
 *               Check for duplicate file                  
 *************************************************************/

bool exists(const char *fname){
    if( access( fname, F_OK ) == 0 ) {
        return true;
    } else {
        return false;
    }
}

/**************************************************************
 *               Open output files routine                   
 *************************************************************/

void openfiles(void) {
  char output_file1[100];
  char teste[100];
  unsigned long identifier = seed;

  #if(INTRANS==0)

    #if(NBINARY==0)

      #if(RESET==0)

        snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        }

        #if(DEBUG==1)
          seed=1111111111;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"binarytrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else

        #if(RESET==1)

          snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          }

          #if(DEBUG==1)
          seed=1111111111;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binarytrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          }

          #if(DEBUG==1)
          seed=1111111111;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #endif    

      #endif

    #else

      #if(RESET==0)

        snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        }

        #if(DEBUG==1)
          seed=1111111111;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"nonbinarytrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else

        #if(RESET==1)

          snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          }

          #if(DEBUG==1)
            seed=1111111111;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinarytrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          }

          #if(DEBUG==1)
            seed=1111111111;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #endif

      #endif

    #endif

  #else

    #if(NBINARY==0)

      #if(RESET==0)

        snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        }

        #if(DEBUG==1)
          seed=1111111111;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"binaryintrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else
      
        #if(RESET==1)

          snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          }

          #if(DEBUG==1)
          seed=1111111111;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binaryintrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          }

          #if(DEBUG==1)
          seed=1111111111;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #endif

      #endif

    #else

      #if(RESET==0)

        snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
        }

        #if(DEBUG==1)
          seed=1111111111;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"nonbinaryintrans-ALPHA%.1f-lg%d-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else

        #if(RESET==1)
          snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,identifier);
          }

          #if(DEBUG==1)
            seed=1111111111;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinaryintrans-ALPHA%.1f-lg%d-FULLRESET-DETA-%.5f-seed%ld.dsf",ALPHA,L,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,identifier);
          }

          #if(DEBUG==1)
            seed=1111111111;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-DETA-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,DETA,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;
          
        #endif

      #endif

    #endif

  #endif

  
}
