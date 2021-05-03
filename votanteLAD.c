/*************************************************************************
*                     Votante Latoski-Arenzon-Dantas                     *
*                             V1.02 03/05/2021                           *
*************************************************************************/

// -DNBINARY [non binary case]
// -DRESET  [full reset case]
// -DGRESET [gamma reset case]
// -DINTRANS [intrans case]
// -DINCA [inc 1e-1], -DINCB [inc 1e-2], -DINCC [inc 1e-3], -DINCD [inc 1e-4], -DINCE [inc 1e-5]
// -DCLUSTER [swap to cluster measures]

// -DSPEEDTEST [speed test]
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

#define L           128   //Lado da Rede
#define N          (L*L)  //Número de Sítios
#define MCS         1000000 //MCS para as medidas
#define THRESHOLD   1. //Conviccao maxima
#define MEASURES    40
#define SAMPLES     1
#define ALPHA       1.  //Probabilidade de Transiçã (pode depender da convicção eventualmente)
#define BETA        1.
#define SEED        1616583043   //Debug purposes
#ifdef INCA
  #define INCREMENTO  0.1
#endif
#ifdef INCB
  #define INCREMENTO  0.01
#endif
#ifdef INCC
  #define INCREMENTO  0.001
#endif
#ifdef INCD
  #define INCREMENTO  0.0001
#endif
#ifdef INCE
  #define INCREMENTO  0.00001
#endif
#define GAMMA       2.

/****************************************************************
 *                            SETTINGS 
 ***************************************************************/

#ifdef INTRANS
  #define INTRANS     1 // 0 --> reversible, 1 --> irreversible
#endif
#ifdef NBINARY
  #define BINARY      0 // 0 --> m(0)=2 ,        
#else
  #define BINARY      1
#endif
#ifdef RESET
  #define RESET       1 // 0 --> without reset, 1 --> full reset, 2 --> gamma reset
#endif
#ifdef GRESET
  #define RESET       2
#endif
#ifdef DEBUG
  #define DEBUG       1 // 1 --> fixed seed
#endif
#ifdef VISUAL
  #define VISUAL      1 // 1 --> run with ./a.out | gnuplot for gifs
#endif
#ifdef SPEEDTEST
  #define SPEED_TEST  1 // 1 --> mcs speed test
#endif
#ifdef SNAPSHOTS
  #define SNAPSHOTS   1 // 1 --> ignore other measures and simply take snapshots of the system       
#endif
#define LOGSCALE    1 // 0 --> measures logaritmically spaced, 1 --> measures in logscale
#define SIMPLIFIED  1

/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/

void inicializacao(void); //Inicia a rede
void openfiles(void); //Abre arquivos de saida (ainda não utilizo)
void sweep(void); //Evolui
void vizualizacao(void); //Habilita vizualizacao
void somas(void);
void states(void); //Vê o numero de estados difereentes no sistema
void measures1(void); //Cria vetor de medidas
void measures2(void);
#ifdef SNAPSHOTS
  void snap(void);  
#endif
#if(CLUSTER==1)
  void hoshen_kopelman(void);
  void connections(int,int);
  int percolates2d(int);
#endif
bool exists(const char*);
bool teste(double);

/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

FILE *fp1,*fp2;
int *spin,**viz,soma,*memory,*measures,*zealot,*right,*left,*up, *down, soma, somaz, conexoes;
#if(CLUSTER==1)
  int *siz, *label, *his0, *his0_perc, cl1, numc, mx1, mx2, aga;
  int probperc0,probperc1;
#endif
unsigned long seed;
double *certainty;

/***************************************************************
 *                          MAIN PROGRAM  
 **************************************************************/
int main(void){

  for (int S=0; S<SAMPLES; S++) {   
    int i,j,k=0;
    seed = time(0);
    
    #if(SNAPSHOTS==0)
      openfiles();
    #endif
    
    inicializacao();
    
    for (j=0;j<=MCS+1;j++)  {
      if (measures[k]==j) {  

        #if(SNAPSHOTS==1)
          snap();   
          k++;        
        #else  
          #if(SPEED_TEST==1)
            clock_t t=clock();
          #endif
          #if(CLUSTER==0)
            states();
            if(somaz==N){
              while(measures[k]!=0){
                fprintf(fp1,"%d %.6f %.6f %.6f\n",measures[k],(double)soma/N,(double)somaz/N,(double)conexoes/N);
                k++;
              }           
              break;
            }
            fprintf(fp1,"%d %.6f %.6f %.6f\n",j,(double)soma/N,(double)somaz/N,(double)conexoes/N);
            k++;
          #else 
            hoshen_kopelman();
            if(numc==1){
              while(measures[k]!=0){
                fprintf(fp1,"%d %d %d %d %d\n", measures[k],numc,probperc0,mx1,mx2);
                k++;
              }
              break;
            }
            fprintf(fp1,"%d %d %d %d %d\n", j,numc,probperc0,mx1,mx2);
            k++;
          #endif
          #if(SPEED_TEST==1)
            t=clock()-t;
            double time_taken = ((double)t)/CLOCKS_PER_SEC;
            printf("one measure: %fs\n",time_taken);
          #endif
        #endif

      }

      #if(SPEED_TEST==1)
        clock_t t=clock();
        sweep();
        t=clock()-t;
        double time_taken = ((double)t)/CLOCKS_PER_SEC;
        printf("onesweep: %fs\n",time_taken);
      #else
        sweep();
      #endif

      #if(VISUAL==1 )
        vizualizacao();
      #endif

    }

    #if(SNAPSHOTS==0)
    fclose(fp1);
    #endif
  
  }

}
/***************************************************************
 *                        INICIALIZAÇÃO  
 **************************************************************/
void inicializacao(void) {
 
  start_randomic(seed);

  #ifdef CLUSTER
    his0 = malloc(N*sizeof(int));
    his0_perc = malloc(N*sizeof(int));
  #endif

  spin = malloc(N*sizeof(int));
  viz = (int**)malloc(N*sizeof(int*));
  memory = malloc(N*sizeof(int));
  measures = malloc(MCS*sizeof(int));
  zealot = malloc(N*sizeof(int));
  right = malloc(N*sizeof(int));
  left = malloc(N*sizeof(int));
  up = malloc(N*sizeof(int));
  down = malloc(N*sizeof(int));
  certainty = malloc(N*sizeof(double));

  for(int i=0; i<N; i++){
    viz[i] = (int*)malloc(4*sizeof(int));
  }

  for(int n=0; n<N; n++) { 
    certainty[n] = 0;
    zealot[n] = 0;
    memory[n] = 0;

    #if(BINARY==1)
      int k=FRANDOM*2;
      spin[n] = k*2 - 1; 
    #else
      spin[n] = n;
    #endif

  } 

  for (int i = 0; i < N; i++) {    
    viz[i][0] = (i+1)%L + (i/L)*L; //dir
    viz[i][1] = (i-1+L)%L + (i/L)*L; //esq
    viz[i][2] = (i-L+N)%N; //cima
    viz[i][3] = (i+L)%N; //baixo
    right[i] = viz[i][0];
    left[i] = viz[i][1];
    up[i] = viz[i][2];
    down[i] = viz[i][3];
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
    int sitio = FRANDOM*N;
    int dir = FRANDOM*4;
    int vizinho = viz[sitio][dir];
    #if(SIMPLIFIED==1)
      if(spin[sitio]!=spin[vizinho]) {
        if(zealot[sitio] == 0){

          #if(BINARY==1)
            memory[sitio]=1;
            spin[sitio] = spin[vizinho];
          #else
            memory[spin[sitio]]--;
            spin[sitio] = spin[vizinho];
            memory[spin[sitio]]++;
          #endif
        }
        certainty[vizinho] += INCREMENTO;

        #if(RESET==2)
          certainty[sitio] = certainty[sitio]/GAMMA;
        #endif
        #if(RESET==1)
          certainty[sitio] = 0;
        #endif
        #if(RESET==0)
          certainty[sitio] -= INCREMENTO;
        #endif     

        #if(INTRANS==0)
          if(certainty[sitio]<=THRESHOLD)zealot[sitio]=0;
        #endif

        if(certainty[vizinho]>=THRESHOLD)zealot[vizinho]=1;
        continue;
      }
      else{
        certainty[sitio] += INCREMENTO;
        certainty[vizinho] += INCREMENTO;
        if(certainty[sitio]>=THRESHOLD)zealot[sitio]=1;
        if(certainty[vizinho]>=THRESHOLD)zealot[vizinho]=1;
        continue;
      }
      
    #else
      bool acc1 = teste(ALPHA);
      bool acc2 = teste(BETA);
      if(spin[sitio]!=spin[vizinho]) {
        if(zealot[sitio] == 0){
          if(acc1==true) {

            #if(BINARY==1)
              memory[sitio]=1;
              spin[sitio] = spin[vizinho];
            #else
              memory[spin[sitio]]--;
              spin[sitio] = spin[vizinho];
              memory[spin[sitio]]++;
            #endif

            #if(RESET==2)
              certainty[sitio] = certainty[sitio]/GAMMA;
            #endif

            #if(RESET==1)
              certainty[sitio] = 0;
            #endif

            #if(RESET==0)
              certainty[sitio] -= INCREMENTO;
            #endif

            certainty[vizinho] += INCREMENTO;
            #if(INTRANS==0)
              if(certainty[sitio]<=THRESHOLD)zealot[sitio]=0;
            #endif
            if(certainty[vizinho]>=THRESHOLD)zealot[vizinho]=1;

          }
          else{
            certainty[sitio] += INCREMENTO;
            certainty[vizinho] -= INCREMENTO;

            #if(INTRANS==0)
              if(certainty[vizinho]<=THRESHOLD)zealot[vizinho]=0;
            #endif

            if(certainty[sitio]>=THRESHOLD)zealot[sitio]=1;
          }
        }

        else {
          if (acc2==true) {

            certainty[vizinho] += INCREMENTO;

            #if(RESET==2)
              certainty[sitio] = certainty[sitio]/GAMMA;
            #endif

            #if(RESET==1)
              certainty[sitio] = 0;
            #endif

            #if(RESET==0)
              certainty[sitio] -= INCREMENTO;
            #endif
            
            #if(INTRANS==0)
              if(certainty[sitio]<=THRESHOLD)zealot[sitio]=0;
            #endif
            if(certainty[vizinho]>=THRESHOLD)zealot[vizinho]=1;

          }

          else  {
            certainty[sitio] -= INCREMENTO;
            certainty[vizinho] += INCREMENTO;

            #if(INTRANS==0)
              if(certainty[sitio]<=THRESHOLD)zealot[sitio]=0;
            #endif

            if(certainty[vizinho]>=THRESHOLD)zealot[vizinho]=1;

          }

        }

      }

      else{
        certainty[sitio] += INCREMENTO;
        certainty[vizinho] += INCREMENTO;
        if(certainty[sitio]>=THRESHOLD)zealot[sitio]=1;
        if(certainty[vizinho]>=THRESHOLD)zealot[vizinho]=1;
        continue;
      }
    #endif    
  }
}

/****************************************************************
 *               Check states numbers
 ***************************************************************/
void states(void) {
  #if(BINARY==1)  
    soma=N;
  #else
    soma=0;
  #endif
  somaz=0;
  conexoes=0;
  for (int i=0; i<N; i++) {
    #if(BINARY==1)
      if (memory[i]!=0) soma--;
    #else
      if (memory[i]>=0) soma++;
    #endif
    
    if (zealot[i]!=0) somaz++;
    if (spin[right[i]]!=spin[i]) conexoes++;
    if (spin[down[i]]!=spin[i]) conexoes++;
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
 *                       Vizualização                   
 *************************************************************/
#ifdef SNAPSHOTS
  void snap(void) {
    int l;
    int identifier = 0;
    char teste[100];

    lat2eps_init(L,L);
    lat2eps_set_color(0,0x00000); //black
    lat2eps_set_color(1,0xFFFFFF); //white
    lat2eps_set_color(2,0x333434); // black gray
    lat2eps_set_color(3,0xC2C2C2); // white gray
    lat2eps_set_color(4,0xFF0000); // red
    lat2eps_set_color(5,0x31a2f2); // blue
    lat2eps_set_color(6,0xff4545); // gray red
    lat2eps_set_color(7,0x81C2EF); // gray blue

    for(l=0; l<N; l++) {
      if(spin[l]==1) {
        if(certainty[l]>=1)lat2eps_set_site(l%L,l/L,6);
        else lat2eps_set_site(l%L,l/L,4);
      }
      else {
        if(certainty[l]>=1)lat2eps_set_site(l%L,l/L,7);
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
 *                       Vizualização                   
 *************************************************************/
void vizualizacao(void) {
  int l;
  printf("pl '-' matrix w image\n");
  for(l = N-1; l >= 0; l--) {
    #if(BINARY==1)
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

/**************************************************************
 *                      Teste
 *************************************************************/

 bool teste(double _ALPHA) {
    double r = FRANDOM;
    if(r<=_ALPHA)return true;
    else return false;
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

#if(CLUSTER==1)

/**************************************************************
 *                    Cluster measures                   
 *************************************************************/
void hoshen_kopelman(void) {
  
  int i,j,temp1,temp2,p0;

  mx1=0;
  mx2=0;
  cl1=0;
  numc=0;
  probperc0=0;
  probperc1=0;
  his0[N] = 0;
  
  label = malloc(N*sizeof(int));
  siz = malloc(N*sizeof(int));

  for (i=0; i<N; ++i) {
    label[i] = i;
    siz[i] = 0;
    his0[i] = 0;
  }

  for (i=0;i < N; ++i) {
    if (spin[i]==spin[right[i]]) connections(i,right[i]);
    if (spin[i]==spin[down[i]]) connections(i,down[i]);
  }

  for (i=0; i<N; ++i) {
    j = i;
    int k=0;
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
      ++his0[siz[i]];
      if (siz[i]>=temp1) {
        temp2 = temp1;
        temp1 = siz[i];
      }
      else if (siz[i]>temp2) temp2 = siz[i];       
      p0 = percolates2d(i);
      if(p0==0)++his0_perc[siz[i]];
      else probperc0 = p0;
      ++(numc); /* Count total number of clusters */
    }
  }

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
while (label[i1] != i1)                   /* while it doesn't point to itself */
      i1 = label[i1];

j1 = label[j];                            /* check where j points to          */
while (label[j1] != j1)                   /* while it doesn't point to itself */
      j1 = label[j1];

if (label[i1] > label[j1]) label[i1] = label[j1];
                  else label[j1] = label[i1];
return;
}

#endif


/**************************************************************
 *               Open output files routine                   
 *************************************************************/
void openfiles(void) {
  char output_file1[100];
  char teste[100];
  unsigned long identifier = seed;

  #if(INTRANS==0)

    #if(BINARY==1)

      #if(RESET==0)

        snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        }

        #if(DEBUG==1)
          seed=SEED;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"binarytrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else

        #if(RESET==1)

          snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
          seed=SEED;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binarytrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
          seed=SEED;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #endif    

      #endif

    #else

      #if(RESET==0)

        snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        }

        #if(DEBUG==1)
          seed=SEED;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"nonbinarytrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else

        #if(RESET==1)

          snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
            seed=SEED;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinarytrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
            seed=SEED;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinarytrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #endif

      #endif

    #endif

  #else

    #if(BINARY==1)

      #if(RESET==0)

        snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        }

        #if(DEBUG==1)
          seed=SEED;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"binaryintrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else
      
        #if(RESET==1)

          snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
          seed=SEED;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binaryintrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"binaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
          seed=SEED;
          #else
          seed=identifier;
          #endif

          sprintf(output_file1,"binaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #endif

      #endif

    #else

      #if(RESET==0)

        snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        while(exists(teste)==true) {
          identifier++;
          snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
        }

        #if(DEBUG==1)
          seed=SEED;
        #else
          seed=identifier;
        #endif

        sprintf(output_file1,"nonbinaryintrans-ALPHA%.1f-lg%d-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
        fp1 = fopen(output_file1,"w");
        fflush(fp1);
        return;

      #else

        #if(RESET==1)
          snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
            seed=SEED;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinaryintrans-ALPHA%.1f-lg%d-FULLRESET-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;

        #else

          snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          while(exists(teste)==true) {
            identifier++;
            snprintf(teste,sizeof teste,"nonbinaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,identifier);
          }

          #if(DEBUG==1)
            seed=SEED;
          #else
            seed=identifier;
          #endif

          sprintf(output_file1,"nonbinaryintrans-ALPHA%.1f-lg%d-GAMMA%.1f-INCREMENTO-%.5f-seed%ld.dsf",ALPHA,L,GAMMA,INCREMENTO,seed);
          fp1 = fopen(output_file1,"w");
          fflush(fp1);
          return;
          
        #endif

      #endif

    #endif

  #endif

  
}
