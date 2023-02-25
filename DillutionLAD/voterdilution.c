/*************************************************************************
*                     Votante Latoski-Arenzon-Dantas 2D                  *
*                             V1.07 18/06/2021                           *
*************************************************************************/

/***************************************************************
 *                    OBRIGATORY DEFINITIONS
 **************************************************************/
// -D L="LATTICE_SIZE"

/***************************************************************
 *                     OPTIONAL DEFINITIONS
 **************************************************************/

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
#define MCS         1E6 //Max evolution time
#define MEASURES    40
#define ALPHA       1.  //Transiten probability 1
#define BETA        1.  //Transiten probability 2
#define GAMMA       1.05  //Gamma reset scale

/****************************************************************
 *                            SETTINGS 
 ***************************************************************/

#define LOGSCALE    1 // 0 --> measures logaritmically spaced, 1 --> measures in logscale.

/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/

void initialize(void); 
void openfiles(void); 
void sweep(void); 
void visualize(double,unsigned long); 
void states(void); 
void measures1(void); 
void measures2(void);
void listremove(int);
void listinclude(int);
void single_update(int);
void updatelist(int);
#ifdef SNAPSHOTS
  void snap(void);  
#endif
void hoshen_kopelman(void);
int biasedwalk(int qual, int *lab);
int delta(int i, int j, int hh);
void connections(int,int);
int percolates2d(int);
bool exists(const char*);
bool probcheck(double);

/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

FILE *fp1,*fp2;
int *spin,*list,*listaux,**neigh,*memory,*measures,*right,*left,*up, *down, sum, sumz, activesum;
int *siz, *label, *his, *qt, cl1, numc, mx1, mx2,CONT,LINKS;
int NACTIVE,probperc0,probperc1;
int hull_perimeter;
char root_name[200];
unsigned long seed;
double tempo;

/***************************************************************
 *                          MAIN PROGRAM  
 **************************************************************/
int main(void){

  #if(DEBUG==0)
    #if(SEED==0)
      seed = time(0);
      if (seed%2==0) ++seed;
    #else
      seed = SEED;
    #endif
  #else
      seed = 1111111111;
  #endif

  openfiles(); 



  int k=0;
  initialize();
  int contagem=0;
  while(tempo <= MCS){
    if( NACTIVE == 0 ){
        states();
        hoshen_kopelman();
        fprintf(fp1,"%.8f %.8f %.8f %.8f %.8f %d %.8f %d %.8f %d\n",tempo,(double)sum/CONT,(double)activesum/LINKS,(double)numc/CONT,(double)mx1/CONT,probperc0,(double)mx2/CONT,probperc1,(double)qt[0]/CONT,NACTIVE);
        fflush(fp1);
        while(measures[k]!=0){
          fprintf(fp1,"%.8f %.8f %.8f %.8f %.8f %d %.8f %d %.8f %d\n",(double)measures[k],(double)sum/CONT,(double)activesum/LINKS,(double)numc/CONT,(double)mx1/CONT,probperc0,(double)mx2/CONT,probperc1,(double)qt[0]/CONT,NACTIVE);
          fflush(fp1);
          k++;
        }           
        break;
    }
    if (tempo >= measures[k]) {  
      states();
      hoshen_kopelman();
      fprintf(fp1,"%.8f %.8f %.8f %.8f %.8f %d %.8f %d %.8f %d\n",tempo,(double)sum/CONT,(double)activesum/LINKS,(double)numc/CONT,(double)mx1/CONT,probperc0,(double)mx2/CONT,probperc1,(double)qt[0]/CONT,NACTIVE);
      fflush(fp1);
      k++;
    }
    int agent = FRANDOM*NACTIVE;
    int site = list[agent];
    //printf("Atualizando agente %d (de %d ativos) localizado no sitio %d com spin %d\n",agent,NACTIVE,site,spin[site]);
    #if(VISUAL==1)
      if(tempo>contagem){
        visualize(tempo,seed);
        contagem++;
      }
    #endif
    single_update(site);
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
  right = malloc(N*sizeof(int));
  left = malloc(N*sizeof(int));
  up = malloc(N*sizeof(int));
  down = malloc(N*sizeof(int));
  list = malloc(N*sizeof(int));
  listaux = malloc(N*sizeof(int));
  #if(NBINARY==0)
    qt = malloc(2*sizeof(int));
  #else
    qt = malloc(N*sizeof(int));
  #endif

  for(int i=0; i<N; i++){
    neigh[i] = (int*)malloc(4*sizeof(int));
    #if(NBINARY==1)
      qt[i] = 0;
    #endif
  }

  CONT=0;
  for(int n=0; n<N; n++) { 
    double RAND = FRANDOM;
    memory[n] = 0;
    if(RAND<=RHO){      
      int k=FRANDOM*2;
      spin[n] = k*2 - 1; 
      qt[k]++;
      CONT++;
    }
    else {
      spin[n] = 0;
    }
  }
  #if(VISUAL==0)
    fprintf(fp1,"# Agents: %d\n",CONT);
    fprintf(fp1,"# Time Persistence Active Clusters Big1 Perc1 Big2 Perc2 qt[0]\n");
    fprintf(fp1,"\n\n");
    fflush(fp1);
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

  LINKS=0;
  NACTIVE=0;
  for (int i=0; i<N; i++) {
    listaux[i]=-1;
    int relative=0;
    if(spin[i]!=0){
      for(int j=0;j<4;j++)if(spin[neigh[i][j]] == -spin[i])relative=1;
      if (spin[right[i]]!=0){
        LINKS++;
        if(spin[right[i]]!=spin[i]){
          activesum++;
        }
      }
      if (spin[down[i]]!=0){
        LINKS++;
        if(spin[down[i]]!=spin[i]){
          activesum++;
        }
      }
      if(relative==1){
        list[NACTIVE]=i;
        listaux[i]=NACTIVE;
        NACTIVE++;
      }
    }
  }
  #if(LOGSCALE==1)
    measures1(); 
  #else
    measures2();
  #endif

}
/****************************************************************
 *               Single_Update
 ***************************************************************/
void single_update(int _site) {
  //Opinion dynamics    
  int dir = FRANDOM*4;
  int neighbour = neigh[_site][dir];
  if(spin[_site]==-spin[neighbour]) {
    int INTERFANTES=0;
    if(spin[up[_site]]==-spin[_site])INTERFANTES++;
    if(spin[right[_site]]==-spin[_site])INTERFANTES++;
    if(spin[down[_site]]==-spin[_site])INTERFANTES++;
    if(spin[left[_site]]==-spin[_site])INTERFANTES++;
    memory[_site]=1;
    qt[(spin[_site] + 1 )/2]--;
    spin[_site] = spin[neighbour];
    qt[(spin[neighbour] + 1 )/2]++;
    int INTERFDEPOIS=0;
    if(spin[up[_site]]==-spin[_site])INTERFDEPOIS++;
    if(spin[right[_site]]==-spin[_site])INTERFDEPOIS++;
    if(spin[down[_site]]==-spin[_site])INTERFDEPOIS++;
    if(spin[left[_site]]==-spin[_site])INTERFDEPOIS++;
    activesum+=(INTERFDEPOIS-INTERFANTES);
    tempo+=1./NACTIVE;
    updatelist(_site);
  }
}

/****************************************************************
 *               Check states numbers
 ***************************************************************/
void states(void) {
  sum=CONT;
  LINKS=0;
  // activesum=0;
  for (int i=0; i<N; i++) {
    if(memory[i]!=0)sum--;
    if(spin[i]!=0){
      if (spin[right[i]]!=0)LINKS++;          
      if (spin[down[i]]!=0)LINKS++;
    }
    //printf("%d\n",sum);
  }
}

/**************************************************************
 *              List management                  
 *************************************************************/
void listremove(int _site){
  int agent = listaux[_site];
  list[agent] = list[NACTIVE];
  listaux[list[agent]]=agent;
  listaux[_site]=-1;
  NACTIVE--;
}

void listinclude(int _site){
  NACTIVE++;
  list[NACTIVE] = _site;
  listaux[_site] = NACTIVE;
}

void updatelist(int _site){
  int cont=0;
  int cont2;
  for(int i=0; i<4; i++){
    int viz = neigh[_site][i];
    if(spin[viz] == -spin[_site])cont=1;
    cont2=0;
    if(spin[viz]!=0){
      for(int j=0; j<4; j++){
        int vizviz = neigh[viz][j];
        if(spin[vizviz] == -spin[viz])cont2=1;
      }
      if(listaux[viz]==-1 && cont2==1)listinclude(viz);
      if(listaux[viz]!=-1 && cont2==0)listremove(viz); 
    }
  }
  if(cont == 0)listremove(_site); 
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
 *                       Vizualização                   
 *************************************************************/
void visualize(double _j,unsigned long _seed) {
  int l;
  printf("pl '-' matrix w image t 'time = %.8f seed = %ld'\n",_j,_seed);
  for(l = N-1; l >= 0; l--) {
    if(spin[l]!=0){  
      printf("%d ", spin[l]);
    }
    else{
      printf("%d ", -2);
    }
    if( l%L == 0 ) printf("\n");
  }
  printf("e\n");
}


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
    if(spin[i]!=0){
      if (spin[i]==spin[right[i]]) connections(i,right[i]);
      if (spin[i]==spin[down[i]]) connections(i,down[i]);
    }

  }

  for (i=0; i<N; ++i) {
    if(spin[i]!=0){
      j = i;
      while (label[j] != j)
      {
        j = label[j];
      }
      ++siz[label[j]];                        
      label[i] = label[j];                    
    }
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
  char output_file1[300];
  char teste[250];
  
  sprintf(root_name,"datavoterdillution_lg%d_rho%.2f",L,RHO);

  unsigned long identifier = seed;
  #if(DEBUG==0)
    sprintf(teste,"%s_sd%ld_1.dsf",root_name,identifier);
    while(exists(teste)==true) {
      identifier+=2;
      sprintf(teste,"%s_sd%ld_1.dsf",root_name,identifier);
    }
  #endif
  sprintf(teste,"%s_sd%ld",root_name,identifier);
  
  seed=identifier;

  sprintf(output_file1,"%s_1.dsf",teste);
  fp1 = fopen(output_file1,"w");
  fprintf(fp1,"# Pure Voter Model 2D Main Output\n");
  fprintf(fp1,"# Seed: %ld\n",seed);
  fprintf(fp1,"# Linear size: %d\n",L);
  fflush(fp1);

  return;
  
}

