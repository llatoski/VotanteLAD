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

#define LOGSCALE    0 // 0 --> measures logaritmically spaced, 1 --> measures in logscale.
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
int *spin,**neigh,*memory,*measures,*zealot,*right,*left,*up, *down, sum, sumz, activesum;
int *siz, *label, *his, *qt, cl1, numc, mx1, mx2;
int probperc0,probperc1;
int hull_perimeter;
char root_name[200];
unsigned long seed;
double *certainty;

/***************************************************************
 *                          MAIN PROGRAM  
 **************************************************************/
int main(void){

  #if(SEED==0)
    seed = time(0);
    if (seed%2==0) ++seed;
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

  for(int n=0; n<N; n++) { 
    certainty[n] = 0;
    zealot[n] = 0;
    memory[n] = 0;

    #if(NBINARY==0)
      int k=FRANDOM*2;
      spin[n] = k*2 - 1; 
      qt[k]++;
    #else
      spin[n] = n;
    #endif

  } 

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
    int dir = FRANDOM*4;
    int neighbour = neigh[site][dir];

    #if(SIMPLIFIED==1)
      if(spin[site]!=spin[neighbour]) {
        if(zealot[site] == 0){

          #if(NBINARY==0)
            memory[site]=1;
            qt[(spin[site] + 1 )/2]--;
            spin[site] = spin[neighbour];
            qt[(spin[neighbour] + 1 )/2]++;
          #else
            memory[spin[site]]--;
            spin[site] = spin[neighbour];
            memory[spin[site]]++;
          #endif

        }

        certainty[neighbour] += DETA;

        #if(RESET==2)
          certainty[site] = certainty[site]/GAMMA;
        #endif
        #if(RESET==1)
          certainty[site] = 0;
        #endif
        #if(RESET==0)
          certainty[site] -= DETA;
        #endif     

        #if(INTRANS==0)
          if(certainty[site]<=THRESHOLD)zealot[site]=0;
        #endif

        if(certainty[neighbour]>=THRESHOLD)zealot[neighbour]=1;
        continue;
      }
      else{
        certainty[site] += DETA;
        certainty[neighbour] += DETA;
        if(certainty[site]>=THRESHOLD)zealot[site]=1;
        if(certainty[neighbour]>=THRESHOLD)zealot[neighbour]=1;
        continue;
      }
    #else
      bool acc1 = probcheck(ALPHA);
      bool acc2 = probcheck(BETA);
      if(spin[site]!=spin[neighbour]) {
        if(zealot[site] == 0){
          if(acc1==true) {

            #if(NBINARY==0)
              memory[site]=1;
              spin[site] = spin[neighbour];
            #else
              memory[spin[site]]--;
              spin[site] = spin[neighbour];
              memory[spin[site]]++;
            #endif

            #if(RESET==2)
              certainty[site] = certainty[site]/GAMMA;
            #endif

            #if(RESET==1)
              certainty[site] = 0;
            #endif

            #if(RESET==0)
              certainty[site] -= DETA;
            #endif

            certainty[neighbour] += DETA;
            #if(INTRANS==0)
              if(certainty[site]<=THRESHOLD)zealot[site]=0;
            #endif
            if(certainty[neighbour]>=THRESHOLD)zealot[neighbour]=1;

          }
          else{
            certainty[site] += DETA;
            certainty[neighbour] -= DETA;

            #if(INTRANS==0)
              if(certainty[neighbour]<=THRESHOLD)zealot[neighbour]=0;
            #endif

            if(certainty[site]>=THRESHOLD)zealot[site]=1;
          }
        }

        else {
          if (acc2==true) {

            certainty[neighbour] += DETA;

            #if(RESET==2)
              certainty[site] = certainty[site]/GAMMA;
            #endif

            #if(RESET==1)
              certainty[site] = 0;
            #endif

            #if(RESET==0)
              certainty[site] -= DETA;
            #endif
            
            #if(INTRANS==0)
              if(certainty[site]<=THRESHOLD)zealot[site]=0;
            #endif
            if(certainty[neighbour]>=THRESHOLD)zealot[neighbour]=1;

          }

          else  {
            certainty[site] -= DETA;
            certainty[neighbour] += DETA;

            #if(INTRANS==0)
              if(certainty[site]<=THRESHOLD)zealot[site]=0;
            #endif

            if(certainty[neighbour]>=THRESHOLD)zealot[neighbour]=1;

          }

        }

      }

      else{
        certainty[site] += DETA;
        certainty[neighbour] += DETA;
        if(certainty[site]>=THRESHOLD)zealot[site]=1;
        if(certainty[neighbour]>=THRESHOLD)zealot[neighbour]=1;
        continue;
      }
    #endif    

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
    lat2eps_set_color(4,0xFF0000); // red
    lat2eps_set_color(5,0x0000FF); // blue
    lat2eps_set_color(6,0xff4545); // gray red
    lat2eps_set_color(7,0x81C2EF); // gray blue

    for(l=0; l<N; l++) {
      if(spin[l]==1) {
        if(certainty[l]>=1)lat2eps_set_site(l%L,l/L,4);
        else lat2eps_set_site(l%L,l/L,6);
      }
      else {
        if(certainty[l]>=1)lat2eps_set_site(l%L,l/L,5);
        else lat2eps_set_site(l%L,l/L,7);
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
    if (spin[i]==spin[right[i]] && zealot[i]==zealot[right[i]]) connections(i,right[i]);
    if (spin[i]==spin[down[i]] && zealot[i]==zealot[down[i]]) connections(i,down[i]);
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


/*************************************************************************
*                     Biased walks along the external hull               *
*                          Last modified: 31/07/2020                     *
* Returns the area enclosed by the hull, despite the presence of smaller *
* internal domains. We depart from the cluster labelling site (smaller   *
* index), and the initial contribution to the area is L+1. Then we   *
* attempt to walk along the left, front, right or backward directions.   *
* The incoming (backward) direction is only accepted if it's not         *
* possible to go on the other directions. We choose the height of the    *
* initial size as 10*L and update the height along the walk.         *
* Percolating clusters are not taken into account since the walls are    *
* disconnected in that case. The area is also updated along the walk     *
* using  the following table (the column shows the precedent step):      *
*                                                                        *
*              up     right    down   left                               *
*  up (0)       0      h+1      h+1     0                                *
*  right (3)    0      h+1      -h      1                                *
*  down (2)    -h      0        0      -h                                *
*  left (1)    -h      1        0      h+1                               *
*                                                                        *
*                                                                        *
*************************************************************************/
int biasedwalk(int qual, int *lab)

{
int i=0,y;
int area;
int dir=0,old,ok,endpoint=1;
int num_visited=0;
long unsigned int *visited;
visited = jmalloc(N*sizeof(int));

/* At first, we should find the starting point for our walk around the hull.
Some domains may cross the horizontal border, and the label site will not be
the top-leftmost site, so we have to find it, starting from the bottom layer up: */
if (qual<L) 
   {
    ok=0;
    y = up[0];
    while (!ok)
          {
           for (i=0; i<L; ++i)
               if (lab[i+y]==lab[qual]) break;
           if (i==L) ok=1;
                    else {
                          qual = y+i;
                          y = up[y];
                         }
          }
   }
/* but perhaps the domain also crosses the vertical border, so we try to find an initial
site more to the left (not necessarily the leftmost): */
while (lab[qual]==lab[left[qual]]) qual = left[qual];

y = 10*L; /* the factor 10 is arbitrary, but should be at least 2 */
area = y + 1;

/* First check whether the starting point has one or two branches going out. If there
   are two branches connected by the starting point, the configuration shoulbe be:  11
                                                                                    10
   The walk will first go through the horizontal branch and, in order to go to the down
   branch it should pass over the initial site, but we should not stop the walk there. 
   To do this, we set endpoint=0. When the walk returns from the horizontal branch, we 
   set it to 1, so the walk can stop the next time it visits the starting point. Notice
   however that the horizontal branch may close the loop and join the down one without
   passing through the starting point. For example:  111
                                                     101
                                                     111                              */
if ((lab[right[qual]]==lab[qual]) && (lab[down[qual]]==lab[qual]) &&
    (lab[down[right[qual]]]!=lab[qual]))
   endpoint = 0;

/* from the starting point, choose the direction to move, there are just 2 possibilities
(from the way we choose it): */
if (lab[right[qual]]==lab[qual]) 
   {
    dir = 3; 
    i = right[qual];
    visited[0] = up[qual];
   }
   else if (lab[down[qual]]==lab[qual]) 
           {
            dir = 2; 
            i = down[qual]; 
            --y;
	    visited[0] = right[qual];
           }
num_visited = 1;

/* start the walk around the cluster, clockwise: */
old = dir;

while ((i!=qual) || ((!endpoint)&&(dir!=0))) 
/* while the walk doesn't return to the starting point, or if it returns, it can continue to the other branch */
   {
    if (i==qual) endpoint = 1; /* next time it returns, the walk is over */
    ok = 0;
    dir = (dir+1)%4;  /* from the incoming direction, try left first */
    while (!ok) /* from the incoming direction: try right, in front, left and backwards */
      {
       switch (dir)
	 {
	  case 0: if (lab[up[i]]==lab[i]) {
	                                   ok = 1; 
					   i = up[i]; 
					   area += delta(old,dir,y); 
					   ++y; 
	                                  }
	          else {
                        visited[num_visited]=up[i];
		        ++num_visited;
                       }
		  break;
	  case 1: if (lab[left[i]]==lab[i]) {
	                                     ok = 1; 
					     i = left[i];
					     area += delta(old,dir,y); 
	                                    }
                  else {
		        visited[num_visited]=left[i];
		        ++num_visited;
                       }
		  break;
	  case 2: if (lab[down[i]]==lab[i]) {
                                             ok = 1; 
					     i = down[i];
					     area += delta(old,dir,y); 
					     --y;
	                                    }
                  else {
		        visited[num_visited]=down[i];
		        ++num_visited;
                       }
		  break;
	  case 3: if (lab[right[i]]==lab[i]) {
                                              ok = 1;
					      i = right[i];
					      area += delta(old,dir,y); 
	                                     }
                  else {
   	                visited[num_visited]=right[i];
		        ++num_visited;
                       }
		  break;
	 }
       if (ok==0) dir = (dir + 3)%4;
      }
    old = dir;
   }

if (dir==1) area -= y;

if (lab[right[i]]!=lab[i]) { /* add the surface sites around the first/last site */
                            visited[num_visited]=right[i];
		            ++num_visited;
                           }
if (lab[left[i]]!=lab[i]) {
                           visited[num_visited]=left[i];
		           ++num_visited;
                          }
if (lab[down[i]]!=lab[i]) {
                           visited[num_visited]=down[i];
		           ++num_visited;
                          }
if (lab[up[i]]!=lab[i]) {
                         visited[num_visited]=up[i];
	                 ++num_visited;
                        }
hull_perimeter = uniq(visited,num_visited);
free(visited);

return area;
}


int delta(int i, int j, int hh)
{
switch (i)
   {
    case 0: switch (j)
              {
	      case 0: return 0;
	      case 1: return 0;
	      case 2: return (hh + 1);
	      case 3: return (hh + 1);
	      }
   case 1: switch (j)
              {
	      case 0: return (-hh);
	      case 1: return (-hh);
	      case 2: return 0;
	      case 3: return 1;
	      }
   case 2: switch (j)
              {
	      case 0: return (-hh);
	      case 1: return (-hh);
	      case 2: return 0;
	      case 3: return 0;
	      }
   case 3: switch (j)
              {
	      case 0: return 0;
	      case 1: return 1;
	      case 2: return (hh + 1);
	      case 3: return (hh + 1);
	      }
   default: printf("Wrong situation!\n"); exit(0);
   }
 return 1;
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
  char teste[300];
  unsigned long identifier = seed;

  #if((INTRANS==0)&&(NBINARY==0)&&(RESET==0))
    sprintf(root_name,"binarytrans-ALPHA%.1f-L%d-DETA%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==0)&&(NBINARY==0)&&(RESET==1))
    sprintf(root_name,"binarytrans-ALPHA%.1f-L%d-FULLRESET-DETA-%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==0)&&(NBINARY==0)&&(RESET==2))
    sprintf(root_name,"binarytrans-ALPHA%.1f-L%d-GAMMA%.1f-DETA-%.5f",ALPHA,L,GAMMA,DETA);
  #endif
  #if((INTRANS==0)&&(NBINARY==1)&&(RESET==0))
    sprintf(root_name,"nonbinarytrans-ALPHA%.1f-L%d-DETA%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==0)&&(NBINARY==1)&&(RESET==1))
    sprintf(root_name,"nonbinarytrans-ALPHA%.1f-L%d-FULLRESET-DETA-%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==0)&&(NBINARY==1)&&(RESET==2))
    sprintf(root_name,"nonbinarytrans-ALPHA%.1f-L%d-GAMMA%.1f-DETA-%.5f",ALPHA,L,GAMMA,DETA);
  #endif
  #if((INTRANS==1)&&(NBINARY==0)&&(RESET==0))
    sprintf(root_name,"binaryintrans-ALPHA%.1f-L%d-DETA%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==1)&&(NBINARY==0)&&(RESET==1))
    sprintf(root_name,"binaryintrans-ALPHA%.1f-L%d-FULLRESET-DETA-%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==1)&&(NBINARY==0)&&(RESET==2))
    sprintf(root_name,"binaryintrans-ALPHA%.1f-L%d-GAMMA%.1f-DETA-%.5f",ALPHA,L,GAMMA,DETA);
  #endif
  #if((INTRANS==1)&&(NBINARY==1)&&(RESET==0))
    sprintf(root_name,"nonbinaryintrans-ALPHA%.1f-L%d-DETA%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==1)&&(NBINARY==1)&&(RESET==1))
    sprintf(root_name,"nonbinaryintrans-ALPHA%.1f-L%d-FULLRESET-DETA-%.5f",ALPHA,L,DETA);
  #endif
  #if((INTRANS==1)&&(NBINARY==1)&&(RESET==2))
    sprintf(root_name,"nonbinaryintrans-ALPHA%.1f-L%d-GAMMA%.1f-DETA-%.5f",ALPHA,L,GAMMA,DETA);
  #endif

  sprintf(teste,"%s_sd%ld.dsf",root_name,identifier);
  while(exists(teste)==true) {
    identifier+=2;
    sprintf(teste,"%s_sd%ld.dsf",root_name,identifier);
  }

  #if(DEBUG==1)
  seed=1111111111;
  #else
  seed=identifier;
  #endif

  sprintf(output_file1,"%s",teste);
  fp1 = fopen(output_file1,"w");
  fflush(fp1);
  return;

}