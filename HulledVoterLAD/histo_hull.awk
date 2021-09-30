#!/bin/bash
# Sintax: awk -f histo_2d_int -v NUMPAR=... -v bin=... 
# LCFL 29/09/2021
{

if (($1=="#")&&($2=="Domain:"))
   {
      DOMAIN=$3;
   }

if (($1=="#")&&($2=="Linear"))
   {
      if(DOMAIN=="Rough")
      {
         L=$4;
        nome=sprintf("HULLrough_%d_t%d_b%d",L,TIME,bin);
      }
      else
      {
         L=$4;
         nome=sprintf("HULLsmooth_%d_t%d_b%d",L,TIME,bin);
      }
   }

if (($1=="#")&&($2="Sample")) samples=$4; 

if (($1!="#") && (NF>0))
   {
    coluna[int($1/bin)] += $3;
    sets += $3;
    if ($1>imax) imax = $1;
   }
}
END {      
     if (bin>1) delta = 0.5*bin - 1;
     for (j=imax; j>=1; --j)
         {
          if (coluna[j]>0)
            printf("%d %d %d %d\n",j*bin+delta,coluna[j],sets,samples) > nome;
	      }

   }
