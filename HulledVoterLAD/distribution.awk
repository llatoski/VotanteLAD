#!/bin/bash
# LCFL 29/09/2021
BEGIN { 
       printL=0;
       files = 0;
       samples = 0;
      }
{

if (($1=="#")&&($2=="LAD"))
   {
      domain=$6;
   }

if (($1=="#")&&($2=="Linear"))
   {
    if (printL==0)
       {
        print "# Domain:",domain;
        print $0; 
        print "# Time:",TIME;
        L = $4; 
        printL=1;
       }
    primeiro = 1;
    ++files;
   }

   if (($1=="#")&&($2=="Time:"))
   {
      tempo=$3;
   }


if ((NF>0)&&($1 !="#")) 
   {
    if (tempo==TIME)
       {
	if (primeiro == 1)
           {
    	    ++sample;
            print "# Sample =",sample;
	    primeiro = 0;
	   }
        print $0;
       }
   }

if ((NF==0)&&(primeiro == 1)) 
   {
    tempo += 1;
    if (tempo==TIME)
       {
	++sample;
        print "# Sample =",sample;
	primeiro = 0;
       }
   }
}
END {
    print "# Files processed:",files," Samples:",sample;
    }
