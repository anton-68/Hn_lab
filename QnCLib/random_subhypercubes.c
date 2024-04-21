/* Hamiltonian Cycles Problem. Subhypercubes test sets.    	       		     */
/* Data format: EDGE_DATA_FORMAT : EDGE_LIST (TSPLIB)					           */
/* Anton Bondarenko	<a_bond@rinet.ru>	                       April 1, 1999  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int   Connected( unsigned, unsigned );

FILE     *output;

main( int argc, char **args)
{
   static char RF[] = "hypercubes.hsc";
   char *ResultsFile = RF;
   unsigned v1, v2, dim, nov;

/* Defining input&output                                                     */

   args++;
   if ( *args != 0 )
   {
      dim = atoi(*args);
      args++;
      if ( *args != 0 ) ResultsFile = *args;
   }
   else
   	printf("\Usage: hypercubes <dimensionality> [<output filename>]");

/* Initalising output                                                        */

   if( (output = fopen( ResultsFile , "w" )) == NULL )
   {
      printf("\nError02: Unable to open output file");
      exit(1);
   }

   nov = 1;
   for ( v1 = 1; v1 <= dim; v1++ ) nov = nov*2;

   fprintf( output, "NAME : %s", ResultsFile );
   fprintf( output, "\nCOMMENT : Hamiltonian cycle problem (Hypercubes)" );
   fprintf( output, "\nTYPE : HCP" );
   fprintf( output, "\nDIMENSION : %8d", nov );
   fprintf( output, "\nEDGE_DATA_FORMAT : EDGE_LIST" );
   fprintf( output, "\nEDGE_DATA_SECTION" );

/* Here we go...																				  */

   for( v1 = nov; v1 > 0; v1-- )
   	for( v2 = (v1 - 1); v2 > 0; v2-- )
      	if( Connected (v1, v2) )
         	fprintf( output, "\n%8d%8d", v1, v2 );

   fprintf( output, "\n-1" );

   return 0;
}

/*---------- Connected ------------------------------------------------------*/

int Connected( unsigned IV1, unsigned IV2 )
{
	unsigned V1, V2;
   V1=--IV1;
   V2=--IV2;
	if (fabs(log(V1 ^ V2)/log(2) - ceil(log(V1 ^ V2)/log(2))) <= 0.0000000001 |
       fabs(log(V1 ^ V2)/log(2) - floor(log(V1 ^ V2)/log(2))) <= 0.0000000001 )
   	return 1;
   else
   	return 0;
}

