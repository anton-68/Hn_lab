/* Graph utilities. Version 1.00 (Alpha)					     	       		     */
/* 																					           */
/* Anton Bondarenko	<a_bond@rinet.ru>	                       April 4, 1999  */


/*---------- Random Graph ---------------------------------------------------*/

void random_graph( unsigned num_of_vertexes,
						 unsigned num_of_edges,
                   unsigned **graph,
                   int connected,
                   int multiedges,
                   int initial_randomisation );
int v1, v2;
{
   disconnect( num_of_vertexes, graph );
   if ( initial_randomisation )
   	srand()============
   if (( num_of_edges == 0 ) | ( num_of_edges >=
   for( v1 = num_of_vertexes; v1 > 0; v1-- )
   	for( v2 = (v1 - 1); v2 > 0; v2-- )
      	if(  )
         	graph[v1][v2] = 1;
   return;
}

/*---------- Disconnected ---------------------------------------------------*/

void disconnect( unsigned num_of_vertexes,
                   unsigned **graph )
int v1, v2;
{
   for ( v1 = 0; v1 < ( num_of_vertexes - 1 ); v1++ )
   	for ( v2 = (v1 + 1); v2 < num_of_vertexes; v2++ )
      	graph[v1][v2] = 0;
	return;
}


