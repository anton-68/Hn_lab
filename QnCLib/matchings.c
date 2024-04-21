// Enumerating non-isomorph perfect colorings in 4-cube
// By Anton Bondarenko (c)
// Version 0.001 Date 01.05.01

// includings

#include<stdio.h>
#include<stdlib.h>

// data definitions and memory allocation

#define Dimensionality = 4
#define Number_Of_Vertexes 32
#define Size_Of_AutH4 = 384
#define Size_Of_Matching = 8

struct GRAPH
{
  GRAPH *previous = NULL;
  GRAPH *next = NULL;
  int vertex_a = -1;
  int vertex_b = -1;
}

struct GRAPH_LIST
{
  GRAPH *previous = NULL;
  GRAPH *next = NULL;
  GRAPH *current = NULL;
}

struct AUTOMORPHISM
{
  int Dimension_Permutation[Dimensionality];
  int Shift_mask[Dimensionality];
}


main()
{

  GRAPH_LIST *Matchings;

  GRAPH_LIST *Pairs;

  GRAPH_LIST *Colorings;

  GRAPH *Pretender;

  int Matching_mask[Number_Of_Vertexes];

  int n[Size_Of_Matching];

  int i, j, k, l, m, p;

  int h_edges[Number_Of_Vertexes][2] = 
  { { 1, 2 }, { 1, 10 }, { 1, 14 }, { 1, 16 }, { 2, 3 }, { 2, 7 }, { 2, 9 }, 
    { 3, 4 }, { 3, 6 }, { 3, 10 }, { 4, 5 }, { 4, 7 }, { 4, 11 }, { 5, 6 }, 
    { 5, 8 }, { 5, 12 }, { 6, 9 }, { 6, 13 }, { 7, 8 }, { 7, 14 }, { 8, 9 }, 
    { 8, 15 }, { 9, 16 }, { 10, 11 }, { 10, 13 }, { 11, 12 }, { 11, 14 }, 
    { 12, 13 }, { 12, 15 }, { 13, 16 }, { 14, 15 }, { 15, 16 } } ;

  int Vertexes_List[Number_Of_Vertexes*2];

  start(timer_01);

  //
  // finding all non-isomorphic perfect matchings in 4-cube
  //

  if ((Matchings = (GRAPH_LIST *) malloc(sizeof(GRAPH_LIST))) == NULL)
    {
      printf ("\nInsufficient Memory\n");
      return 0;
    }
  if ((Pretender = (GRAPH *) malloc(sizeof(GRAPH))) == NULL)
    {
      printf ("\nInsufficient Memory\n");
      return 0;
    }
  for (i = 0; i < Size_Of_Matching; i++) n[i] = i;
  p = Size_Of_Matching - 1;
  do
    {
      for (j = 0; j < Size_Of_Matching; j++)
	{
	  Vertexes_List[2 * j] = h_edges[n[j]][0];
	  Vertexes_List[(2 * j) + 1] = h_edges[n[j]][1];
	}
      if (IsMatchingVL(Vertexes_list))
	{

	  for (k = 0; k < Size_Of_Matching; k++)
	    Add(h_edges[n[k]], Pretender);
	  if (! ExistIn(Pretender, Matchings))
	    {
	      Add(Pretender, Matchings);
	      if ((Pretender = (GRAPH *) malloc(sizeof(GRAPH))) == NULL)
		{
		  printf ("\nInsufficient Memory\n");
		  return 0;
		}
	    }
	  else
	    Empty(Pretender);
	}
      if (n[Size_Of_Matching - 1] == (Number_Of_Vertexes - 1)) p--;
      else p = Size_Of_Matching - 1;
      if (p >= 0)
	for (i = (Size_Of_Matching - 1); i >= p; i--)
	  n[i] = n[p] + i - p + 1;
    } while (p >= 0);
  free(Pretender);

  //
  // finding complement pairs of matchings
  //

  if ((Pairs = (GRAPH_LIST *) malloc(sizeof(GRAPH_LIST))) == NULL)
    {
      printf ("\nInsufficient Memory\n");
      return 0;
    }
  if ((Pretender = (GRAPH *) malloc(sizeof(GRAPH))) == NULL)
    {
      printf ("\nInsufficient Memory\n");
      return 0;
    }
  for (i = 1; i < Size(Matchings); i++)
    {
      for (j = (i+1); j < Size(Matchings); j++)
	{
	  for (k = 1; k < Size_Of_AutH4; k++)
	    {
	      Pretender = Union(ElementByNumber(Matchings, i), 
				Image(ElementByNumber(Matchings, j), AutH4(k)));
	      if (IsSimpleGraph(Pretender))
		if (!EsistIn(Pretender, Pairs))
		  {
		    Add(Pretender, Pairs);
		    if ((Pretender = (GRAPH *) malloc(sizeof(GRAPH))) == NULL)
		      {
			printf ("\nInsufficient Memory\n");
			return 0;
		      }
		  }
	      else
		Empty(Pretender);
	    }
	}
    }
  free(Pretender);

  //
  // finding perfect colorings
  //

  if ((Colorings = (GRAPH_LIST *) malloc(sizeof(GRAPH_LIST))) == NULL)
    {
      printf ("\nInsufficient Memory\n");
      return 0;
    }
  if ((Pretender = (GRAPH *) malloc(sizeof(GRAPH))) == NULL)
    {
      printf ("\nInsufficient Memory\n");
      return 0;
    }
  for (i = 1; i < Size(Pairs); i++)
    {
      for (j = (i+1); j < Size(Pairs); j++)
	{
	  for (k = 1; k < Size_Of_AutH4; k++)
	    {
	      Pretender = Union(ElementByNumber(Pairs, i), 
				Image(ElementByNumber(Pairs, j), AutH4(k)));
	      if (IsSimpleGraph(Pretender))
		if (!EsistIn(Pretender, Colorings))
		  {
		    Add(Pretender, Colorings);
		      if ((Pretender = (GRAPH *) malloc(sizeof(GRAPH))) == NULL)
			{
			  printf ("\nInsufficient Memory\n");
			  return 0;
			}
		  }
	      else
		Empty(Pretender);
	    }
	}
    }
  free(Pretender);

  stop(timer_01);
  Output(Colorings);
  Print(timer_01);

}





