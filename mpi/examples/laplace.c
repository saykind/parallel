#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NR 50   
#define NC 100

int main (int argc, char *argv[]) 
{
int	i, j, k;
double	old[NR][NC],
	new[NR][NC],
	delta;
delta = 100.0;

/*** Initialization and boundary conditions ***/

  for (i=0; i<NR; i++)
    for (j=0; j<NC; j++)
      old[i][j]= 0.0;

  for (j=0; j<NC; j++)
      old[0][j]=9.49*sin(3.14159265*j/(NC-1.0));

/*
printf("******************************************************\n");
printf("Initial Matrix:\n");
for (i=0; i<NR; i++)
  {
  for (j=0; j<NC; j++) 
    printf("%1.0f", old[i][j]);
  printf("\n"); 
  }
printf("******************************************************\n");
*/

  while (delta>0.000001) 
  {

    for (i=1; i<NR-1; i++)    
     for (j=1; j<NC-1; j++)
      new[i][j]= 0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]);
  
    delta=0.0;
    for (i=1; i<NR-1; i++)    
     for (j=1; j<NC-1; j++)
     {
	delta+=(new[i][j]-old[i][j])*(new[i][j]-old[i][j]);
        old[i][j]=new[i][j];
     }
    delta=sqrt(delta);
 }

/*** Print results ***/
printf("******************************************************\n");
printf("Delta = %10.7f\n",delta);
printf("Result Matrix:\n");
for (i=0; i<NR; i++)
  {
  for (j=0; j<NC; j++) 
    printf("%1.0f", old[i][j]);
  printf("\n"); 
  }
printf("******************************************************\n");
printf ("Done.\n");

}
