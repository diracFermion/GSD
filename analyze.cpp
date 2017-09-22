#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

double bendingEner[NMAX];
double bondHarmonicEner[NMAX];
double backbone;
double total_DHE,total_BHE;


double bond_length (int i,int j)
{
  return (sqrt(pow(position[3*i]-position[3*j],2)+pow(position[3*i+1]-position[3*j+1],2)+\
					pow(position[3*i+2]-position[3*j+2],2)));
}

double backbone_length(int frame)
{
  int node;
  backbone=0;
  for(int i=2;i<NX-1;i++)
  {
	node = (NY/2)*NX + i;
	backbone += bond_length(node,node-1);
	//printf("%d\t%d\t%lf\t%.8f\t%.8f\n",node,node-1,bond_length(node,node-1),position[3*node],position[3*(node-1)]);	  
  }
  //fprintf (fp,"%d\t%lf\n",frame,backbone);
  return (backbone);
}


/*	Cross product	*/
double* cross_product(double u[3],double v[3])
{
  static double u_cross_v[3];
  u_cross_v[0] = u[1]*v[2] - u[2]*v[1];
  u_cross_v[1] = u[2]*v[0] - u[0]*v[2];
  u_cross_v[2] = u[0]*v[1] - u[1]*v[0];

  double mod_u_cross_v = sqrt(u_cross_v[0]*u_cross_v[0] + u_cross_v[1]*u_cross_v[1] + u_cross_v[2]*u_cross_v[2]);

  u_cross_v[0] = u_cross_v[0]/mod_u_cross_v;
  u_cross_v[1] = u_cross_v[1]/mod_u_cross_v;
  u_cross_v[2] = u_cross_v[2]/mod_u_cross_v;

  //printf("CP:\t%lf,%lf,%lf\n",u_cross_v[0],u_cross_v[1],u_cross_v[2]);
  return u_cross_v ;
}


/*	Function evaluating Dihedral Harmonic Energy	*/
double bending_energy()
{
  double vec_cb[3],vec_ab[3],vec_dc[3];
  double *A,*B,V_A[3],V_B[3];
  double be,dot_AB;
  double total_DHE = 0;

  for(int i=0;i<Nd;i++)
  {
        for(int j=0;j<3;j++)
        {
                vec_cb[j] = position[3*dihedralGroup[4*i+2]+j] - position[3*dihedralGroup[4*i+1]+j];
                vec_ab[j] = position[3*dihedralGroup[4*i]+j] - position[3*dihedralGroup[4*i+1]+j];
                vec_dc[j] = position[3*dihedralGroup[4*i+3]+j] - position[3*dihedralGroup[4*i+2]+j];
        }
        //printf ("Dihedral %d:\t%d %d %d %d\n",i,dihedralGroup[4*i],dihedralGroup[4*i+1],dihedralGroup[4*i+2],dihedralGroup[4*i+3]);
        //printf("vec_cb:\t%lf,%lf,%lf\n",vec_cb[0],vec_cb[1],vec_cb[2]);
        //printf("vec_ab:\t%lf,%lf,%lf\n",vec_ab[0],vec_ab[1],vec_ab[2]);
        //printf("vec_dc:\t%lf,%lf,%lf\n",vec_dc[0],vec_dc[1],vec_dc[2]);

        A = cross_product(vec_cb,vec_ab);
	for(int k=0;k<3;k++)
	{
		V_A[k]=*(A+k);		
	}
        B = cross_product(vec_cb,vec_dc);
	for(int k=0;k<3;k++)
        {
                V_B[k]=*(B+k);
        }
	//printf("vec_cb X vec_ab : %.8f,%.8f,%.8f\n",V_A[0],V_A[1],V_A[2]);
	//printf("vec_cb X vec_dc : %.8f,%.8f,%.8f\n",V_B[0],V_B[1],V_B[2]);

        dot_AB = V_A[0]*V_B[0]+V_A[1]*V_B[1]+V_A[2]*V_B[2];
        //printf("dot_AB = %lf\n",dot_AB);
        be = 0.5 * KAPPA * (1+dot_AB);//Using HOOMD kappa

	total_DHE += be;
        //printf("BE = %lf\n",be);

  }
  return(total_DHE);
}

/*	Function evaluating Bond Harmonic Energy	*/
double bond_harmonic_energy()
{
  total_BHE = 0;
  double l;//current length of bond
  double se;
  for(int i=0;i<Nb;i++)
  {
	l=0;
  	for(int j=0;j<3;j++)
  	{
		l = l + (position[3*bondGroup[2*i]+j] - position[3*bondGroup[2*i+1]+j]) * (position[3*bondGroup[2*i]+j] - position[3*bondGroup[2*i+1]+j]);
	}
	l = sqrt(l);
	se = 0.5 * EPSILON * (l-a) * (l-a);
	//printf("Bond %d %d , se = %lf\n",bondGroup[2*i],bondGroup[2*i+1],se);
	total_BHE = total_BHE + se;
  } 
  return (total_BHE);
}

/*	Function evaluating Average Z height above z=0 plane	*/
double avg_hgt()
{
   double hgt=0;
   for(int i=0;i<N;i++)
   {
	hgt+=position[3*i+2];
   }
   return (hgt/N);
}

/*	Average <h^2> = 1/N * Sum_i(z_i - <z>)^2	*/
double avg_hgt_sq()
{
   double hgtSq=0;
   double h_avg = avg_hgt();
   for(int i=0;i<N;i++)
   {
        hgtSq+=pow((position[3*i+2]-h_avg),2);
   }
   return (hgtSq/N);
   
}





