#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

float bendingEner[NMAX];
float bondHarmonicEner[NMAX];


float* cross_product(float u[3],float v[3])
{
  static float u_cross_v[3];
  u_cross_v[0] = u[1]*v[2] - u[2]*v[1];
  u_cross_v[1] = u[2]*v[0] - u[0]*v[2];
  u_cross_v[2] = u[0]*v[1] - u[1]*v[0];

  float mod_u_cross_v = sqrt(u_cross_v[0]*u_cross_v[0] + u_cross_v[1]*u_cross_v[1] + u_cross_v[2]*u_cross_v[2]);

  u_cross_v[0] = u_cross_v[0]/mod_u_cross_v;
  u_cross_v[1] = u_cross_v[1]/mod_u_cross_v;
  u_cross_v[2] = u_cross_v[2]/mod_u_cross_v;

  return u_cross_v ;
}


int bending_energy()
{
  float vec_cb[3],vec_ab[3],vec_dc[3];
  float *A,*B;
  float be,dot_AB;
  for(int i=0;i<Nd;i++)
  {
        for(int j=0;j<3;j++)
        {
                vec_cb[j] = position[3*dihedralGroup[4*i+2]+j] - position[3*dihedralGroup[4*i+1]+j];
                vec_ab[j] = position[3*dihedralGroup[4*i]+j] - position[3*dihedralGroup[4*i+1]+j];
                vec_dc[j] = position[3*dihedralGroup[4*i+3]+j] - position[3*dihedralGroup[4*i+2]+j];
        }
        printf ("Dihedral %d:\t%d %d %d %d\n",i,dihedralGroup[4*i],dihedralGroup[4*i+1],dihedralGroup[4*i+2],dihedralGroup[4*i+3]);
        printf("vec_cb:\t%lf\t%lf\t%lf\n",vec_cb[0],vec_cb[1],vec_cb[2]);
        printf("vec_ab:\t%lf\t%lf\t%lf\n",vec_ab[0],vec_ab[1],vec_ab[2]);
        printf("vec_dc:\t%lf\t%lf\t%lf\n",vec_dc[0],vec_dc[1],vec_dc[2]);

        A = cross_product(vec_cb,vec_ab);
        B = cross_product(vec_cb,vec_dc);
        dot_AB = (*A) * (*B) + (*(A+1)) * (*(B+1)) + (*(A+2)) * (*(B+2));
        printf("dot_AB = %lf\n",dot_AB);
        be = 0.5 * KAPPA * (1+dot_AB);

        printf("BE = %lf\n",be);

  }
  return 0;
}

