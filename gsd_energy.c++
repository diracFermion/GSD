#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "gsd.h"
#include "stdint.h"
#include "gsd_tools.h"
#include "gsd_fn.h"
#define KAPPA		5.0
#define EPSILON		800.0
#define a		1
#define M		1

int N,Nb,Nd,i,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
float velocity[NMAX*3],acceleration[NMAX*3];
float position[NMAX*3];
uint32_t particleID[NMAX];
char particleType[5][2];
float u_cross_v[3];
float total_BE,total_SE;
float accln_stretch_x[NMAX],accln_stretch_y[NMAX],accln_stretch_z[NMAX];
float accln_dihedral[3][NMAX];


void print_and_exit(char *format, ...)
{
    va_list list;

    va_start(list,format);
    vprintf(format,list);
    va_end(list);
    exit(1);
}

void load_gsd( char fname[30], uint64_t frame)
{

  struct gsd_handle h;
  
  //Open gsd file with handle
  gsd_open(&h,fname,GSD_OPEN_READONLY);
  //Read number of particles
  gsd_read_chunk(&h,&N,gsd_find_chunk(&h,frame,"particles/N"));
  //Read particle type
  gsd_read_chunk(&h,particleType,gsd_find_chunk(&h,frame,"particles/types"));
  //Read particle TypeId
  gsd_read_chunk(&h,particleID,gsd_find_chunk(&h,frame,"particles/typeid"));
  //Read positions
  gsd_read_chunk(&h,position,gsd_find_chunk(&h,frame,"particles/position"));
  //Read velocity
  gsd_read_chunk(&h,velocity,gsd_find_chunk(&h,frame,"particles/velocity"));
  //Read acceleration
  gsd_read_chunk(&h,acceleration,gsd_find_chunk(&h,frame,"particles/acceleration"));


  if(frame==0)
  {
	  //Read number of bonds
	  gsd_read_chunk(&h,&Nb,gsd_find_chunk(&h,frame,"bonds/N"));
	  //Read the bond groups
	  gsd_read_chunk(&h,bondGroup,gsd_find_chunk(&h,frame,"bonds/group"));
	  //Read number of dihedrals
	  gsd_read_chunk(&h,&Nd,gsd_find_chunk(&h,frame,"dihedrals/N"));
	  //Read dihedral group
	  gsd_read_chunk(&h,dihedralGroup,gsd_find_chunk(&h,frame,"dihedrals/group"));
  }

  //printf("# particles = %d\n",N);

  if(frame==0)
  {
	  printf("\n\n\nSystem Attributes read from Frame 0\n");
	  for(int i=0;i<5;i++)
  	  {
          	printf("%s\n",particleType[i]);
  	  }
	  printf("\n# bonds = %d\n",Nb);
	  for(int i=0;i<Nb;i++)
	  {
		printf("%d %d\n",bondGroup[2*i],bondGroup[2*i+1]);
	  }
	  printf("# dihedrals = %d\n",Nd);
	  for(int i=0;i<Nd;i++)
	  {
		printf("%d %d %d %d\n",dihedralGroup[4*i],dihedralGroup[4*i+1],dihedralGroup[4*i+2],dihedralGroup[4*i+3]);
	  }
	  printf("Particle TypeIDs\n");
  	  for(int i=0;i<N;i++)
  	  {
        	printf("%u\n",particleID[i]);
  	  }
	  printf("************************************************\n");
          printf("************************************************\n\n\n\n");
  }

  //printf("Particle Positions at Frame %d\n",frame);
  for(int i=0;i<N;i++)
  {
        //printf("%f %f %f\n",position[3*i],position[3*i+1],position[3*i+2]);
  }

  //printf("Particle Velocity\n");
  for(int i=0;i<N;i++)
  {
        //printf("%lf %lf %lf\n",velocity[3*i],velocity[3*i+1],velocity[3*i+2]);
  }
  
  //printf("Particle Accelaration\n");
  for(int i=0;i<N;i++)
  {
        //printf("%lf %lf %lf\n",acceleration[3*i],acceleration[3*i+1],acceleration[3*i+2]);
  }
  printf("\n\n");

  return;
}

  

int cross_product_energy(float u[3],float v[3])
{
  u_cross_v[0] = u[1]*v[2] - u[2]*v[1];
  u_cross_v[1] = u[2]*v[0] - u[0]*v[2];
  u_cross_v[2] = u[0]*v[1] - u[1]*v[0];

  float mod_u_cross_v = sqrt(u_cross_v[0]*u_cross_v[0] + u_cross_v[1]*u_cross_v[1] + u_cross_v[2]*u_cross_v[2]);
  
  u_cross_v[0] = u_cross_v[0]/mod_u_cross_v;
  u_cross_v[1] = u_cross_v[1]/mod_u_cross_v;
  u_cross_v[2] = u_cross_v[2]/mod_u_cross_v;

  return 0;
}

int bending_energy()
{
  total_BE=0;
  float vec_cb[3],vec_ab[3],vec_dc[3];
  float A[3],B[3];
  float be,se,dot_AB;
  for(int i=0;i<Nd;i++)
  {
	for(int j=0;j<3;j++)
        {
		vec_cb[j] = position[3*dihedralGroup[4*i+2]+j] - position[3*dihedralGroup[4*i+1]+j];
		vec_ab[j] = position[3*dihedralGroup[4*i]+j] - position[3*dihedralGroup[4*i+1]+j];
		vec_dc[j] = position[3*dihedralGroup[4*i+3]+j] - position[3*dihedralGroup[4*i+2]+j];
	}
	//printf ("Dihedral %d:\t%d %d %d %d\n",i,dihedralGroup[4*i],dihedralGroup[4*i+1],dihedralGroup[4*i+2],dihedralGroup[4*i+3]);
	//printf("vec_cb:\t%lf\t%lf\t%lf\n",vec_cb[0],vec_cb[1],vec_cb[2]);
	//printf("vec_ab:\t%lf\t%lf\t%lf\n",vec_ab[0],vec_ab[1],vec_ab[2]);
	//printf("vec_dc:\t%lf\t%lf\t%lf\n",vec_dc[0],vec_dc[1],vec_dc[2]);
	
	cross_product_energy(vec_cb,vec_ab);
	for(int j=0;j<3;j++)
        {
		A[j] = u_cross_v[j];
	}
	cross_product_energy(vec_cb,vec_dc);
	for(int j=0;j<3;j++)
        {
                B[j] = u_cross_v[j];
        }
	dot_AB = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
	//printf("dot_AB = %lf\n",dot_AB);
	be = 0.5 * KAPPA * (1+dot_AB);
	total_BE = total_BE + be;
	//printf("BE = %lf\n",be);

  }
  return 0;
}

int bond_harmonic_energy()
{
  total_SE = 0;
  float l;//current length of bond
  float se;
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
	total_SE = total_SE + se;
  } 
  return 0;
}

/***********************************************************************************/
/*	Function for calculating accleration from Bond Harmonic Potential	   */
/***********************************************************************************/

int accelaration_bondstretch()
{
  /*	Initializing the accln_stretch_coordinate array	*/
  for(int i=0;i<NMAX;i++)
  {
	accln_stretch_x[i]=0;
	accln_stretch_y[i]=0;
	accln_stretch_z[i]=0;
  }

  float l;//current length of bond
  float ax,ay,az;
  for(int i=0;i<Nb;i++)
  {
        l=0;
        for(int j=0;j<3;j++)
        {
                l = l + (position[3*bondGroup[2*i]+j] - position[3*bondGroup[2*i+1]+j]) * (position[3*bondGroup[2*i]+j] - position[3*bondGroup[2*i+1]+j]);
        }
        l = sqrt(l);
	ax = (-EPSILON/M) * (l-a)/a * (position[3*bondGroup[2*i]] - position[3*bondGroup[2*i+1]]);
	ay = (-EPSILON/M) * (l-a)/a * (position[3*bondGroup[2*i]+1] - position[3*bondGroup[2*i+1]+1]);
	az = (-EPSILON/M) * (l-a)/a * (position[3*bondGroup[2*i]+2] - position[3*bondGroup[2*i+1]+2]);

	accln_stretch_x[bondGroup[2*i]] = accln_stretch_x[bondGroup[2*i]] + ax;
	accln_stretch_y[bondGroup[2*i]] = accln_stretch_y[bondGroup[2*i]] + ay;
	accln_stretch_z[bondGroup[2*i]] = accln_stretch_z[bondGroup[2*i]] + az;

	accln_stretch_x[bondGroup[2*i+1]] = accln_stretch_x[bondGroup[2*i+1]] - ax;
	accln_stretch_y[bondGroup[2*i+1]] = accln_stretch_y[bondGroup[2*i+1]] - ay; 
	accln_stretch_z[bondGroup[2*i+1]] = accln_stretch_z[bondGroup[2*i+1]] - az;
  }
  return 0;
}

/***********************************************************************************/
/*      Functions for calculating accleration from Dihedral Harmonic Potential      */
/***********************************************************************************/
float cross_product(float u[3],float v[3])
{
  u_cross_v[0] = u[1]*v[2] - u[2]*v[1];
  u_cross_v[1] = u[2]*v[0] - u[0]*v[2];
  u_cross_v[2] = u[0]*v[1] - u[1]*v[0];

  float mod_u_cross_v = sqrt(u_cross_v[0]*u_cross_v[0] + u_cross_v[1]*u_cross_v[1] + u_cross_v[2]*u_cross_v[2]);
  
  return mod_u_cross_v;
}

float dot_product(float u[3],float v[3])
{
  return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}

int accelaration_dihedral()
{
  /*	Initializing the accln_dihedral array	*/
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<NMAX;j++)
      {
	accln_dihedral[i][j]=0;
      }
  }

  /*	Generating the Levi-Civita Tensor	*/
  float levi_civita[3][3][3];
  for(int i=0;i<3;i++)
  {
	for(int j=0;j<3;j++)
	{
		for(int k=0;k<3;k++)
		{
			if((i==0 && j==1 && k==2)||(i==1 && j==2 && k==0)||(i==2 && j==0 && k==1))
				levi_civita[i][j][k]=1.0;
			else if((i==0 && j==2 && k==1)||(i==2 && j==1 && k==0)||(i==1 && j==0 && k==2))
				levi_civita[i][j][k] = -1.0;
			else
				levi_civita[i][j][k] = 0;
		}
	}
  } 

  float R_ij[3],R_kj[3],R_lk[3];
  float Rkj_x_Rlk[3],Rkj_x_Rij[3];
  float mod_Rkj_x_Rlk,mod_Rkj_x_Rij;
  for(int i=0;i<Nd;i++)
  {
	/*	Generating the vectors making the dihedral	*/
	for(int b=0;b<3;b++)
	{
		R_kj[b] = position[3*dihedralGroup[4*i+2]+b] - position[3*dihedralGroup[4*i+1]+b];
                R_ij[b] = position[3*dihedralGroup[4*i]+b] - position[3*dihedralGroup[4*i+1]+b];
                R_lk[b] = position[3*dihedralGroup[4*i+3]+b] - position[3*dihedralGroup[4*i+2]+b];
	}
	/*	Generating the two cross products and modulus that go into the dihedral potential calculation	*/
 	mod_Rkj_x_Rlk = cross_product(R_kj,R_lk);
        for(int b=0;b<3;b++)
        {
                Rkj_x_Rlk[b] = u_cross_v[b];
        }
	mod_Rkj_x_Rij = cross_product(R_kj,R_lk);
        for(int b=0;b<3;b++)
        {
                Rkj_x_Rij[b] = u_cross_v[b];
        }
	/*	Accelaration derived from derivative of dihedral harmonic potential	*/
	for(int b=0;b<3;b++)
	{
	  for(int c=0;c<3;c++)
	  {
	    for(int d=0;d<3;d++)
	      {
		/*	Accelaration when node i is wiggled	*/
	        accln_dihedral[b][dihedralGroup[4*i]] += -(KAPPA/(2*M*mod_Rkj_x_Rlk*mod_Rkj_x_Rij)) * (levi_civita[b][c][d] *Rkj_x_Rlk[c]*R_kj[d])  + (KAPPA/(M*mod_Rkj_x_Rlk*pow(mod_Rkj_x_Rij,2))) * (levi_civita[b][c][d]*Rkj_x_Rij[c]*R_kj[d]*dot_product(Rkj_x_Rlk,Rkj_x_Rij));
	
		/*	Accelaration when node j is wiggled   NOT CORRECT one of the 2nd terms goes to 0 CHECK */
		//accln_dihedral[b][dihedralGroup[4*i+1]] += ((KAPPA/(2*M*mod_Rkj_x_Rlk*mod_Rkj_x_Rij)) * (levi_civita[b][c][d] * (Rkj_x_Rij[c] *(R_kj[d]-R_ij[d])))) + (((-KAPPA*dot_product(Rkj_x_Rlk,Rkj_x_Rij))/(2*M*pow(mod_Rkj_x_Rij,3)*mod_Rkj_x_Rlk)*levi_civita[b][c][d] * Rkj_x_Rij[c]*(R_kj[d]-R_ij[d])) + ((KAPPA/(2*M*mod_Rkj_x_Rlk*mod_Rkj_x_Rij)) * (levi_civita[b][c][d] * (Rkj_x_Rlk[c] *(R_kj[d]-R_lk[d])))) + (((-KAPPA*dot_product(Rkj_x_Rlk,Rkj_x_Rij))/(2*M*pow(mod_Rkj_x_Rlk,3)*mod_Rkj_x_Rij)*levi_civita[b][c][d] * Rkj_x_Rlk[c]*(R_kj[d]-R_lk[d]));

		/*      Accelaration when node k is wiggled  NOT CORRECT   */
		//accln_dihedral[b][dihedralGroup[4*i+2]] += ((KAPPA/(2*M*mod_Rkj_x_Rlk*mod_Rkj_x_Rij)) * (levi_civita[b][c][d] * (Rkj_x_Rij[c] *(R_kj[d]-R_ij[d])))) + (((-KAPPA*dot_product(Rkj_x_Rlk,Rkj_x_Rij))/(2*M*pow(mod_Rkj_x_Rij,3)*mod_Rkj_x_Rlk)*levi_civita[b][c][d] * Rkj_x_Rij[c]*(R_kj[d]-R_ij[d])) + ((KAPPA/(2*M*mod_Rkj_x_Rlk*mod_Rkj_x_Rij)) * (levi_civita[b][c][d] * (Rkj_x_Rlk[c] *(R_kj[d]-R_lk[d])))) + (((-KAPPA*dot_product(Rkj_x_Rlk,Rkj_x_Rij))/(2*M*pow(mod_Rkj_x_Rlk,3)*mod_Rkj_x_Rij)*levi_civita[b][c][d] * Rkj_x_Rlk[c]*(R_kj[d]-R_lk[d]));
			
		/*	Accelaration when node l is wiggles	*/
		accln_dihedral[b][dihedralGroup[4*i+3]] += -(KAPPA/(2*M*mod_Rkj_x_Rlk*mod_Rkj_x_Rij)) * (levi_civita[b][c][d] *Rkj_x_Rij[c]*R_lk[d])  + (KAPPA/(M*pow(mod_Rkj_x_Rlk,2)*mod_Rkj_x_Rij,2)) * (levi_civita[b][c][d]*Rkj_x_Rlk[c]*R_kj[d]*dot_product(Rkj_x_Rlk,Rkj_x_Rij));
	      }
	  }
        }
  }
  return 0;
}


int main(int argc, char **argv)
{
  printf("Reading GSD file: %s\n",argv[1]);
  load_gsd(argv[1],0);
  bending_energy();
  bond_harmonic_energy();
  printf("\n\nSystem Bending Energy = %lf\n",total_BE);
  printf("System Bond Harmonic Energy = %lf\n\n",total_SE);
 /* load_gsd("../Sim_dump/trajectory.gsd",atoi(argv[2]));
  bending_energy();
  bond_harmonic_energy();
  printf("\n\nSystem Bending Energy = %lf\n",total_BE);
  printf("System Bond Harmonic Energy = %lf\n",total_SE);
  printf("System Potential Energy = %lf\n\n",total_BE+total_SE);
  *///accelaration_bondstretch();
/*  for(int i=0;i<N;i++)
  {
	printf("%lf %lf %lf\n",accln_stretch_x[i],accln_stretch_y[i],accln_stretch_z[i]);
  }*/
  FILE *fp;
  fp = fopen("../Sim_dump_wall/oscillator.txt", "wa");  
  for(uint64_t i=1000;i<2000;i++)
  {
	load_gsd("../Sim_dump_wall/trajectory.gsd",i);
        fprintf(fp,"%d %lf %lf\n",i,position[0],position[1]);
	printf("%d %lf %lf\n",i,position[0],position[1]);
  }
  fclose(fp);
  return 0;
}
