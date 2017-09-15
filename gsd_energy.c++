#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "gsd.h"
#include "stdint.h"
#include "gsd_tools.h"
#include "gsd_fn.h"
#define KAPPA 5
#define FRAMES 201
#define NX 50
#define LEN 2900

int N,Nb,Nd,i,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
float position[NMAX*3];
uint32_t particleID[NMAX];
char particleType[3][2];
float bendingEner[NMAX];
float bondHarmonicEner[NMAX];



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
//  printf("# particleType = %u\n",particleType);
  for(int i=0;i<N;i++)
  {
	//printf("%lf %lf %lf\n",position[3*i],position[3*i+1],position[3*i+2]);
  }
  for(int i=0;i<3;i++)
  {
        //printf("%s\n",particleType[i]);
  }

  if(frame==0)
  {
	  printf("\n# bonds = %d\n",Nb);
	  for(int i=0;i<Nb;i++)
	  {
		//printf("%d %d\n",bondGroup[2*i],bondGroup[2*i+1]);
	  }
	  printf("# dihedrals = %d\n",Nd);
	  for(int i=0;i<Nd;i++)
	  {
		//printf("%d %d %d %d\n",dihedralGroup[4*i],dihedralGroup[4*i+1],dihedralGroup[4*i+2],dihedralGroup[4*i+3]);
	  }
  }

  //printf("Particle TypeIDs\n");
  for(int i=0;i<N;i++)
  {
        //printf("%u\n",particleID[i]);
  } 
  return;
}

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


int main(int argc, char **argv)
{
  printf("Reading GSD file: %s\n",argv[1]);
  load_gsd(argv[1],0);
  //bending_energy();
  //bond_harmonic_energy();
  for(int frames=1;frames<FRAMES;frames++)
  {
	load_gsd(argv[2],frames);
	printf("%d\t%lf\t%lf\n",frames,position[3*(NX-1)+2],position[3*(LEN-NX)+2]);
  }
  return 0;
}
