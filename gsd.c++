#include "stdio.h"
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "gsd.h"
#include "stdint.h"
//#include ".../System.h"
//#include "../Simulation.h"
#include "gsd_tools.h"
#include "gsd_fn.h"

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
  int N,Nb,Nd,i,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
  float position[NMAX*3];
  uint32_t particleID[NMAX];
  int8_t particleType[5];

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
  //Read number of bonds
  gsd_read_chunk(&h,&Nb,gsd_find_chunk(&h,frame,"bonds/N"));
  //Read the bond groups
  gsd_read_chunk(&h,bondGroup,gsd_find_chunk(&h,frame,"bonds/group"));
  //Read number of dihedrals
  gsd_read_chunk(&h,&Nd,gsd_find_chunk(&h,frame,"dihedrals/N"));
  //Read dihedral group
  gsd_read_chunk(&h,dihedralGroup,gsd_find_chunk(&h,frame,"dihedrals/group"));

  printf("# particles = %d\n",N);
//  printf("# particleType = %u\n",particleType);
  for(int i=0;i<N;i++)
  {
	printf("%lf %lf %lf\n",position[3*i],position[3*i+1],position[3*i+2]);
  }
  for(int i=0;i<5;i++)
  {
        printf("%u\n",particleType[i]);
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
  return;
}

int main(int argc, char **argv)
{
  printf("Reading GSD file: %s\n",argv[1]);
  load_gsd(argv[1],0);
  return 0;
}
