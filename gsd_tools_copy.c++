#include "stdio.h"
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "gsd.h"
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

double orientation2angle(float q1, float q2)
{ // function to convert HOOMD's 2d (x,y) quaternion values to radians
  static double theta;
  if ( q1 >= 0 ) {
    if ( q2 >= 0 )
      theta = atan(q2/q1);
    else
      theta = atan(q2/q1) + 2*M_PI;
  } else
    theta = atan(q2/q1) + M_PI;
  
  theta -= M_PI/2.;
  return fmod( theta + 2*M_PI, 2*M_PI );
}


void load_gsd( char fname[30], uint64_t frame)
{

  struct gsd_handle h;
  int N,i;
  float position[NMAX*3];
  //float orientation[NMAX*4];
  //int image[NMAX*3];
  //double angle[NMAX];

  //Open gsd file with handle
  gsd_open(&h,fname,GSD_OPEN_READONLY);
  //Read number of particles
  gsd_read_chunk(&h,&N,gsd_find_chunk(&h,frame,"particles/N"));
  //Read positions
  gsd_read_chunk(&h,position,gsd_find_chunk(&h,frame,"particles/position"));
  //Read orientations
  //gsd_read_chunk(&h,orientation,gsd_find_chunk(&h,frame,"particles/orientation"));
  //Read images
  //gsd_read_chunk(&h,image,gsd_find_chunk(&h,frame,"particles/image"));

/*
  for(i=0; i<N;i++)
  {
    q[i].x = double(position[3*i]  ) + 0.5*sys.Lx();
    q[i].y = double(position[3*i+1]) + 0.5*sys.Ly();
    q[i].x_abs = q[i].x + double(image[3*i]  )*sys.Lx();
    q[i].y_abs = q[i].y + double(image[3*i+1])*sys.Ly();
    q[i].theta = orientation2angle(orientation[4*i+1],orientation[4*i+2]);
  }
*/
//  for(i=0; i <N;i++)
//    printf("%lf %lf %lf %d %d\n", q[i].x,q[i].y,q[i].theta,image[3*i],image[3*i+1]);
  printf("# particles = %d\n",N);
  return;
}

int main()
{
  load_gsd(init_strip.gsd"",0);
  return 0;
}
