#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
//#include "gsd_tools.h"
//#include "gsd_fn.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"


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
