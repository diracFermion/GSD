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

  FILE *fp;
  char filepath[256];
  double dhe,bhe;

  sprintf(filepath,argv[3]);
  printf("Filename of analyzed data: %s\n",filepath);  
  fp = fopen(filepath, "w");
  if (fp == NULL)
   {
        print_and_exit("Could Not Open File to write analyzed data");
   }

  printf("Reading GSD file: %s\n",argv[1]);
  load_gsd(argv[1],0);
  //backbone_length(0,fp);
  fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\n",0,bending_energy(),bond_harmonic_energy(),bending_energy()+bond_harmonic_energy(),backbone_length(0));

  for(int frames=1;frames<FRAMES;frames++)
  {
	load_gsd(argv[2],frames);
	//backbone_length(frames,fp);
	//printf("%d\t%lf\t%lf\n",frames,position[3*(NX-1)+2],position[3*(LEN-NX)+2]);
	dhe = bending_energy();
	bhe = bond_harmonic_energy();
	fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\n",frames,dhe,bhe,dhe+bhe,backbone_length(frames));
  }
  fclose(fp);
  return 0;
}
