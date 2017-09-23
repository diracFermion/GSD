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
  char filepath[256],init_strip[256],trajectory_file[256];
  double dhe,bhe;
  double backbone_T0;

  //sprintf(filepath,argv[3]);
  // Output filepath
  sprintf(filepath,"../Sim_dump_ribbon/TE_L%dW_%d_k%.1f_r%d.log",NX,NY,KAPPA,RUN);
  printf("Filename of analyzed data: %s\n",filepath);
  
  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",NX,NY);  

  // Trajectory.gsd filepath
  sprintf(trajectory_file,"../Sim_dump_ribbon/traj_L%dW_%d_k%.1f_r%d.gsd",NX,NY,KAPPA,RUN);


  fp = fopen(filepath, "w");
  if (fp == NULL)
   {
        print_and_exit("Could Not Open File to write analyzed data");
   }

  //printf("Reading GSD file: %s\n",argv[1]);
  
  //load_gsd(argv[1],0);
  load_gsd(init_strip,0);
  backbone_T0 = backbone_length(0);
  //backbone_length(0,fp);
  fprintf(fp,"Frames\tDihedral_Bending_Energy\tBond_Harmonic_Energy\tPotential_Energy\tDelta_Backbone\tAvg_hgt\tAvg_hgt_Sq\tAvg_Slider_Pos\n");  
  fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",0,bending_energy(),bond_harmonic_energy(),bending_energy()+bond_harmonic_energy(),backbone_length(0),avg_hgt(),avg_hgt_sq(),avg_slider_pos());

  for(int frames=1;frames<FRAMES;frames++)
  {
	//load_gsd(argv[2],frames);
        load_gsd(trajectory_file,frames);
	//backbone_length(frames,fp);
	//printf("%d\t%lf\t%lf\n",frames,position[3*(NX-1)+2],position[3*(LEN-NX)+2]);
	dhe = bending_energy();
	bhe = bond_harmonic_energy();
	fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",frames,dhe,bhe,dhe+bhe,backbone_length(frames)-backbone_T0,avg_hgt(),avg_hgt_sq(),avg_slider_pos());
  }
  fclose(fp);
  return 0;
}
