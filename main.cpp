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

  FILE *fp,*hgt,*wid;
  char filepath[256],init_strip[256],trajectory_file[256],hgt_profile_file[256],hgt_width_file[256];
  double dhe,bhe;
  double backbone_T0,slider_T0;;
  int frame_cnt=0;

  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",NX,NY);
  printf("Init_strip.gsd : %s\n",init_strip);

  // Avg, Height Squared ribbon profile path
  sprintf(hgt_profile_file,"../Sim_dump_ribbon/hgt_prof_L%d_W%d_k%.1f.dat",NX,NY,KAPPA);
  printf("Height Profile File: %s\n",hgt_profile_file);

  hgt = fopen(hgt_profile_file, "w");
  if (hgt == NULL)
   {
	print_and_exit("Could Not Open File to write height profile data");
   }


  for(int run=1;run<=RUN;run++)
  {

	  // Output filepath 
	  sprintf(filepath,"../Sim_dump_ribbon/L%d_W%d_k%.1f_r%d.log",NX,NY,KAPPA,run);
	  printf("Filename of analyzed data: %s\n",filepath);
	  
	  // Trajectory.gsd filepath
	  sprintf(trajectory_file,"../Sim_dump_ribbon/traj_L%d_W%d_k%.1f_r%d.gsd",NX,NY,KAPPA,run);
	  printf("Trajectory File : %s\n",trajectory_file);

	  //Avg Width height of the ribbon
	  sprintf(hgt_width_file,"../Sim_dump_ribbon/width_L%d_W%d_k%.1f_r%d.dat",NX,NY,KAPPA,run);
	  printf("Height width File: %s\n",hgt_width_file);


	  fp = fopen(filepath, "w");
	  if (fp == NULL)
	   {
		print_and_exit("Could Not Open File to write analyzed data");
	   }

	  wid = fopen(hgt_width_file, "w");
  	  if (wid == NULL)
   	  {
        	print_and_exit("Could Not Open File to write height width data");
   	  }

	  //printf("Reading GSD file: %s\n",argv[1]);
	  //load_gsd(argv[1],0);
	  load_gsd(init_strip,0);
	  backbone_T0 = backbone_length(0);
	  //backbone_length(0,fp);
<<<<<<< HEAD
	  fprintf(fp,"Frames\tDihedral_Bending_Energy\tBond_Harmonic_Energy\tPotential_Energy\tDelta_Backbone\tAvg_hgt\tAvg_hgt_Sq\tAvg_Slider_Pos\n");  
	  fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",0,bending_energy(),bond_harmonic_energy(),bending_energy()+bond_harmonic_energy(),backbone_length(0),avg_hgt(),avg_hgt_sq(),avg_slider_pos());
	  
          int c=0;//count of frames > FRAMES/2
=======
	  slider_T0 = avg_slider_pos();
	  fprintf(fp,"Frames\tDihedral_Bending_Energy\tBond_Harmonic_Energy\tPotential_Energy\tDelta_Backbone\tAvg_hgt\tAvg_hgt_Sq\tAvg_Slider_Pos\tDelta_Slider\n");  
	  fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",0,bending_energy(),bond_harmonic_energy(),bending_energy()+bond_harmonic_energy(),backbone_length(0),avg_hgt(),avg_hgt_sq(),avg_slider_pos(),slider_T0-avg_slider_pos());

>>>>>>> 8c70cf4aebd2e11e2657e4b0d590c2555aa844fd
	  for(int frames=1;frames<FRAMES;frames++)
	  {
		//load_gsd(argv[2],frames);
		load_gsd(trajectory_file,frames);
		//backbone_length(frames,fp);
		//printf("%d\t%lf\t%lf\n",frames,position[3*(NX-1)],position[3*(LEN-NX)]);
		dhe = bending_energy();
		bhe = bond_harmonic_energy();
<<<<<<< HEAD
		fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",frames,dhe,bhe,dhe+bhe,backbone_length(frames)-backbone_T0,avg_hgt(),avg_hgt_sq(),avg_slider_pos());
		
=======
		fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",frames,dhe,bhe,dhe+bhe,backbone_length(frames)-backbone_T0,avg_hgt(),avg_hgt_sq(),avg_slider_pos(),slider_T0-avg_slider_pos());
>>>>>>> 8c70cf4aebd2e11e2657e4b0d590c2555aa844fd
		if(frames>FRAMES/2)
		{
			frame_cnt++;
			sum_hgt_node();
			width_hgt(c);
			c++;
			//printf("%d\t",frames - (FRAMES/2 + 1));
		}
		//if(frames == FRAMES/2 + 1)
		//width_hgt(0);			
	  }
	  print_width(wid); 
	  fclose(fp);
	  fclose(wid);
  }

  avg_hgt_node(frame_cnt);

  frame_cnt=0;
  for(int run=1;run<=RUN;run++)
  {
	for(int frames=FRAMES/2;frames<FRAMES;frames++)
	{
		
		// Trajectory.gsd filepath
         	sprintf(trajectory_file,"../Sim_dump_ribbon/traj_L%d_W%d_k%.1f_r%d.gsd",NX,NY,KAPPA,run);
		load_gsd(trajectory_file,frames);
		hgt_profile();
		frame_cnt++;
	}
  }

  avg_hgt_profile(hgt,frame_cnt);
  fclose(hgt);

  return 0;
}
