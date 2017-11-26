#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include "variables.h"
#include "jack_knife.h"


int DATA_COUNT,BIN_SIZE,RUNS;
int NX,NY,STEPS,LOGGING,LOG2BIN;
double EPSILON,KAPPA,k_bT;

int main( int argc, char **argv )
{

  //File containing the runs for which Error Estimation will be done
  char validrunfile[256];
  FILE *fp;

  /*    Logging toggle to print to terminal     */
  switch (argc){
   case 7:
       sscanf(argv[1],"%d",&NX);
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%lf",&KAPPA);
       sscanf(argv[4],"%d",&STEPS);
       sscanf(argv[5],"%lf",&k_bT);
       sscanf(argv[6],"%d",&LOGGING);
       break;
   default:
       print_and_exit("Usage Pass command line arguments: NX NY KAPPA STEPS kBT 0 or 1 for LOGGING\n");
   }

  sprintf(validrunfile,"../Sim_dump_ribbon/L%d/W%d/k%.1f/valid_slider_runs.log",NX,NY,KAPPA);

  if(NULL==(fp=fopen(validrunfile,"r")))
        print_and_exit("I could not open file with simulation run numbers %s\n",validrunfile);

  /*    Epsilon */
  EPSILON = 720.0 * KAPPA;

  /*    Total MD Steps and Jack Knife Bin Count */
  DATA_COUNT = STEPS/PERIOD;

  /*    Printing read paramater.dat data        */
  if (LOGGING == 1)
  {
          printf("NX = %d\n",NX);
          printf("NY = %d\n",NY);
          printf("STEPS\t%d \n",STEPS);
          printf("PERIOD\t%d \n",PERIOD);
          printf("EPSILON\t%f \n",EPSILON);
          printf("KAPPA\t%f \n",KAPPA);
          printf("DATA_COUNT\t%d \n",DATA_COUNT);
          printf("runfile\t%s\n",validrunfile);
  }

 log2_single_observable_time_evolution(DATA_COUNT,fp);

 fclose(fp);
 return 0;

}

