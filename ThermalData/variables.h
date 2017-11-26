#define MAXLOG 200
#define MAXRUNS 501
#define MAXMES  100000
#define NUMOBSER 8
#define PERIOD 10000

/*      Observable variables    */
extern int DATA_COUNT,BIN_SIZE,JK_BIN_COUNT,RUNS;
extern char parameter_file[1024],result_file[1024];
extern FILE *Finput,*Fresult;
extern int NX,NY,STEPS,LOGGING,LOG2BIN;
extern double EPSILON,KAPPA,k_bT;

typedef struct
{
 int step;
 double temp;
 double potential_energy;
 double bond_harmonic_energy;
 double kinetic_energy;
 double dihedral_harmonic_energy;
} Observable;

extern int raw_data_marker,log_bin_size;
extern unsigned long long raw_step[MAXMES];
extern char observable_file[1024],time_evolution_dir[1024],time_evolution_file[1024];
extern double raw_data[NUMOBSER][MAXRUNS][MAXMES],jk_blocks[MAXLOG][MAXRUNS];
extern double sum_log2_bin;
extern double log2_bin[MAXRUNS][MAXLOG],step_log2[MAXLOG];
extern double jk_avg[NUMOBSER][MAXLOG];
extern double error[NUMOBSER][MAXLOG],error_term1[NUMOBSER][MAXLOG],error_term2[NUMOBSER][MAXLOG];
extern char read_line[1024];
