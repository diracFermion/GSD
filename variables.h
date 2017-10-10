#define KAPPA 5.0
#define STEPS 100000000
#define PERIOD 10000
#define FRAMES STEPS/PERIOD
#define NX 101
#define NY 51
#define LEN NX*NY
#define NMAX 50000
#define EPSILON 3600.0
#define a 1.0
#define RUN 10

extern int N,Nb,Nd,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
//N:#particles, Nb:#bonds, Nd:#dihedrals
extern float position[NMAX*3];
extern uint32_t particleID[NMAX];
extern char particleType[3][2];
extern double bendingEner[NMAX];
extern double bondHarmonicEner[NMAX];
extern double total_DHE,total_BHE;
extern double hgt_fluctuation[NMAX];
extern double h_avg_node[NMAX];
extern double hgt_fluctuation[NMAX];
extern double h_width[FRAMES/2][2*NX];
