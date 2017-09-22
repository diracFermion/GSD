#define KAPPA 2.0
#define FRAMES 10
#define NX 101
#define NY 21
#define LEN 1000
#define NMAX 50000
#define EPSILON 1440.0
#define a 1.0


extern int N,Nb,Nd,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
//N:#particles, Nb:#bonds, Nd:#dihedrals
extern float position[NMAX*3];
extern uint32_t particleID[NMAX];
extern char particleType[3][2];
extern double bendingEner[NMAX];
extern double bondHarmonicEner[NMAX];
extern double total_DHE,total_BHE;
