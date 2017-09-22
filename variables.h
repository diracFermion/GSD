#define KAPPA 5
#define FRAMES 20001
#define NX 100
#define NY 11
#define LEN 1100
#define NMAX 50000


extern int N,Nb,Nd,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
extern float position[NMAX*3];
extern uint32_t particleID[NMAX];
extern char particleType[3][2];
extern float bendingEner[NMAX];
extern float bondHarmonicEner[NMAX];
