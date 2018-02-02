#define NSAM 32
#define NCHS 48
#define MIN_ADCL -50  // used in getH
#define MAX_ADCL 1000 // used in getH
#define MIN_ADC 50   // used in getP
#define MAX_ADC 750  // used in getP
#define NBINS  256
#define MAX_WIDTH 10 // about merging. Depending on drift velocity: HV and Gas

#define NLAY 9
//#define NLAY  11
//#define NLAY  19

#if NLAY == 9
#define NLAYp2 11
#define NCEL  11
#define NCELA 99
#define NBRD 2
#define NCHT 96
#define NZXP 8
#define NITERSMAX 100
#elif NLAY == 11
#define NLAYp2 13
#define NCEL  12
#define NCELA 132
#define NBRD 3
#define NCHT 144
#define NZXP 10
#define NITERSMAX 100
#elif NLAY == 19
#define NLAYp2 21
#define NCEL  12
#define NCELA 228
#define NBRD 5
#define NCHT 240
#define NZXP 18
#define NITERSMAX 25
#endif

#define NCAND 4
