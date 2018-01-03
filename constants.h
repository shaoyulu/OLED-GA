//**************//
// Important   //
//  Constants //
//***********//

// Number of parallel MOPAC
// (one node operation only)
#define MAXPARA 12
// Number of cores per Q-Chem instance (not used)
#define DFTPARA 1
// Number of iterations to wait for queue to finish
#define MAX_TIME_WAIT 25000
#define MAX_BRANCH_WAIT 10000

// Safe Mode Options
#define SKIPMOPAC 0
#define SKIPDFT 1
#define SKIPGSM 1
#define DO_NOT_WAIT 1
#define GEN_SPECIFIC 0

// Screening Thresholds 
#define DFT_EMAX 50     //vs. starting int
#define GSM_EMAX_TS 20  // Ea max = DFT_EMAX + GSM_EMAX
#define PM6_EMAX 80     //was 70
#define MINTSE 10
#define EVSORIGTOL 4.5

// Correction for PM6 overestimating H2 stability
#define H2PENALTY 15

// Connection thresholds for PM6 structures
#define MMVSMOPAC 2 //at least 2

//MOPAC/DFT Settings
#define UNRESTRICTED 0
#define DFTB3LYP 1
#define DFTB97 0
#define DFTwB97 0
#define B631PGSS 0
#define B631GSS 0
#define B631GS  1
#define LANL2DZ 0
#define DFT_OPT_CYCLES 200
#define CATION 1
#define TRIPLET 0
//enter dielectric constant (not implemented)
#define USE_CPCM 0
#define CPCM_DIEL 7.6

