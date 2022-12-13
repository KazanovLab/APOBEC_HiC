
#ifndef _MUTAS_H__
#define _MUTAS_H__

#include "vector"
using namespace std;

#define RUN_WITH_RANDOM_CYCLES

//#define DEBUG_TRACE_RTSEG

#define DEBUG_TRACE_TSLD

//#define RANDOM_TRACE

#define TREE_TRACER

#define FOLDER_HUMAN "/Users/gena/mutas/humanXRO/"
//#define FOLDER_MUT_DATA "/Users/gena/mutas/mut_set/"
//#define FOLDER_CONTACT_DATA "/Users/gena/mu3D/IMR90_1mb/"
//#define FOLDER_OUT_DATA "/Users/gena/mu3D/APO3d/outd/"

#define FOLDER_MUT_DATA     "/Users/gena/mu3D/test_1mb/3/"
#define FOLDER_CONTACT_DATA "/Users/gena/mu3D/test_1mb/3/"
#define FOLDER_OUT_DATA     "/Users/gena/mu3D/test_1mb/3/"

#if defined ( RUN_WITH_RANDOM_CYCLES )
#define NO_RANDOM_CYCLES 10000
#else
#define NO_RANDOM_CYCLES 0
#endif

#define XRO_SET_FILENAME "xromo_set.txt"

#define NO_HUMAN_XRO 24
#define SEGM_SIZE 1000000     // _1mb
#define IN_POWER_THRESHOLD 1000.0
#define EX_POWER_THRESHOLD 100.0

#define NO_CANCER_ID 6
//#define NO_ZONES_ID 8


int LoadDNKset( const char *pFolder, const char *fSetName );
int LoadContactSet(const char *pFolder);
int loadMutation( int indCancer );
//int proSampMut( SAMPLE *pSampl );
void printTableCap( );

void printTouchSet(int step);       //test
void putMutation(int iSam, int iXro, int maxPos);

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );
void tst_fgets();


#endif
