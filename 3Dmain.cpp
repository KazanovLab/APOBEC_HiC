// dnk.cpp : Defines the entry point for the console application.
//
//#include <iostream>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	FILE *tsld=NULL;
#endif

	FILE *fRezult=NULL;
//    FILE *fRezAPO=NULL;

/*
 const char CancerID[NO_CANCER_ID][8] = { "BLCA-US", "BOCA-UK", "BRCA-EU", "BRCA-UK", "BRCA-US",
 "BTCA-SG", "CESC-US", "CLLE-ES", "CMDI-UK", "COAD-US",
 "DLBC-US", "EOPC-DE", "ESAD-UK", "GACA-CN", "GBM-US",
 "HNSC-US", "KICH-US", "KIRC-US", "KIRP-US", "LAML-KR",
 "LGG-US",  "LICA-FR", "LIHC-US", "LINC-JP", "LIRI-JP",
 "LUAD-US", "LUSC-US", "MALY-DE", "MELA-AU", "ORCA-IN",
 "OV-AU",   "OV-US",   "PACA-AU", "PACA-CA", "PAEN-AU",
 "PAEN-IT", "PBCA-DE", "PRAD-CA", "PRAD-UK", "PRAD-US",
 "READ-US", "RECA-EU", "SARC-US", "SKCM-US", "STAD-US",
 "THCA-US", "UCEC-US"
 };
*/

const char CancerID[NO_CANCER_ID][8] = { "BLCA", "BRCA", "CESC", "HNSC", "LUAD", "LUSC"};
//const char ZoneID[NO_ZONES_ID][4] = { "apr", "dr", "gq", "ir", "mr", "str", "z", "cpg"};

long XROSOMA::maxXsize = 0;
extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

//void testBVin ();
void testCNCT ();
void tstLower_Bound();
void tst_L_B();

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) 
//
{
    
		char buffr[1024];
        vector <XROSOMA>::iterator curXRO;
        int indCancer;

	if ( argc < 2 ) {
		printf ("\nUsage: %s  'Cancer_ID'\n", argv[0]);
		return -1;
	}

    for ( indCancer=0; indCancer<NO_CANCER_ID; indCancer++)
        if ( strcmp(argv[1], CancerID[indCancer])==0 )
            break;
    if ( indCancer >=NO_CANCER_ID ) {
        printf ("INV.Cancer_ID = '%s'\n", argv[1]);
        return -1;
    }
    
#ifdef DEBUG_TRACE_TSLD
    sprintf(buffr, "%smutrace.txt", FOLDER_OUT_DATA);
    tsld = fopen(buffr, "w");
#endif

//    testCNCT ();
//    tstLower_Bound();
//    tst_fgets();
//    return 0;
//
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )
		return -1;
    
    if ( LoadContactSet(FOLDER_CONTACT_DATA) < 0 )
        return -1;
    
// printTouchSet();
//    tst_L_B();
//    return 0;
//---------------------------------------------
	if ( loadMutation( indCancer ) < 0 )
		return -1;
/*
    putMutation(0, 0, 50);
    putMutation(0, 1, 50);
    putMutation(0, 2, 50);
    return 0;
*/
    sprintf(buffr, "%sdist_%s.txt", FOLDER_OUT_DATA, argv[1]);
    fRezult = fopen(buffr, "w");
    printTableCap( );
    
    clock_t start = clock();
    for ( int iSamp=0; iSamp<vecSAMPL.size(); iSamp++ ) {
//        fprintf(tsld, "\n*********Sample #%d mutations:\n", iSamp );
//        vecSAMPL[iSamp].tstSampMut(  );
        
        vecSAMPL[iSamp].proSampMut( );
        
    }
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("proSampl dT=%5.2f sec\n", duration);
    
    fclose (fRezult);
    
	return 0;
}
//////////////////////////////////////////////////////////////////////

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in )
{
    char buff[4096];
    int lng;
    
    if ( fgets(shortRec, sizeRec, f_in)==NULL )
        return 0;
    lng = (int)strlen(shortRec);
    if ( *(shortRec+lng-1) != '\n' )
        while ( 1 ) {
            if ( fgets(buff, sizeof(buff), f_in)==NULL )
                break;
            if ( buff[strlen(buff)-1] == '\n' )
                break;
        }
    
    return lng;
}
////////////////////////////////////////////////////////

