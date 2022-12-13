


#include <string.h>
#include <time.h>
#include <ctype.h>
#include <algorithm>
#include <functional>      // For greater<int>( )

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif

vector < XROSOMA > vecDNK;

char complment[4][2] = { {'C','G'},  {'G','C'}, {'A','T'}, {'T','A'} };
char NucID[5] = "CGAT";
int  cmpInd[4] = { _G, _C, _T, _A };

///////////////////////////////////////////////////////

int LoadDNKset( const char *pFolder, const char *fSetName )
{
    FILE *fiSet=NULL;
    char fBuff[1024];
    int    noDNK = 0;
    long   noTcx;
    
    sprintf(fBuff, "%s%s", pFolder, fSetName);
    printf("\nLoadDNK from list '%s'.....\n", fBuff);
    clock_t start = clock();
    
    if ( ( fiSet=fopen(fBuff, "r"))==NULL )    {
        printf("File '%s' : ERR_OPENop\n",fBuff);
        return -1;
    }
    
    while ( 1 )    {
        if( fgets(fBuff, sizeof(fBuff), fiSet) == NULL )
            break;
        if ( fBuff[0]=='/' )
            continue;
        char *pn = strchr(fBuff, '\r');
        if ( pn )
            *pn = '\0';
        else    {
            pn = strchr(fBuff, '\n');
            if ( pn )    *pn = '\0';
        }
        vecDNK.push_back(XROSOMA());
        if ( (noTcx=vecDNK.back().read_DNK(pFolder, fBuff)) < 0 )    {
            vecDNK.pop_back( );
            continue;
        }
        noDNK++;
        printf("\tX.%02d [%ld]\tNoTcx=%ld\n",
               vecDNK.back().chrNum, vecDNK.back().Xsize, noTcx);
    }
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("===endLoadDNK=[%d]\tdT=%5.2f sec\n", noDNK, duration);
    
    fclose(fiSet);
    
    if ( noDNK != NO_HUMAN_XRO ) {
        printf("\n Mismatch NO_HUMAN_XRO=%d :: Must be = %d\n", noDNK, NO_HUMAN_XRO);
        return -1;
    }
    
    return noDNK;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

long XROSOMA::read_DNK( const char *pFolder, char *fName )
{
		FILE *fileDnk;
        char fPath[1024];

    sprintf(fPath, "%s%s", pFolder, fName);
	if ( !( fileDnk=fopen(fPath, "r")) )	{
		printf("File '%s': ERR_OPEN_op\n",fPath);
		return -1;
	}
	
	XfName.assign( fPath );

	fseek(fileDnk, 0, SEEK_END); 
	Xsize = ftell(fileDnk); 
    Xbody =  new char[Xsize+1];//+3];
    if ( Xsize > maxXsize )
        maxXsize = Xsize;
    

	fseek(fileDnk, 0, SEEK_SET);
		
	fgets(Xbody, Xsize, fileDnk );  // header :

	char *pn = strstr(Xbody,"ref");
	if ( pn )
		pn = strstr(Xbody,"NC_");

	if ( !pn || Xbody[0]!='>' ) {
		printf("File '%s': not Found 1_st control Rec '%40s'\n",fPath, Xbody);
		return -1;
	}

	pn+=3;
	while( *pn && (*pn)!= '|' )  Xref_NC += *pn++;
	chrNum = atoi( Xref_NC.c_str() );	
	Xref_NC.insert(0,"NC_");

	char *pbody = Xbody;
	long count= Xsize;
	
	while ( 1 )	{
        if ( fgets(pbody, count, fileDnk) == NULL )
				break;
            if ( (pn = strchr(pbody, '\n')) )
					*pn = '\0';
			while ( *pbody ) pbody++;
			count = Xsize - (pbody-Xbody);
    }
	*(pbody+1) = '\0';
    
    Xtag  =  new char[Xsize+1];
    memset((void *)Xtag, '\0', Xsize+1);
    int noTCx=0;
//    noTcx = markTCX_TAG ( );       ///// TCx
    
    noTCx = makeTCxRef( );
	
    fclose(fileDnk);

    return noTCx;
}
////////////////////////////////////////////////////////////////////////////////////

int XROSOMA::makeTCxRef( )
{
    char *pStop = Xbody + Xsize -1;
    int pos = 0;
    
    for ( char *pXb=Xbody; pXb<=pStop; pXb++, pos++)  {
        if ( (*pXb=='C' && *(pXb-1)=='T') || (*pXb=='G' && *(pXb+1)=='A')  )    {
            vTCxRef.push_back(pos);
        }
    }
    
    return (int) vTCxRef.size();
}
/////////////////////////////////////////////////////////////////////////////

int XROSOMA::APOtest( long Pos, char chREF, char chALT )
{
    // APOBEC:
    //    a) TCx ---> 'C' > ['G' | 'T' ]
    //    b) xGA ---> 'G' > ['C' | 'A' ]
    // Returns: 1=APOBEC; 2=APOBEC in ZONE; 0= any other
    //
    
    char *pB = Xbody+Pos-1;
    int retC=0;
    
    if (  *pB != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%d[%ld]='%c';  MUTATION: '%c' > '%c'\n",
               chrNum, Pos, *pB, chREF, chALT);
        return 0;
    }
    
    switch ( chREF ) {
        case 'C':
            if ( ! (chALT == 'G' || chALT == 'T') )
                return 0;
            if ( *(pB-1) != 'T' )
                return 0;
            retC = 1;
            break;
        case 'G':
            if ( ! (chALT == 'C' || chALT == 'A') )
                return 0;
            if ( *(pB+1) != 'A' )
                return 0;
            retC = 1;
            break;
        default:
            return 0;
    }
    
    return retC;
}
////////////////////////////////////////////////////////////////////////////

int XROSOMA::testValidDNK( int Pos, char chREF, char chALT )
{
    //    char *pB = Xbody+Pos-1;
    
    if (  *(Xbody+Pos-1) != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%d[%d]='%c';  MUTATION: '%c' > '%c'\n",
               chrNum,Pos,*(Xbody+Pos-1), chREF, chALT);
        return 0;
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////////

int  findXroByID( int xID )
{
    int indX;
    
    for ( indX=0; indX<vecDNK.size(); indX++ )
        if ( vecDNK[indX].chrNum == xID )
            return indX;
    printf("ERR: XroID=%d NOT FOUND\n", xID);
    return -1;
}
////////////////////////////////////////////////////// //////////////////////////////

int getNucID ( const char Nuc )
{
    //    char *pp = strchr(NucID, Nuc);
    //    return ( ( !pp ) ? -1 : (int)(pp-NucID) );
    switch (Nuc) {
        case 'C':
            return _C;
        case 'G':
            return _G;
        case 'A':
            return _A;
        case 'T':
            return _T;
        default:
            break;
    }
    return -1;
}
////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////

bool lesser_TOUCH ( const TOUCH_SET &x1, const TOUCH_SET &x2 )
{
    
 return x1.ownTAR < x2.ownTAR;
    
};
//------------------------------
bool lesser_PAIR ( const pair <int, int> &x1, const pair <int, int> &x2 )
{
    
    return x1.first < x2.first;
    
};
///////////////////////////////////////////////////////

int LoadContactSet(const char *pFolder)
{
    char fPath[1024];
    int indXr;

//  ===========================  load INTRA contacts  ==============================
    for ( int IdX1=1; IdX1 <= NO_HUMAN_XRO; IdX1++ )    {
        sprintf(fPath, "%schr%d_chr%d_1mb.txt", pFolder, IdX1, IdX1);
        indXr = findXroByID(IdX1);
        vecDNK[indXr].read_CNTCs(fPath, -1);
    }
//    printTouchSet(0);
//  ===========================  load EXTER contacts  ==============================
    for ( int IdX1=1; IdX1 <= NO_HUMAN_XRO; IdX1++ )    {
        for ( int IdX2=IdX1+1; IdX2 <= NO_HUMAN_XRO; IdX2++ )   {
            sprintf(fPath, "%schr%d_chr%d_1mb.txt", pFolder, IdX1, IdX2);
            indXr = findXroByID(IdX1);
            vecDNK[indXr].read_CNTCs(fPath, IdX2);
        }
    }
//    printTouchSet(1);
    return 1;
}
///////////////////////////////////////////////////////
 
long XROSOMA::read_CNTCs(const char *fPath, int IdXro2)
{

    FILE *fiCnct=NULL;
    char Buffr[1024];
    vector <TOUCH_SET>::iterator ItC;

    int iXr2 =       ( IdXro2 <= 0 ) ? (int)(this-&vecDNK[0]) : findXroByID(IdXro2);
    int PowerThres = ( IdXro2 <= 0 ) ? IN_POWER_THRESHOLD : EX_POWER_THRESHOLD;
    
    if ( !( fiCnct=fopen(fPath, "r")) )    {
        printf("File '%s': ERR_OPEN_op\n",fPath);
        return -1;
    }

    while ( 1 ) {
        if ( fgets(Buffr, sizeof(Buffr), fiCnct) == NULL )
            break;
        if ( Buffr[0]=='/' && Buffr[1]=='/')
            continue;
        float power;
        int ownT, extT;
         if ( sscanf(Buffr, "%d%d%f", &ownT, &extT, &power) != 3 )    {
            printf("File '%s': ERR_Rec_FRMT='%s'\n", fPath, Buffr);
            continue;
        }
        if ( (power < PowerThres) || (IdXro2 <= 0 && ownT==extT)  )
            continue;
        
        ownT /= SEGM_SIZE;
        extT /= SEGM_SIZE;
        makeNewTouch (ownT, extT, iXr2);
        if ( IdXro2 <= 0 )  {                       // for INTRA contact
            makeNewTouch (extT, ownT, iXr2);        // when (ownT >> extT) then (extT >> ownT)
            continue;
        }
// -----EXTERN contact  and now insert mirror [chr2 chr1]
        vecDNK[iXr2].makeNewTouch (extT, ownT, (int)(this-&vecDNK[0]));
    }
    fclose(fiCnct);

    return 0;
}
///////////////////////////////////////////////////////////////////

int XROSOMA:: makeNewTouch (int ownT, int extT, int indXro2)
{
    vector <TOUCH_SET>::iterator ItTset;    //ItTar1
    vector < pair <int, int> >::iterator ItFork;
    
    if ( vTouchTree.empty() || vTouchTree.back().ownTAR < ownT )    {
        vTouchTree.push_back( TOUCH_SET(ownT) );
        ItTset = vTouchTree.end() - 1;
    }
    else {
        ItTset = lower_bound(vTouchTree.begin(), vTouchTree.end(), TOUCH_SET(ownT), lesser_TOUCH);
        if ( ItTset->ownTAR != ownT )
            ItTset = vTouchTree.insert(ItTset, TOUCH_SET(ownT) );
    }
    
//  ---- now fork_set
    if ( ItTset->vFork.empty() )
        ItTset->vFork.push_back( pair <int, int> (indXro2, extT) );
    else {
        for ( ItFork=ItTset->vFork.begin(); ItFork != ItTset->vFork.end(); ItFork++) {
            if ( ItFork->first == indXro2 && ItFork->second == extT )
                break;
        }
        if ( ItFork == ItTset->vFork.end() )
            ItTset->vFork.push_back(  pair <int, int> (indXro2, extT) );
        else
            printf("DUBL_Touch: X.%d:%d X.%d:%d\n",
                   this->chrNum,  ownT, vecDNK[indXro2].chrNum, extT);
    }


    return 1;
}
////////////////////////////////////////////////////////////////////////////////////
/*
int CNCT_Dubl( vector <CNCT_SET>::iterator ItCN, CNCT_SET cnct,
              vector <CNCT_SET>::iterator ItLast)
{
    vector <CNCT_SET>::iterator It;
    if ( ItCN->nCnct1 != cnct.nCnct1 )
        return 0;
    for ( It=ItCN; It != ItLast; It++ ) {
        if ( It->nCnct1 != cnct.nCnct1)
            break;
        if ( It->nCnct1==cnct.nCnct1 && It->nCnct2==cnct.nCnct2 && It->indXr2==cnct.indXr2 )
            return 1;
    }
    
    return 0;
}
 */
////////////////////////////////////////////////////////////////////////////////////

