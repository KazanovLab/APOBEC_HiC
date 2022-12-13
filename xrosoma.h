#ifndef _XROSOMA_H__
#define _XROSOMA_H__


#include "vector"
#include "string"
using namespace std;

#define _C 0
#define _G 1
#define _A 2
#define _T 3

//////////////////////////////////////////////////////////

struct MUTANT
{
		long nucPos;
		char nucREF;
		char nucALT;
		MUTANT() {nucPos=-1; nucREF='?'; nucALT='?';};
        MUTANT(long P) {nucPos=P; };
		MUTANT(long P, char R, char A) {nucPos=P; nucREF=R; nucALT=A;};
};
bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 );

//////////////////////////////////////////////////////////

struct TOUCH_SET {
    int ownTAR;                       //  own touch area
    vector < pair <int, int> > vFork;    // < indXro2, foreignTAR  >
    TOUCH_SET( int c1 ) { ownTAR=c1; vFork.clear(); };
    TOUCH_SET( int x1, int x2, int tar) { ownTAR=x1; vFork.clear();
        vFork.push_back(pair <int, int>(x2, tar) ); };
};
bool lesser_TOUCH ( const TOUCH_SET &x1, const TOUCH_SET &x2 );
bool lesser_PAIR ( const pair <int, int> &x1, const pair <int, int> &x2 );

//////////////////////////////////////////////////////////

class SAMPLE {
public:
    string SampName;
    vector < MUTANT > vMutAPO[NO_HUMAN_XRO]; // only APOBEC set mutarions in XRO's
    
    SAMPLE(string S) { SampName=S; };
    
    int proSampMut( );
    int tstSampMut(  );
    int cnt_ALL_Mut(int ixXr, int ixMu, vector <TOUCH_SET>::iterator ItRoot);
    int cntMutInSeg(int ixXr, int ixMu);
    int cntMutInRoot(vector <TOUCH_SET>::iterator ItRoot);
    int cntMutInRoot(vector <TOUCH_SET>::iterator ItRoot, vector<pair<int, int>> vMinS);
    int getNearMutUP(int iXro, vector <TOUCH_SET>::iterator ItRoot); //, int iMut); //, int OUTway);
    int getNearMutDOWN(int iXro, vector <TOUCH_SET>::iterator ItRoot); //, int iMut); //, int OUTway);
    void printMutLines( int iXro, int iCurMu, const char *comment );

};

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

class  XROSOMA {
private:
//		vector <HUGEN> :: iterator begSearch;
public:
    static long maxXsize;
	string XfName;
	string Xref_NC;
	int  chrNum;				// get from Xref_NC
	string version;
//
//  ===== real body of XRO will be read from files =============================
    long Xsize;
    char *Xbody;
    char *Xtag;     //
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
/*  ===== Targets (TCx, xGA) for APOBEC mutations  =====================
    long TCxMemSize;    // = sizeof (TCXtag)
    char *TCXtag;       //contains tags from '*Xtag' with TCx_TAG
                        //                      for current segment
    long TCXsize;      // = real no Targets in *TCXtag  (<=TCxMemSize)
*/
    vector < int > vTCxRef;
    int makeTCxRef( );
    vector < TOUCH_SET >  vTouchTree;  // contacts set only with power_threshold
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

///////////
    XROSOMA() { chrNum=0; Xsize=0L; Xbody=NULL; };
    XROSOMA(int cN) { chrNum=cN; Xsize=0L; Xbody=NULL; };
    

	long read_DNK( const char *pFolder, char *fName );
    long read_CNTCs(const char *pFolder, int indXro2);
    int  makeNewTouch (int ownTAR, int extTAR, int indXro2);
    int TreeUP( int Rlevel, vector <TOUCH_SET>::iterator ItRoot, int iNextMut, int OUTway );
    int TreeDOWN( int Rlevel, vector <TOUCH_SET>::iterator ItRoot, int iNextMut, int OUTway );
    int addLock(int iXro, int segRoot, int dirMove, int OUTway);
//    int isIt_bLoop(int iXro, int segRoot, int segRoot_1);

    int testValidDNK( int Pos, char chREF, char chALT );
    int APOtest( long Pos, char chREF, char chALT );
};
/////////////////////////////////////////////////////////

//#define NUCID(p_XBODY)         ( *(p_XBODY) & 0x7F )
#define MUZONE_TAG  0x80
//#define TCx_TAG     0x40
#define APO_TAG     0x20        // MUT_TAG
#define MUZTWO_TAG  0x08        // for save all zone tags together

#define SET_MUZONE_TAG(p_TAG)    *(p_TAG) |= MUZONE_TAG
#define GET_MUZONE_TAG(p_TAG)  (((*(p_TAG) & MUZONE_TAG)==0) ? 0 : 1 )

#define SET_MUZTWO_TAG(p_TAG)    *(p_TAG) |= MUZTWO_TAG
#define GET_MUZTWO_TAG(p_TAG)  (((*(p_TAG) & MUZTWO_TAG)==0) ? 0 : 1 )

//#define SET_TCx_TAG(p_TAG)    *(p_TAG) |= TCx_TAG
//#define GET_TCx_TAG(p_TAG)  (((*(p_TAG) & TCx_TAG)==0) ? 0 : 1 )

#define SET_APO_TAG(p_TAG)    *(p_TAG) |= APO_TAG
#define GET_APO_TAG(p_TAG)  (((*(p_TAG) & APO_TAG)==0) ? 0 : 1 )


//void  switchToXRO( );    // pointXRO( ); //
//void  switchToRT( );     // pointRT( );  //

///////////////////////////////////////////////

int findXroByID( int xID );
int getNucID(const char Nuc);
//int CNCT_Dubl( vector <CNCT_SET>::iterator ItCN, CNCT_SET cnct,
//              vector <CNCT_SET>::iterator ItLast);

void printTableCap( FILE *fRez );

#endif
