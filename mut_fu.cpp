

#include <string.h>
#include <time.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif

//extern FILE *fRezult;

extern vector < XROSOMA > vecDNK;
vector < SAMPLE >  vecSAMPL;

const char mutaFiles[NO_CANCER_ID][32] = {
    "snv_mnv_BLCA-US.tsv", "snv_mnv_BRCA-US.tsv", "snv_mnv_CESC-US.tsv",
    "snv_mnv_HNSC-US.tsv", "snv_mnv_LUAD-US.tsv", "snv_mnv_LUSC-US.tsv"
};

int getMutApoPart(FILE *f_Mutas);
int skipMutationPart(FILE *f_Mutas, string &sample, int XroID);
void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT,
                string &samp );

//void testfseekOP (FILE *f_Mutas);

////////////////////////////////////////////////////////////////////////

bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 )
{
    return x1.nucPos < x2.nucPos;
}
/////////////////////////////////////////////////////////////////////////////////

int loadMutation( int indCancer )
{
		int ch;
		int cnt1, 
            cntMut=0;
        char buff[1024];
        clock_t start = clock();
    
    sprintf(buff, "%s%s", FOLDER_MUT_DATA, mutaFiles[indCancer]);
    printf("\nLoadMutation from %s\n", buff);
    
	FILE *f_Mutas=fopen( buff, "rb");
    if ( f_Mutas==NULL ) {
			printf ("\nMutations_File =""%s"" not found\n", buff);
			return -1;
    }
	
    while ( (ch=fgetc( f_Mutas )) != 0x0A && ( feof( f_Mutas )) ==0 )
        { }	//skip header
		
		while ( ( cnt1=getMutApoPart(f_Mutas )) > 0 )
			cntMut += cnt1;
		
		if ( cnt1 < 0 )
			return cnt1;
			
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\tloaded [%d] dT=%5.2f sec\n", cntMut, duration);

    return cntMut;
}
///////////////////////////////////////////////////////

int getMutApoPart(FILE *f_Mutas)
{ 
		long int f_Pos=ftell(f_Mutas);
		char cBuffer[4096];
		char chREF, chALT;
		string sample;
		int XroID;
		int NucPos;
		int nXro, nSamp;
        SAMPLE *pSam;
		int retC;   // = noMutations

						
	if ( fgets_ShortRec(cBuffer, sizeof(cBuffer), f_Mutas) == 0 )
        return 0;
    
    while  ( cBuffer[0]=='/' && cBuffer[1]=='/')    {
        if ( fgets_ShortRec(cBuffer, sizeof(cBuffer), f_Mutas) == 0 )
            return 0;
    }
	parsMutRec( cBuffer,  &XroID, &NucPos,  &chREF, &chALT, sample);

    for ( nSamp=0; nSamp<vecSAMPL.size(); nSamp++ )
        if ( sample.compare(vecSAMPL[nSamp].SampName)==0 )
            break;
	if ( nSamp >= vecSAMPL.size() )
        vecSAMPL.push_back(SAMPLE(sample));
    pSam = &vecSAMPL[nSamp];
    
    nXro = findXroByID( XroID );
	if ( nXro < 0 )	{
		retC = 1 + skipMutationPart(f_Mutas, sample, XroID);
		printf("\tINV.chrID '%.16s ...' >>> Ignored %d Recs\n", cBuffer, retC);
        return	retC;
	}
    if ( pSam->vMutAPO[nXro].size() > 0 ) {
        printf("\nBAD MUTATION ORDER: CHRO_ID=%d  sample=%s\n",
                vecDNK[nXro].chrNum, sample.c_str ( ));
        return (-1);
	}
    
    vecDNK[nXro].testValidDNK(NucPos, chREF, chALT );
    
    if ( vecDNK[nXro].APOtest(NucPos, chREF, chALT) > 0 )
        if ( pSam->vMutAPO[nXro].empty() || pSam->vMutAPO[nXro].back().nucPos <= NucPos )
            pSam->vMutAPO[nXro].push_back( MUTANT(NucPos, chREF, chALT ));
        else {
            vector <MUTANT>::iterator Iter =
            lower_bound( pSam->vMutAPO[nXro].begin( ), pSam->vMutAPO[nXro].end( ),
                        MUTANT(NucPos),  lesser_MUT);
            pSam->vMutAPO[nXro].insert( Iter, MUTANT(NucPos, chREF, chALT ) );
        }
    
	f_Pos=ftell(f_Mutas);

	while ( 1 ) {
		if ( fgets_ShortRec(cBuffer, sizeof(cBuffer), f_Mutas )==0 )
			break;
		int pars_xroNum;
		string pars_samp;
        if ( cBuffer[0]=='/' && cBuffer[1]=='/')
            continue;
		parsMutRec( cBuffer, &pars_xroNum, &NucPos,  &chREF, &chALT, pars_samp);
		if ( pars_xroNum != XroID || pars_samp.compare(sample)!= 0 )	{		
			fseek(f_Mutas, f_Pos, SEEK_SET);
			break;
		}
		vecDNK[nXro].testValidDNK(NucPos, chREF, chALT );
        
        if ( vecDNK[nXro].APOtest(NucPos, chREF, chALT) > 0 )
            if ( pSam->vMutAPO[nXro].empty() || pSam->vMutAPO[nXro].back().nucPos <= NucPos )
                pSam->vMutAPO[nXro].push_back( MUTANT(NucPos, chREF, chALT ));
            else {
                vector <MUTANT>::iterator Iter =
                lower_bound( pSam->vMutAPO[nXro].begin( ), pSam->vMutAPO[nXro].end( ),
                            MUTANT(NucPos),  lesser_MUT);
                pSam->vMutAPO[nXro].insert( Iter, MUTANT(NucPos, chREF, chALT ) );
            }

		f_Pos = ftell(f_Mutas);
	}
    printf("X.%d #S=%d\tnAPOmut=%d\n",
           vecDNK[nXro].chrNum, (int)(pSam-&vecSAMPL[0]), (int)pSam->vMutAPO[nXro].size() );
	return (int)pSam->vMutAPO[nXro].size();
}
//////////////////////////////////////////////////////////

void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT, string &samp )
{

	if ( *pB == 'X' || *pB=='x' )	
		*xrNum = 23;
	else
	if ( *pB == 'Y' || *pB=='y' )
		*xrNum = 24;
	else
		*xrNum = atoi( pB );

	while ( *pB && *pB!='\t' ) pB++;	//skip CHROM
	if (*pB)	pB++;

	*Xpos = atoi( pB );								//POS
	while ( *pB && *pB!='\t' ) pB++;	//skip POS
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip ID
	if (*pB)	pB++;

	*chREF = *pB;											//REF
	while ( *pB && *pB!='\t' ) pB++;	//skip REF
	if (*pB)	pB++;

	*chALT = *pB;											//ALT
	while ( *pB && *pB!='\t' ) pB++;	//skip ALT
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip QUAL
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip FILTER
	if (*pB)	pB++;
	
	while ( *pB && *pB!='\t' ) pB++;	//skip INFO
	if (*pB)	pB++;

	char *pBend = pB;
	while ( *pBend && *pBend!='\t' ) pBend ++;
	*pBend = '\0';
	samp = pB;

	return;
}
///////////////////////////////////////////////////////

int skipMutationPart(FILE *f_Mutas, string &sample, int XroID)
{
		long int f_Pos=ftell(f_Mutas);
		char cBuffer[4096];
		int retC = 0;

	while ( 1 ) {
		if ( fgets(cBuffer, sizeof(cBuffer), f_Mutas )==NULL )
			break;
		int pars_xroNum;
		string pars_samp;
		int NucPos;
		char chREF, chALT;
		parsMutRec( cBuffer, &pars_xroNum, &NucPos,  &chREF, &chALT, pars_samp);
		if ( pars_xroNum != XroID || pars_samp.compare(sample)!= 0 )	{
			fseek(f_Mutas, f_Pos, SEEK_SET);
			return (retC);
		}
		retC++;
		f_Pos = ftell(f_Mutas);
	}

	return retC;
}

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
/*
void SAMPLE::makeRndCycle( int nS, const char *canc_id)//f_in)
{
	int nX;

#ifdef RANDOM_TRACE
//		makeRezultName(f_in, nS, fName, "?");
//    fName += ".txt");
    sprintf(buffr, "%s%s_smp%2d.txt", FOLDER_OUT_DATA, canc_id, nS);
    FILE *fRandTrace = fopen( buffr, "w");
    
        fprintf(fRandTrace, "%s\nnoMut:\t%d\t", SampName.c_str(), sumMut_All );
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ )    {
            if ( nX >= vecDNK.size() )
                break;
            fprintf(fRandTrace, "\t%d\t", mutX[nX].noMut );
        }
        fprintf(fRandTrace, "\n");
    
        fprintf(fRandTrace, "\nseqNo\tSUM_IN\tcntLess");
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ )	{
            if ( nX >= vecDNK.size() )
                break;
            fprintf(fRandTrace, "\tX.%d_INrang\t%d_cntL",
                                vecDNK[nX].chrNum, vecDNK[nX].chrNum );
        }
        fprintf(fRandTrace, "\n");
    
		fprintf(fRandTrace, "REAL\t%d\t",sumREAL_IN);
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ ) {
            if ( nX >= vecDNK.size() )
                break;
            fprintf(fRandTrace, "\t%d\t", mutX[nX].noREAL_IN);
        }
        fprintf(fRandTrace, "\n\n");
#endif
//   .........................................................
//
    printf("Samp_#%d [nMu=%d : apo=%d]...\t", nS+1, sumXr.nMut, sumXrAPO.nMut);
    clock_t start = clock();
    
    unsigned int nM;
    unsigned long rnd_pos;
    MUTATION *pmutX, *pmutAPO;
    
#if ! defined ( RUN_WITH_RANDOM_CYCLES )
    clearRANDstat( );
#endif
    
    for ( int nCycl=0; nCycl<NO_RANDOM_CYCLES; nCycl++ )	{
        clearRANDstat( );
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ )	{
            vector <XROSOMA> :: pointer pXRO = &vecDNK[nX];
            
            pmutX = &muXr[nX];
            char *pPos = pXRO->pTagBody;
            long sizeBody = pXRO->TagSize;
            for ( nM=0; nM<pmutX->nMut; nM++ )	{
                rnd_pos = ( (double)rand() / (double)RAND_MAX ) * sizeBody;
                if ( GET_MUZONE_TAG(pPos+rnd_pos) )
                    pmutX->nRAND_IN += 1;
                else {
                    pmutX->nRAND_OUTcur++;
                    if ( ! GET_MUZTWO_TAG(pPos+rnd_pos) )
                        pmutX->nRAND_OUTall++;
                }
            } //of noMut
            sumXr.nRAND_IN     += pmutX->nRAND_IN;
            sumXr.nRAND_OUTcur += pmutX->nRAND_OUTcur;
            sumXr.nRAND_OUTall += pmutX->nRAND_OUTall;
            
            pmutAPO = &muXrAPO[nX];
            pPos = pXRO->TCXtag;
            sizeBody = pXRO->TCXsize;
            for ( nM=0; nM<pmutAPO->nMut; nM++ )    {
                rnd_pos = ( (double)rand() / (double)RAND_MAX ) * sizeBody;
                if ( GET_MUZONE_TAG(pPos+rnd_pos) )
                    pmutAPO->nRAND_IN += 1;
                else {
                    pmutAPO->nRAND_OUTcur++;
                    if ( ! GET_MUZTWO_TAG(pPos+rnd_pos) )
                        pmutAPO->nRAND_OUTall++;
                }
            } //of noMut_APO
            sumXrAPO.nRAND_IN     += pmutAPO->nRAND_IN;
            sumXrAPO.nRAND_OUTcur += pmutAPO->nRAND_OUTcur;
            sumXrAPO.nRAND_OUTall += pmutAPO->nRAND_OUTall;
            
            if ( pmutX->nREAL_IN >0 && pmutX->nRAND_IN <= pmutX->nREAL_IN )
                    pmutX->cntLessReal += 1;
            if ( pmutAPO->nREAL_IN >0 && pmutAPO->nRAND_IN <= pmutAPO->nREAL_IN )
                    pmutAPO->cntLessReal += 1;
        }	// of Xromo
			
        if ( sumXr.nREAL_IN >0 && sumXr.nRAND_IN <= sumXr.nREAL_IN )
                sumXr.cntLessReal += 1;
        if ( sumXrAPO.nREAL_IN >0 && sumXrAPO.nRAND_IN <= sumXrAPO.nREAL_IN )
                sumXrAPO.cntLessReal += 1;
//========================================
			
    }	// of limon cycles
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("dT=%5.2f sec\n", duration);

	return;
}
*/
////////////////////////////////////////////////////////

/*
void testfseekOP (FILE *f_Mutas)
{
	char cBuffer[4096];
	int XroID;
	long posBefore, posAfter;
	long int f_pos[10];
	int cnt=2;

	fgets(cBuffer, sizeof(cBuffer), f_Mutas );
	XroID =atoi(cBuffer);
	fprintf(tsld, "%d.\t0\t->\t%d\n", cnt, XroID);

	while ( 1 ) {
		cnt++;
		posBefore = ftell(f_Mutas);
		if ( fgets(cBuffer, sizeof(cBuffer), f_Mutas ) == NULL )
			break;

		int nextID = atoi(cBuffer); 
		if ( nextID != XroID ) {
			fprintf(tsld, "%d.\t%d\t->\t%d\n", cnt, XroID, nextID);
			fseek(f_Mutas, posBefore, SEEK_SET);
			cnt--;
			XroID = nextID;
		}
 
	}
	return;
}
*/
///////////////////////////////////////////////////////
