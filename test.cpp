//
//  test.cpp
//  APO3d
//
//  Created by Ирина Пономарева on 31/07/2022.
//  Copyright © 2022 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>

#include "mutas.h"
#include "xrosoma.h"

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE > vecSAMPL;

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

/////////////////////////////////////////////////////////////////

void testCNCT ()
{
    
    vector < TOUCH_SET > vTouchTree;
    vector <TOUCH_SET>::iterator ItTar1;
    vector < pair <int, int> >::iterator ItFork;
    int Cin[15][3] = { {1,0,0}, {4,0,0}, {6,0,0}, {1,1,0}, {3,0,0},
        {1,1,1}, {4,0,1}, {7,0,0}, {7,1,0}, {7,0,1},
        
        {1,1,1}, {1,1,0}, {7,1,0}, {7,0,0}, {7,0,1}
    };
    vector < pair <int, int> > fork;    // < tpX2, indXr2 >
    
    for ( int n=0; n<15; n++ )  {
        if ( vTouchTree.empty() || vTouchTree.back().ownTAR < Cin[n][0] )    {
            vTouchTree.push_back( TOUCH_SET(Cin[n][0]) );
            ItTar1 = vTouchTree.end() - 1;
        }
        else {
            ItTar1 = lower_bound(vTouchTree.begin(), vTouchTree.end(), Cin[n][0], lesser_TOUCH);
            if ( ItTar1->ownTAR != Cin[n][0] )
                ItTar1 = vTouchTree.insert(ItTar1, TOUCH_SET(Cin[n][0]) );
        }
        //  ---- now fork_set
        if ( ItTar1->vFork.empty() )
            ItTar1->vFork.push_back(  pair <int, int> (Cin[n][1], Cin[n][2]) );
        else {
            for ( ItFork=ItTar1->vFork.begin(); ItFork != ItTar1->vFork.end(); ItFork++)
                if ( ItFork->first == Cin[n][1] && ItFork->second == Cin[n][2])
                    break;
            if ( ItFork == ItTar1->vFork.end() )
                ItTar1->vFork.push_back(  pair <int, int> (Cin[n][1], Cin[n][2]) );
            else
                printf("DUBL_Touch: %d %d %d\n", Cin[n][0], Cin[n][1], Cin[n][2]);
            
        }
    }
    for ( int n=0; n<vTouchTree.size(); n++ )   {
        printf("%d\t", vTouchTree[n].ownTAR );
        for ( int f=0; f<vTouchTree[n].vFork.size(); f++ )
            printf("[%d:%d] ", vTouchTree[n].vFork[f].first, vTouchTree[n].vFork[f].second);
        printf("\n");
    }
    printf("------\n");
    vTouchTree.clear();
    
    for ( int n=14; n>=0; n-- )  {
        if ( vTouchTree.empty() || vTouchTree.back().ownTAR < Cin[n][0] )    {
            vTouchTree.push_back( TOUCH_SET(Cin[n][0]) );
            ItTar1 = vTouchTree.end() - 1;
        }
        else {
            ItTar1 = lower_bound(vTouchTree.begin(), vTouchTree.end(), Cin[n][0], lesser_TOUCH);
            if ( ItTar1->ownTAR != Cin[n][0] )
                ItTar1 = vTouchTree.insert(ItTar1, TOUCH_SET(Cin[n][0]) );
        }
        //  ---- now fork_set
        if ( ItTar1->vFork.empty() )
            ItTar1->vFork.push_back(  pair <int, int> (Cin[n][1], Cin[n][2]) );
        else {
            for ( ItFork=ItTar1->vFork.begin(); ItFork != ItTar1->vFork.end(); ItFork++)
                if ( ItFork->first == Cin[n][1] && ItFork->second == Cin[n][2])
                    break;
            if ( ItFork == ItTar1->vFork.end() )
                ItTar1->vFork.push_back(  pair <int, int> (Cin[n][1], Cin[n][2]) );
            else
                printf("DUBL_Touch: %d %d %d\n", Cin[n][0], Cin[n][1], Cin[n][2]);
            
        }
    }
    for ( int n=0; n<vTouchTree.size(); n++ )   {
        printf("%d\t", vTouchTree[n].ownTAR );
        for ( int f=0; f<vTouchTree[n].vFork.size(); f++ )
            printf("[%d:%d] ", vTouchTree[n].vFork[f].first, vTouchTree[n].vFork[f].second);
        printf("\n");
    }
    printf("------\n");
    vTouchTree.clear();
    
    
}
////////////////////////////////////////////////////////

void printTouchSet(int step)
{
    XROSOMA *pXro;
    FILE *fTs;
    char fName[1024];
    
    if ( step==0)   {
        sprintf(fName, "%stTouch_OWN.txt", FOLDER_OUT_DATA);
        fTs = fopen(fName, "w");
        fprintf(fTs,"=========OWN touches ========================\n");
    }
    else    {
        sprintf(fName, "%stTouch_ALL.txt", FOLDER_OUT_DATA);
        fTs = fopen(fName, "w");
        fprintf(fTs,"\n=========ALL touches ========================\n");
    }
    
    for ( int nXr=0; nXr < NO_HUMAN_XRO; nXr++ ) {
        pXro = &vecDNK[nXr];
        fprintf(fTs,"X.%d [%d]\t=> %d ===========\n", pXro->chrNum,
                (int)(pXro->Xsize/SEGM_SIZE), (int)pXro->vTouchTree.size() );
        
        for ( int nTs=0; nTs < pXro->vTouchTree.size(); nTs++ )
            fprintf(fTs,"  #%d\t=> %d\n",
                    pXro->vTouchTree[nTs].ownTAR,
                    (int)pXro->vTouchTree[nTs].vFork.size() );
        fprintf(fTs,"======\n");
    }
    
    fclose (fTs);
    
    return;
}
//////////////////////////////////////////////////////////////////////

void tstLower_Bound()
{
    
    int n, m;
    vector < TOUCH_SET >  vT_Tree;
    int token[10] = { 0, 5, 10, 11, 19, 40, 43, 49, 50, 55};
    vector < TOUCH_SET >::iterator itLow, itUpp;
    char cUp[8], cLow[8];
    
    for ( n=0; n<=5; n++ )  {
        vT_Tree.push_back( TOUCH_SET( (n*10) ) );
        vT_Tree[n].vFork.push_back( pair<int,int>(n+1,n+2) );
    }
    
    for ( n=0; n<vT_Tree.size(); n++ )
        fprintf(tsld, "%d ",  vT_Tree[n].ownTAR);
    
    fprintf(tsld, "\ntoken\tLow\tUpper\n");
    for ( n=0; n<10; n++ )  {
        itLow= lower_bound(vT_Tree.begin(), vT_Tree.end(), TOUCH_SET(token[n]), lesser_TOUCH);
        itUpp= upper_bound(vT_Tree.begin(), vT_Tree.end(), TOUCH_SET(token[n]), lesser_TOUCH);
        m=0;
        if ( itUpp==vT_Tree.begin() )
            sprintf(cUp, "Beg");
        else    if ( itUpp==vT_Tree.end() )
                    sprintf(cUp, "End");
                else
                    sprintf(cUp, "%d", itUpp->ownTAR);
        
        if ( itLow==vT_Tree.begin() )
            sprintf(cLow, "Beg");
        else    if ( itLow==vT_Tree.end() )
                    sprintf(cLow, "End");
                else
                    sprintf(cLow, "%d", itLow->ownTAR);

        fprintf(tsld, "%d\t%s\t%s\n",  token[n], cLow, cUp);
        
    }
    fprintf(tsld, "End() - Beg() =  %d\n", (int)(vT_Tree.end() - vT_Tree.begin()));
    fprintf(tsld, "End()-1 =[%d]= %d\n", (int)(vT_Tree.end()-1 - vT_Tree.begin()),
                    (vT_Tree.end()-1)->ownTAR);
    
    return;
}
//////////////////////////////////////////////////////////////////////

void putMutation(int iSam, int iXro, int maxSegm)
{
    FILE *fMut;
    char fName[1024];
    
    sprintf(fName, "%sSMPL%d_X%d.txt", FOLDER_OUT_DATA, iSam, iXro);
    fMut = fopen(fName, "w");
    fprintf(fMut,"CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tsample\t1000genomes_AF\n");
    
    SAMPLE *pSamp = &vecSAMPL[iSam];
    int Pos = (maxSegm+1) * SEGM_SIZE;
    for ( int nM=0; nM<pSamp->vMutAPO[iXro].size(); nM++ ) {
        if ( pSamp->vMutAPO[iXro][nM].nucPos > Pos )
            break;
        fprintf(fMut,"//%d\t%ld\t.\t%c\t%c\t.\t.\t1000geno, \t%s\t.\n",
                iXro+1, pSamp->vMutAPO[iXro][nM].nucPos, pSamp->vMutAPO[iXro][nM].nucREF,
                pSamp->vMutAPO[iXro][nM].nucALT, pSamp->SampName.c_str() );
    }
    fclose (fMut);
    return;
}
//////////////////////////////////////////////////////////////////////

void tst_fgets()
{
    int lng;
    FILE *fMut;
    char fName[1024], *p;
    
    
    sprintf(fName, "%ssnv_mnv_BLCA-US.tsv", FOLDER_OUT_DATA);
    fMut = fopen(fName, "r");
    
    while ( 1 ) {
        if ( fgets(fName, sizeof(fName), fMut)==NULL )
            break;
        p = strchr (fName, '\n');
        lng = (int)strlen(fName);
    }
}
