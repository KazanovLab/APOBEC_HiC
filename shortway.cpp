//
//  shortway.cpp
//  APO3d
//
//  Created by Ирина Пономарева on 25/07/2022.
//  Copyright © 2022 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>
//#include <string.h>

#include "mutas.h"
#include "xrosoma.h"

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;
extern FILE *tsld;
extern FILE *fRezult;

char ntabs[256];        // testTrace
SAMPLE *pCurSmpl;
vector <MUTANT>::iterator ItPrimeMut;
//int iCurMu;   // index current mutation 'vMutNuc[iCurXr][curMut]'

int MinLway;                            //  min Length of way  between 'curMut' & 'nearMut'
vector < pair <int, int> > vMinSet;     // <iNearXr, iNearMu> all case for 'MinLway'
vector < pair <int, int> > LockSet;     // <key,OUTway>  key = +-(segRoot*SEGM_SIZE+iRootXr)
//vector < int > LockSet;              // UP =(segRoot*SEGM_SIZE+iRootXr); DOWN =-(segRoot*1000+iRootXr)

int doublTest(int iXr, int iMut);
/*
struct Trip {
    int root;
    int r_1;
    int iX;
    Trip(int x, int y, int z){ root=x; r_1=y; iX=z;};
};
vector < Trip > bLoop;
/////////////////////////////////////////////////////////////

 int XROSOMA:: isIt_bLoop(int iXro, int segRoot, int segRoot_1)
 {
 //    int key = 1000 * segRoot + iXro;
 for ( int n=0; n<bLoop.size(); n++ )
 if ( bLoop[n].root==segRoot && bLoop[n].r_1==segRoot_1 && bLoop[n].iX==iXro )
 return 1;
 bLoop.push_back( Trip (segRoot, segRoot_1, iXro));
 
 return 0;
 }
 */
/////////////////////////////////////////////////////////////

int XROSOMA:: addLock(int iXro, int segRoot, int dirMove, int OUTway)
{
//    int segRoot= ( ItRoot != vTouchTree.end() && ItRoot >= vTouchTree.begin() )
//                ?  ItRoot->ownTAR : (vTouchTree.end()-1)->ownTAR + 1;
    char dir[8];
    if (dirMove>0)
        strcpy(dir, "UP");
    else strcpy(dir, "DOWN");
    
    int key = dirMove*(SEGM_SIZE * segRoot + iXro);
    for ( int n=0; n<LockSet.size(); n++ )
        if ( LockSet[n].first==key )    {
            if ( LockSet[n].second > OUTway )   {
                fprintf(tsld, "  chng_LOCK_%s=%d:%d->L%d\n", dir, key, LockSet[n].second, OUTway);
                LockSet[n].second = OUTway;
                return n;
            }
            fprintf(tsld, "  %s_LOCKed=%d:%d L%d\n",dir, key, LockSet[n].second, OUTway);
            return -(n+1);
        }
    LockSet.push_back(pair <int, int> (key, OUTway));
    fprintf(tsld, "  add_LOCK_%s=%d:%d\n", dir, key, OUTway);
    
    return (int)LockSet.size()-1;
}
/////////////////////////////////////////////////////////////

int doublTest(int iXr, int iMut)
{
    for ( int n=0; n<vMinSet.size(); n++ )  {
        if ( vMinSet[n].first==iXr && vMinSet[n].second==iMut )
            return 1;
    }
    return 0;
}
/////////////////////////////////////////////////////////////

int SAMPLE:: proSampMut(  )
{
        XROSOMA *pXro;
        vector <TOUCH_SET>::iterator ItRoot;
//        vector <MUTANT>::iterator ItMutSet;
        int state, Lway, nSegm;
    
    pCurSmpl = this;
    memset(ntabs, '\t', sizeof(ntabs)-1);        // testTrace
    
    for ( int iCurXr=0; iCurXr < NO_HUMAN_XRO; iCurXr++ ) {
        pXro = &vecDNK[iCurXr];
        
        fprintf(tsld, "\n*****Sample #%d of X.%d:\n", (int)(this - &vecSAMPL[0]), iCurXr );
        
        if ( pXro->vTouchTree.size()==0 )   // testCase  for limited noXRO
            continue;
        
        for ( int iCurMu=0; iCurMu < vMutAPO[iCurXr].size(); iCurMu++ ) {
            ItPrimeMut = vMutAPO[iCurXr].begin() + iCurMu;
            int segCurMut = (int)ItPrimeMut->nucPos / SEGM_SIZE;// vMutAPO[iCurXr][iCurMu].nucPos/SEGM_SIZE;
            vMinSet.clear();
            LockSet.clear();
            MinLway     = (int)(pXro->Xsize / SEGM_SIZE);
            ItRoot = lower_bound(pXro->vTouchTree.begin(), pXro->vTouchTree.end(),
                                TOUCH_SET(segCurMut), lesser_TOUCH);
            
            fprintf(tsld, "\n==M.%d[%d] Root[%d]", iCurMu, segCurMut,
                    (ItRoot==pXro->vTouchTree.end()) ? -1 : ItRoot->ownTAR );
                
            if ( cnt_ALL_Mut(iCurXr, iCurMu, ItRoot) > 0 ) {  //в сегменте iCurMu есть еще мутации
                            // и если CurMu==Root, проверяем сегменты всех Fork на мутации в них
                fprintf(tsld, "\t\tMultyMut={%ld} minL=0\n", vMinSet.size() );

                MinLway = 0;
                printMutLines( iCurXr, iCurMu, "SameSeg" );
                continue;
            };
//-------------------
//---- UP way
            state = 0;
            if ( iCurMu+1 >= vMutAPO[iCurXr].size() )
                state = ( ItRoot == pXro->vTouchTree.end() ) ? 1 : 2;
            else    {
                if ( ItRoot == pXro->vTouchTree.end() )
                    state = 3;
                else
                    state = ( vMutAPO[iCurXr][iCurMu+1].nucPos / SEGM_SIZE > ItRoot->ownTAR ) ? 2 : 3;
            }
            
            fprintf(tsld, "\n  UPstate=%d\n", state);
            
            switch (state) {//                cM=last  R=end()
                case 1:     // ---------------#--------X----------->
                            // последняя мутация и выше нее пересечений нет
                    break;
                    
                case 2:     // a) следеющая мутация ВЫШЕ, чем ItRoot->ownTAR
                            //       cM    R      cM+1
                            //  -----#-----X------#--------------->
                            // б) последняя мутация, но выше есть пересечения
                            //        cM==last    R
                            //    -----#----------X-----X--X------>
                    
 //                   LockSet.push_back(ItRoot->ownTAR*SEGM_SIZE+iCurXr);
//                    LockSet.push_back(-(ItRoot->ownTAR*SEGM_SIZE+iCurXr));  //оттуда пришли
                    pXro->addLock(iCurXr, ItRoot->ownTAR, 1, ItRoot->ownTAR - segCurMut);
                    pXro->TreeUP( 0, ItRoot, iCurMu+1, ItRoot->ownTAR - segCurMut );
                    break;
                    
                case 3:     // а) следеющая мутация есть, но пересечений нет
                            //   R==end()  cM     cM+1
                            //  --X--------#------#--------------->
                            // б) следеющая мутация НЕ ВЫШЕ, чем Root
                            //    cM     cM+1    R
                            //  --#--------#-----X--------------->
                            //  --#--------X--------------------->
                    nSegm = (int)vMutAPO[iCurXr][iCurMu+1].nucPos / SEGM_SIZE;
                    Lway = nSegm - segCurMut;
                    
                    fprintf(tsld, "    nearMut.%d[%d] L=%d:$%d minSet{%ld}->",
                            iCurMu+1, nSegm, Lway, MinLway, vMinSet.size() );
                    
                    if ( Lway <= MinLway ) {
                        if ( Lway < MinLway)   {
                            vMinSet.clear();
                            MinLway = Lway;
                        }
                        vMinSet.push_back(pair <int, int>(iCurXr, iCurMu+1));
                        pCurSmpl->cnt_ALL_Mut(iCurXr, iCurMu+1, ItRoot);
                    }
                    
                    fprintf(tsld, " {%ld} $%d\n", vMinSet.size(), MinLway);
                    
                    break;
                    
                default:
                    printf("BAD UPstate: iXr=%d iCurMu.%d[%ld] Root=%d\n", iCurXr, iCurMu,
                           vMutAPO[iCurXr][iCurMu].nucPos,
                           (ItRoot == pXro->vTouchTree.end()) ? -1 : ItRoot->ownTAR );
                    continue;
            }

////---- DOWN way
            state = 0;
            if ( iCurMu == 0 )
                state = ( ItRoot == pXro->vTouchTree.begin() ) ? 1 : 2;
            else    {
                if ( ItRoot-1 < pXro->vTouchTree.begin() )
                    state = 3;
                else
                    state = ( (ItRoot-1)->ownTAR > vMutAPO[iCurXr][iCurMu-1].nucPos/SEGM_SIZE ) ? 2 : 3;
            }
            
            fprintf(tsld, "  DOWNstate=%d\n", state);
            
            switch (state) {
                case 1:     // первая мутация и ниже нее пересечений нет
                            //      cM=0     R=begln()
                    break;  // <-----#-------X----------------
                    
                case 2:     // а) предыдущая мутация НИЖЕ, чем предыдущее пересечение
                            //      M-1      R-1    cM  R
                            // <----#----(x)--X-----#---X-----------
                            // б) первая мутация, но ниже есть пересечения
                            //                R-1   cM=0  R
                            // <---------(x)--X-----#-----X-----
//                    LockSet.push_back( -( (ItRoot-1)->ownTAR*SEGM_SIZE+iCurXr) );
//                    LockSet.push_back(    (ItRoot-1)->ownTAR*SEGM_SIZE+iCurXr  );  //оттуда пришли
                    if ( pXro->addLock(iCurXr, (ItRoot-1)->ownTAR, -1, segCurMut-(ItRoot-1)->ownTAR) >=0 )
                        pXro->TreeDOWN( 0, (ItRoot-1), iCurMu-1, segCurMut-(ItRoot-1)->ownTAR );
                    break;
                    
                case 3:     // а) предыдущая мутация есть, но пересечений нет
                            //           M-1    cM    R=begln()
                            // <----------#-----#-----X---------------
                            // б) предыдущая мутация НЕ НИЖЕ, чем предыдущее пересечение
                            //        R-1  M-1   cM    R
                            //  <-----X----#-----#-----X---------
                            //  <----------#X----#-----X---------
                    nSegm = (int)vMutAPO[iCurXr][iCurMu-1].nucPos / SEGM_SIZE;
                    Lway = segCurMut - nSegm;
                    
                    fprintf(tsld, "    nearMut.%d[%d] L=%d:$%d minSet{%ld}",
                            iCurMu-1, nSegm, Lway, MinLway, vMinSet.size() );
                    
                    if ( Lway <= MinLway ) {
                        if ( Lway < MinLway )   {
                            vMinSet.clear();
                            MinLway = Lway;
                        }
//  может быть дубликат (Lway==MinLway), если по ветке 'UP way' эту мутацию уже определили
//  дубликат откидывается при печати
                        if ( ! doublTest(iCurXr, iCurMu-1) )    {
                            vMinSet.push_back(pair <int, int>(iCurXr, iCurMu-1));
                            pCurSmpl->cnt_ALL_Mut(iCurXr, iCurMu-1, ItRoot-1);
                        }
                        else
                            fprintf(tsld, " DOUBLE ");
                    }
                    
                    fprintf(tsld, " ->{%ld} $%d\n", vMinSet.size(), MinLway);
                    
                    break;
                    
                default:
                    printf("BAD DOWNstate: iXr=%d iCurMu.%d[%ld] Root=%d\n", iCurXr, iCurMu,
                           vMutAPO[iCurXr][iCurMu].nucPos,
                           (ItRoot == pXro->vTouchTree.end()) ? -1 : ItRoot->ownTAR );
                    continue;
            }
            
            fprintf(tsld, "--$%d {%ld}:", MinLway, vMinSet.size());
            for ( int n=0; n<vMinSet.size(); n++ )
                fprintf(tsld, "  (X%d %d)", vMinSet[n].first, vMinSet[n].second);
            fprintf(tsld, "\n");
            
            printMutLines( iCurXr, iCurMu, " " );
        }
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////
                
int XROSOMA:: TreeUP( int Rlevel, vector <TOUCH_SET>::iterator ItRoot, int iMut, int OUTway )
//  Работает в пределах поддерева от ItRoot до мутации iMut (до конца, если мутаций нет)
//           пока длина пути <= уже найденной ранее минимальной (MinLway)
//  Н.У. : 1) ItRoot < end();  2) ItRoot < iMut; 3) iMut <= size() м.б. больше последней
//         4) в Fork_ах ItRoot мутаций нет
//  OUTway длина пути от начальной (iCurMu) мутации до данного корня (ItRoot)
// a) следеющая мутация ВЫШЕ, чем ItRoot->ownTAR
//       cM        OUTway      R                 iM
//  -----#---------------------X------(X)--(X)---#------>
// б) последняя мутация, но выше м.б. пересечения
//    -----#-------------------X------(X)---(X)--------->
//      cM==last    OUTwayR    R
{
    
    int iRootXr = (int)(this - &vecDNK[0]);
    int Lway;

    
//    int segNearMut =  1+ (int)(pCurSmpl->vMutAPO[iRootXr].end( )-1)->nucPos / SEGM_SIZE; // max
    int segNearMut = 1+ (vTouchTree.end()-1)->ownTAR; // max
    if ( iMut >= (int)pCurSmpl->vMutAPO[iRootXr].size() )
        iMut = -1;  // no more muts
    
    fprintf(tsld, "%.*sTreeUP(%d [%d] iM.%d oW%d): %s",
            Rlevel+1, ntabs,  Rlevel, ItRoot->ownTAR, iMut, OUTway,
            ((iMut<0)? "  NO more UP_Muts\n" : " " ));
    
    if ( iMut >= 0 ) {      // a) следеющая мутация ВЫШЕ, чем ItRoot->ownTAR
        segNearMut = (int) (pCurSmpl->vMutAPO[iRootXr][iMut].nucPos / SEGM_SIZE);
        Lway = OUTway + (segNearMut - ItRoot->ownTAR);
        
        fprintf(tsld, "  iM.%d[%d] L%d::$%d MinSet{%ld} %s", iMut, segNearMut,
                Lway, MinLway, vMinSet.size(), ((Lway > MinLway) ? " skip(L>$)\n" : " ") );
        
        if ( Lway <= MinLway )  {
            if ( Lway < MinLway )   {
                fprintf(tsld, " Clear(L<$) ");
                vMinSet.clear();
                MinLway = Lway;
            }
            vMinSet.push_back(pair <int, int>(iRootXr, iMut));
            pCurSmpl->cnt_ALL_Mut(iRootXr, iMut, ItRoot);
            fprintf(tsld, " Add $%d:{%ld}\n", MinLway, vMinSet.size());
            //      Не покидаем, т.к. от ItRoot до iMut м.б. пересечения: -X--(X)--(X)-#->
        }
        
    }
//    else {  // выше ItRoot мутаций нет !!!!!
//        fprintf(tsld, "  NO more Muts\n");
//    }
    
//  проход по всем разветвлениям до сегмента найденной ближайшей мутации
//  и поиск более близкой мутации в ветках: -X--(X)--(X)-#------->

//    fprintf(tsld, "%.*s  LockSet{%ld} =", Rlevel+1, ntabs, LockSet.size() );
//    for ( int n=0; n<LockSet.size(); n++ ) fprintf(tsld," %d:%d", LockSet[n].first,LockSet[n].second);
//    fprintf(tsld, "\n");
    
    vector <TOUCH_SET>::iterator ItNextTT, ItForkTT;
    vector <MUTANT >::iterator ItForkMut;
    vector<pair<int, int>> vMinS;
    for ( ItNextTT=ItRoot; ItNextTT != vTouchTree.end(); ItNextTT++ )    {
        Lway = OUTway + (ItNextTT->ownTAR - ItRoot->ownTAR);
        
        fprintf(tsld, "%.*s  nexTT[%d] L=%d:$%d Fork{%ld}=", Rlevel+1, ntabs,
                    ItNextTT->ownTAR, Lway, MinLway, ItNextTT->vFork.size());
        for ( int n=0; n<ItNextTT->vFork.size(); n++)
            fprintf(tsld, " {X%d %d}", ItNextTT->vFork[n].first, ItNextTT->vFork[n].second);
        fprintf(tsld, "\n");
        
        if ( ItNextTT->ownTAR >= segNearMut || Lway > MinLway ) {
            if ( ItNextTT->ownTAR >= segNearMut )
                fprintf(tsld, "%.*s  STOP_UPway: nxTT[%d]>=Mut[%d]\n",
                        Rlevel+1, ntabs, ItNextTT->ownTAR, segNearMut);
            else
                fprintf(tsld, "%.*s  STOP_UPway: Lway=%d > MinLway=%d\n",
                        Rlevel+1, ntabs, Lway, MinLway);
            break;
        }
// -----тест на наличие мутаций в Fork_ах точки ветвления  NextTT
        if ( ItNextTT > ItRoot )    {
            if ( pCurSmpl->cntMutInRoot(ItNextTT, vMinS) > 0 )  {
                if ( Lway < MinLway )   {
                    vMinSet.clear();
                    MinLway = Lway;
                }
                vMinSet.insert(vMinSet.end(), vMinS.begin(), vMinS.end());
                fprintf(tsld, "%.*s  STOP_UPway: MultiMut inRoot{%ld}\n",
                        Rlevel+1, ntabs, vMinS.size());
                break;
            }
        }
        for ( int nf=0; nf<ItNextTT->vFork.size(); nf++ )    {    // по веткам
            int iFrkXr = ItNextTT->vFork[nf].first;
            int FrkSeg = ItNextTT->vFork[nf].second;
            int iFrMut, nL;
            ItForkTT = lower_bound(vecDNK[iFrkXr].vTouchTree.begin(), vecDNK[iFrkXr].vTouchTree.end(),
                                  TOUCH_SET(FrkSeg), lesser_TOUCH);
            ItForkMut = lower_bound(pCurSmpl->vMutAPO[iRootXr].begin(),
                                    pCurSmpl->vMutAPO[iRootXr].end(),
                                    MUTANT(FrkSeg*SEGM_SIZE), lesser_MUT);
            fprintf(tsld, "%.*s   Fork_UP.%d {X%d [%d]} :", Rlevel+1, ntabs,
                        nf, iFrkXr, FrkSeg);
            if ( ItForkMut==ItPrimeMut ) {
                fprintf(tsld, "\tLOOP for Upway - goto DOWNway !!!!!\n");
                goto DOWNway;//continue;               //  LOOP !!!!!!
            }
// -----UPway
            if (  (nL=addLock(iFrkXr, FrkSeg, 1, Lway)) >= 0 )   {
//                fprintf(tsld, " add_UP_LOCK=%d:%d\n", LockSet[nL].first,LockSet[nL].second);
                iFrMut = (ItForkMut==pCurSmpl->vMutAPO[iRootXr].end())
                        ? (int)pCurSmpl->vMutAPO[iRootXr].size()
                        : (int)(ItForkMut-pCurSmpl->vMutAPO[iRootXr].begin());
                vecDNK[iFrkXr].TreeUP( Rlevel+1, ItForkTT, iFrMut, Lway );
            }  else {
//                fprintf(tsld, " UP_LOCKed=%d:%d\n",LockSet[-nL-1].first,LockSet[-nL-1].second);
                }
// -----DOWNway
DOWNway:
            ItForkMut -= 1;
            fprintf(tsld, "%.*s   Fork_DOWN.%d {X%d [%d]} :", Rlevel+1, ntabs,
                    nf, iFrkXr, FrkSeg);
            if ( ItForkMut==ItPrimeMut )    {
                fprintf(tsld, "\tLOOP for DOWNway - goto next Fork !!!!!\n");
                continue;               //  LOOP !!!!!!
            }
            if ( ( nL = addLock(iFrkXr, FrkSeg, -1, Lway) ) >= 0 )   {
//                fprintf(tsld, " add_DOWN_LOCK=%d:%d\n",LockSet[nL].first,LockSet[nL].second);
                iFrMut = (ItForkMut < pCurSmpl->vMutAPO[iRootXr].begin())  ? -1
                                : (int)(ItForkMut-pCurSmpl->vMutAPO[iRootXr].begin());
                vecDNK[iFrkXr].TreeDOWN( Rlevel+1, ItForkTT, iFrMut, Lway );
            } else  {
//               fprintf(tsld, " DOWN_LOCKed=%d:%d\n",LockSet[-nL-1].first,LockSet[-nL-1].second);
            }
        }
    }

    return 0;
}
/////////////////////////////////////////////////////////////

int SAMPLE:: getNearMutUP(int iXro, vector <TOUCH_SET>::iterator ItRoot) //, int iMut) //, int OUTway)
{
//  поиск ближайшей мутации ВЫШЕ базового сегмента
//  Returns: индекс мутации;
//           -1 если ВЫШЕ базового сегмента мутаций нет
     
//     int segRoot= ( ItRoot == vecDNK[iXro].vTouchTree.end() )
//                ?  (vecDNK[iXro].vTouchTree.end()-1)->ownTAR : ItRoot->ownTAR;
     
     vector <MUTANT >::iterator ItNearMut = lower_bound( vMutAPO[iXro].begin( ), vMutAPO[iXro].end( ),
                                                        MUTANT(ItRoot->ownTAR*SEGM_SIZE),  lesser_MUT);
     if ( ItNearMut == vMutAPO[iXro].end( ) )
         return -1;

     return (int)( ItNearMut - vMutAPO[iXro].begin( ) );
}
/////////////////////////////////////////////////////////////

int XROSOMA:: TreeDOWN( int Rlevel, vector <TOUCH_SET>::iterator ItRoot, int iMut, int OUTway )
{
//  Н.У. : 1) ItRoot >= begin();  2) ItRoot > iMut; 3) iMut <= 0 м.б. меньше первой
//         4) в Fork_ах ItRoot мутаций нет
//  OUTway длина пути от начальной (iCurMu) мутации до данного корня (ItRoot)
// а) мутация НИЖЕ, чем пересечение
//      iM        R    OUTway       cM
// <----#----(X)--X-----------------[#]--------
// б) первая мутация, но ниже есть пересечения
//   iM=-1        R    OUTway     cM=0
// <---------(X)--X-----..........[#]-----------

    int iRootXr = (int)(this - &vecDNK[0]);
    int Lway;
    int segNearMut = - 1;
    
    fprintf(tsld, "%.*sTreeDOWN(%d [%d] iM.%d oW%d): %s",
            Rlevel+1, ntabs,  Rlevel, ItRoot->ownTAR, iMut, OUTway,
            ((iMut<0)? "  NO more DOWN_Muts\n" : " " ));
    
    if ( iMut >= 0 ) {
        segNearMut = (int) (pCurSmpl->vMutAPO[iRootXr][iMut].nucPos / SEGM_SIZE);
        Lway = OUTway + (ItRoot->ownTAR - segNearMut);
        
        fprintf(tsld, "  iM.%d[%d] L%d::$%d MinSet{%ld} %s", iMut, segNearMut,
                Lway, MinLway, vMinSet.size(), ((Lway > MinLway) ? " skip(L>$)\n" : " ") );
        
        if ( Lway <= MinLway )  {
            if ( Lway < MinLway )   {
                fprintf(tsld, " Clear(L<$) ");
                vMinSet.clear();
                MinLway = Lway;
            }
            vMinSet.push_back(pair <int, int>(iRootXr, iMut));
            pCurSmpl->cnt_ALL_Mut(iRootXr, iMut, ItRoot);
            fprintf(tsld, " Add $%d:{%ld}\n", MinLway, vMinSet.size());
            //      Не покидаем, т.к. от ItRoot до iMut м.б. пересечения:  <--#--(X)---X--
        }
        
    }
//    else {  // ниже ItRoot мутаций нет !!!!!
//            fprintf(tsld, "  NO more Muts\n");
//        }

//  проход по всем разветвлениям до сегмента найденной ближайшей мутации
    //  и поиск более близкой мутации в ветках: <--#--(X)---X--
    
//    fprintf(tsld, "%.*s  LockSet{%ld} =", Rlevel+1, ntabs, LockSet.size() );
//    for ( int n=0; n<LockSet.size(); n++ ) fprintf(tsld," %d:%d",LockSet[n].first,LockSet[n].second);
//    fprintf(tsld, "\n");
    
    vector <TOUCH_SET>::iterator ItNextTT, ItForkTT;
    vector <MUTANT >::iterator ItForkMut;
    vector<pair<int, int>> vMinS;
    for ( ItNextTT=ItRoot; ItNextTT >= vTouchTree.begin(); ItNextTT-- )    {
        Lway = OUTway + (ItRoot->ownTAR - ItNextTT->ownTAR);
        
        fprintf(tsld, "%.*s  nexTT[%d] Fork{%ld}=", Rlevel+1, ntabs,
                ItNextTT->ownTAR, ItNextTT->vFork.size());
        for ( int n=0; n<ItNextTT->vFork.size(); n++)
            fprintf(tsld, " {X%d %d}", ItNextTT->vFork[n].first, ItNextTT->vFork[n].second);
        fprintf(tsld, "\n");
        
        if ( ItNextTT->ownTAR <= segNearMut || Lway > MinLway ) {
            if ( ItNextTT->ownTAR <= segNearMut )
                fprintf(tsld, "%.*s  STOP_DOWNway: nxTT[%d]<=Mut[%d]\n",
                        Rlevel+1, ntabs, ItNextTT->ownTAR, segNearMut);
            else
                fprintf(tsld, "%.*s  STOP_DOWNway: Lway=%d > MinLway=%d\n",
                        Rlevel+1, ntabs, Lway, MinLway);
            break;
        }
// -----тест на наличие мутаций в Fork_ах точки ветвления  NextTT
        if ( ItNextTT < ItRoot )    {
            if ( pCurSmpl->cntMutInRoot(ItNextTT, vMinS) > 0 )  {
                if ( Lway < MinLway )   {
                    vMinSet.clear();
                    MinLway = Lway;
                }
                vMinSet.insert(vMinSet.end(), vMinS.begin(), vMinS.end());
                fprintf(tsld, "%.*s  STOP_DOWNway: MultiMut inRoot{%ld}\n",
                        Rlevel+1, ntabs, vMinS.size());
                break;
            }
        }

        for ( int nf=0; nf<ItNextTT->vFork.size(); nf++ )    {    // по веткам
            int iFrkXr = ItNextTT->vFork[nf].first;
            int FrkSeg = ItNextTT->vFork[nf].second;
            int iFrMut, nL;
            ItForkTT = lower_bound(vecDNK[iFrkXr].vTouchTree.begin(), vecDNK[iFrkXr].vTouchTree.end(),
                                   TOUCH_SET(FrkSeg), lesser_TOUCH);
            ItForkMut = lower_bound(pCurSmpl->vMutAPO[iRootXr].begin(),
                                    pCurSmpl->vMutAPO[iRootXr].end(),
                                    MUTANT(FrkSeg*SEGM_SIZE), lesser_MUT);
            fprintf(tsld, "%.*s   Fork_UP.%d {X%d [%d]} :", Rlevel+1, ntabs,
                    nf, iFrkXr, FrkSeg);
            if ( ItForkMut==ItPrimeMut )    {
                fprintf(tsld, "\tLOOP for UPway - goto DOWNway !!!!!\n");
                goto DOWNway;   //continue;               //  LOOP !!!!!!
            }
// -----UPway
            if ( (nL=addLock(iFrkXr, FrkSeg, 1, Lway)) >= 0 )   {
//                fprintf(tsld, " add_UP_LOCK=%d:%d\n",LockSet[nL].first,LockSet[nL].second);
                iFrMut = (ItForkMut==pCurSmpl->vMutAPO[iRootXr].end())
                        ? (int)pCurSmpl->vMutAPO[iRootXr].size()
                        : (int)(ItForkMut-pCurSmpl->vMutAPO[iRootXr].begin());
                vecDNK[iFrkXr].TreeUP( Rlevel+1, ItForkTT, iFrMut, Lway );
            }  else {
//                fprintf(tsld, " UP_LOCKed=%d:%d\n",LockSet[-nL-1].first,LockSet[-nL-1].second);
            }
// -----DOWNway
DOWNway:
            ItForkMut -= 1;
            fprintf(tsld, "%.*s   Fork_DOWN.%d {X%d [%d]} :", Rlevel+1, ntabs,
                    nf, iFrkXr, FrkSeg);
            if ( ItForkMut==ItPrimeMut )    {
                fprintf(tsld, "\tLOOP for DOWNway - goto next Fork !!!!!\n");
                continue;               //  LOOP !!!!!!
            }
            if ( (nL=addLock(iFrkXr, FrkSeg, -1, Lway) ) >= 0 )  {
//                fprintf(tsld, " add_DOWN_LOCK=%d:%d\n", LockSet[nL].first,LockSet[nL].second);
                iFrMut = (ItForkMut < pCurSmpl->vMutAPO[iRootXr].begin())  ? -1
                                : (int)(ItForkMut - pCurSmpl->vMutAPO[iRootXr].begin());
                vecDNK[iFrkXr].TreeDOWN( Rlevel+1, ItForkTT, iFrMut, Lway );
            }  else {
//                fprintf(tsld, " DOWN_LOCKed=%d:%d\n", LockSet[-nL-1].first,LockSet[-nL-1].second);
            }
        }
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////

int SAMPLE:: getNearMutDOWN(int iXro, vector <TOUCH_SET>::iterator ItRoot)
{
//  поиск ближайшей мутации НИЖЕ базового сегмента
    
    int segRoot = ItRoot->ownTAR;
    vector <MUTANT >::iterator ItNearMut= lower_bound( vMutAPO[iXro].begin( ), vMutAPO[iXro].end( ),
                                                      MUTANT(segRoot*SEGM_SIZE),  lesser_MUT);
    if ( ItNearMut==vMutAPO[iXro].end( ) )
        ItNearMut--;
    
    while ( ItNearMut >= vMutAPO[iXro].begin( ) ) {
        if ( ItNearMut->nucPos/SEGM_SIZE < segRoot )
            break;
        ItNearMut--;
    }
    if ( ItNearMut < vMutAPO[iXro].begin( ) )
        return -1;
    
    return (int)( ItNearMut - vMutAPO[iXro].begin( ) );
    
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

int SAMPLE:: cnt_ALL_Mut(int ixXr, int ixMu, vector <TOUCH_SET>::iterator ItRoot)
//  1. Для мутации ixMu собирает все мутации, которые лежат в этом же сегменте
//          (сама мутация ixMu в этот перечень не входит)
// 2. Если мутация ixMu, лежит в сегменте пересечения хромосом (ItNearRoot),
//      то для всех сегментов пересечения (vFork) собираются мутации, попадающие в них.
//  Returns: количество найденных мутаций,
//           которые дополняются в глобальный vMinSet
{
// 1.
    int RetC = cntMutInSeg(ixXr, ixMu);
//----------
    if (    ItRoot == vecDNK[ixXr].vTouchTree.end() ||
            ItRoot <  vecDNK[ixXr].vTouchTree.begin() )
        return RetC;
    
    if ( ItRoot->ownTAR != (int)(vMutAPO[ixXr][ixMu].nucPos / SEGM_SIZE) )
        return RetC;
// 2.
//----------
    RetC += cntMutInRoot(ItRoot);
/*
    for ( int nFr=0; nFr<ItRoot->vFork.size(); nFr++ )  {
        int iXroFork  = ItRoot->vFork[nFr].first;
        int segmFork  = ItRoot->vFork[nFr].second;
        vector <MUTANT>::iterator ItForkMut =                       // first Mut near segmFork
                lower_bound(vMutAPO[iXroFork].begin( ), vMutAPO[iXroFork].end( ),
                            MUTANT(segmFork*SEGM_SIZE),  lesser_MUT);
        if ( ItForkMut == vMutAPO[iXroFork].end( ) || ItForkMut->nucPos/SEGM_SIZE != segmFork )
            continue;                               // next Fork
        
        for ( ItForkMut=ItForkMut; ItForkMut != vMutAPO[iXroFork].end( ); ItForkMut++ )    {
            if ( (int)(ItForkMut->nucPos / SEGM_SIZE) != segmFork )
                break;
            int nMu = (int)( ItForkMut - vMutAPO[iXroFork].begin() );
            vMinSet.push_back(pair <int, int>(iXroFork, nMu));
            RetC++;
        }
        
    }
*/
    return RetC;
}
/////////////////////////////////////////////////////////////

int SAMPLE:: cntMutInSeg(int ixXr, int ixMu)
{
//  Для мутации ixMu считает все мутации, которые лежат в этом же сегменте
//      сама мутация ixMu в этот перечень не входит
//  Returns: количество найденных мутаций, которые дополняются в глобальный vMinSet
    
    int nMu;
    int cnt =0;
    int segMut  = (int)(vMutAPO[ixXr][ixMu].nucPos / SEGM_SIZE);
    
    for ( nMu=ixMu+1; nMu < vMutAPO[ixXr].size();  nMu++)   {
        if ( (int)(vMutAPO[ixXr][nMu].nucPos / SEGM_SIZE) != segMut )
            break;
        vMinSet.push_back(pair <int, int>(ixXr, nMu));
        cnt++;
    }
    for ( nMu=ixMu-1; nMu >= 0;  nMu--)   {
        if ( (int)(vMutAPO[ixXr][nMu].nucPos / SEGM_SIZE)!=segMut )
            break;
        vMinSet.push_back(pair <int, int>(ixXr, nMu));
        cnt++;
    }
    
    return cnt;
}
/////////////////////////////////////////////////////////////

int SAMPLE:: cntMutInRoot(vector <TOUCH_SET>::iterator ItRoot)
{
//  Для сегмента пересечения хромосом (ItRoot) подсчитываются мутации,
//      попавшие в сегменты пересечения (vFork).
//  Returns: количество найденных мутаций, которые дополняются в глобальный vMinSet
//
    int cnt =0;
    int iXroFork, segmFork, iMut;
    vector <MUTANT>::iterator ItForkMut;
    
    for ( int nFr=0; nFr<ItRoot->vFork.size(); nFr++ )  {
        iXroFork  = ItRoot->vFork[nFr].first;
        segmFork  = ItRoot->vFork[nFr].second;
        ItForkMut =  lower_bound(vMutAPO[iXroFork].begin( ), vMutAPO[iXroFork].end( ),
                                    MUTANT(segmFork*SEGM_SIZE),  lesser_MUT);   // first Mut
        if ( ItForkMut == vMutAPO[iXroFork].end( ) ||
            ItForkMut->nucPos/SEGM_SIZE != segmFork )
            continue;
        
        iMut = (int)(ItForkMut - vMutAPO[iXroFork].begin( ));
        vMinSet.push_back(pair <int, int>(iXroFork, iMut));
        cnt++;
//------test 'segmFork' to multyMUT
        for ( iMut=iMut+1; iMut<vMutAPO[iXroFork].size(); iMut++ )  {
            if ( (vMutAPO[iXroFork][iMut].nucPos / SEGM_SIZE) != segmFork )
                break;
            vMinSet.push_back(pair <int, int>(iXroFork, iMut));
            cnt++;
        }
        
    }
    return cnt;
}
/////////////////////////////////////////////////////////////

int SAMPLE:: cntMutInRoot(vector <TOUCH_SET>::iterator ItRoot, vector<pair<int, int>> vMinS)
{
    //  Для сегмента пересечения хромосом (ItRoot) подсчитываются мутации,
    //      попавшие в сегменты пересечения (vFork).
    //  Returns: количество найденных мутаций, которые дополняются в глобальный vMinSet
    //
    int cnt =0;
    int iXroFork, segmFork, iMut;
    vector <MUTANT>::iterator ItForkMut;
    
    vMinS.clear();
    for ( int nFr=0; nFr<ItRoot->vFork.size(); nFr++ )  {
        iXroFork  = ItRoot->vFork[nFr].first;
        segmFork  = ItRoot->vFork[nFr].second;
        ItForkMut =  lower_bound(vMutAPO[iXroFork].begin( ), vMutAPO[iXroFork].end( ),
                                 MUTANT(segmFork*SEGM_SIZE),  lesser_MUT);   // first Mut
        if ( ItForkMut == vMutAPO[iXroFork].end( ) ||
            ItForkMut->nucPos/SEGM_SIZE != segmFork )
            continue;
        
        iMut = (int)(ItForkMut - vMutAPO[iXroFork].begin( ));
        vMinS.push_back(pair <int, int>(iXroFork, iMut));
        cnt++;
        //------test 'segmFork' to multyMUT
        for ( iMut=iMut+1; iMut<vMutAPO[iXroFork].size(); iMut++ )  {
            if ( (vMutAPO[iXroFork][iMut].nucPos / SEGM_SIZE) != segmFork )
                break;
            vMinS.push_back(pair <int, int>(iXroFork, iMut));
            cnt++;
        }
        
    }
    return cnt;
}
/////////////////////////////////////////////////////////////

void printTableCap( )
{
    fprintf(fRezult, "Sampl\t#S\tXro_1\t#M1\tMutSeg_1\tMutPos_1\tXro_2\t#M2\tMutSeg_2\tMutPos_2\t" );
    fprintf(fRezult, "shortWay\tmREF\tmALT\t#CrossX1\tCrossXx\n" );
}
////////////////////////////////////////////////////////

void SAMPLE:: printMutLines( int iXro, int iCurMu, const char *comment  )
{
    if ( iCurMu==0 )
        fprintf(fRezult, "%s\t%ld\t", this->SampName.c_str(), this-&vecSAMPL[0] );
    else
        fprintf(fRezult, "\t%ld\t", this-&vecSAMPL[0] );
    
    int segCurMut = (int)vMutAPO[iXro][iCurMu].nucPos / SEGM_SIZE;
    for ( int n=0; n<vMinSet.size(); n++ )  {
        if ( n>0 )
            fprintf(fRezult, "\t%ld\t", this-&vecSAMPL[0] );
        int iNearXro = vMinSet[n].first;
        int iNearMut = vMinSet[n].second;
        int segNearMut = (int)vMutAPO[iNearXro][iNearMut].nucPos / SEGM_SIZE;
        fprintf(fRezult, "%d\t%d\t%d\t%ld\t%d\t%d\t%d\t%ld\t%d\t%c\t%c\t%d\t%d\t%s\n",
            iXro, iCurMu, segCurMut, vMutAPO[iXro][iCurMu].nucPos,
            iNearXro, iNearMut, segNearMut, vMutAPO[iNearXro][iNearMut].nucPos,
            MinLway, vMutAPO[iXro][iCurMu].nucREF, vMutAPO[iXro][iCurMu].nucALT,
            0, 0, comment);   // cntForks, cntForksOWN, comment
    }
    if ( vMinSet.size()==0 )
        fprintf(fRezult, "\n");
    return;
}
/////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

void tst_L_B()
{
    XROSOMA *pXro;
    vector < TOUCH_SET >::iterator itTT;
    
    pXro = &vecDNK[0];
    for ( int n=0; n<=250; n++ )  {
        itTT = lower_bound(pXro->vTouchTree.begin(), pXro->vTouchTree.end(),
                           TOUCH_SET(n), lesser_TOUCH);
        printf("Finding=%d...\t", n);
        if ( itTT==pXro->vTouchTree.end() )
            printf("NOT FOUND==end()\n");
        else
            printf("Finded=%d ix=%d\n", itTT->ownTAR, (int)(itTT-pXro->vTouchTree.begin()));
        continue;
    }
    
    return;
}
/////////////////////////////////////////////////////////////

int SAMPLE:: tstSampMut(  )
{
    XROSOMA *pXro;
    vector <TOUCH_SET>::iterator ItRoot;
    vector <MUTANT>::iterator ItMutSet;
    
    pCurSmpl = this;
    
    for ( int iCurXr=0; iCurXr <= 0; iCurXr++ ) {
        pXro = &vecDNK[iCurXr];
        ItMutSet = vMutAPO[iCurXr].begin( );
        if ( pXro->vTouchTree.size()==0 )
            continue;
        //        fprintf(tsld, "\nXRO.%d =====================\n", ItSMP->iCurXr);
        for ( int iCurMu=0; iCurMu < vMutAPO[iCurXr].size(); iCurMu++ ) {
            int segCurMu = (int)vMutAPO[iCurXr][iCurMu].nucPos / SEGM_SIZE;
            MinLway     = (int)(pXro->Xsize / SEGM_SIZE);
            vMinSet.clear();
            ItRoot = lower_bound(pXro->vTouchTree.begin(), pXro->vTouchTree.end(),
                                 TOUCH_SET(segCurMu), lesser_TOUCH);
            printf("Mut.%d = %d [%d] ...\t", iCurMu, (int)vMutAPO[iCurXr][iCurMu].nucPos, segCurMu);
            if ( ItRoot==pXro->vTouchTree.end() )
                printf("NOT FOUND==end()\n");
            else
                printf("TT.%d TAR=%d\n", (int)(ItRoot-pXro->vTouchTree.begin()),
                       ItRoot->ownTAR );
            continue;
        }
    }
    return 0;
}
//////////////////////////////////////////////////////////////////////



