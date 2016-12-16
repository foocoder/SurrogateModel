// ---- Program Info Start----
//FileName:     MOPM.cpp
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-01-03 10:40:45
// ---- Program Info End  ----

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "mopm.h"

using namespace std;

#define EPSINON 1e-18

template <typename T>
vector<size_t> fnSortIndex(const vector<T> & v, int iType)
{
    //iType : 0 for ascend and 1 for descend
    vector<size_t> idx(v.size());
    iota(idx.begin(),idx.end(),0);
    if(iType == 0)
    {
        sort(idx.begin(),idx.end(),[&v](size_t i1, size_t i2){ return v[i1] < v[i2]; });
    }
    else if(iType == 1)
    {
        sort(idx.begin(),idx.end(),[&v](size_t i1, size_t i2){ return v[i1] > v[i2]; });
    }
    return idx;
}

NSGAII::NSGAII(){}
NSGAII::~NSGAII(){}

NSGAII::NSGAII
(
 int iPopSize,
 int iPopDims,
 int iDataSize,
 int iDataDims,
 int iObjDims,
 const vector<myBitSet<M>> & vbitDatabase
 ):
    _vbitDatabase(vbitDatabase)
{
    _iPopSize  = iPopSize;
    _iPopDims  = iPopDims;
    _iDataSize = iDataSize;
    _iDataDims = iDataDims;
    _iObjDims  = iObjDims;

    _vbitItemDatabase.reserve(iDataDims);
    _viItemLength.reserve(iDataDims);
    _viTransLength.reserve(iDataSize);
    _vbitDatabase.reserve(iDataSize);

    //vector<int> viItemLenTemp(iDataDims,0);
    for(int i=0; i<_iDataDims; i++)
    {
        myBitSet<N> bitItemTemp;
        _vbitItemDatabase.push_back(bitItemTemp);
    }
    for(int i=0; i<_iDataSize; i++)
    {
        int tmpTransLength = 0;
        for(int j=0; j<_iDataDims; j++)
        {
            if(_vbitDatabase[i].test(j))
            {
                _vbitItemDatabase[j].set(i);
                tmpTransLength ++;
                //viItemLenTemp[j] ++;
            }
        }
        _viTransLength.push_back(tmpTransLength);
    }
    //cout<<"Init _vbitDatabase"<<endl;
    for(int i=0; i<_iDataDims; i++)
    {
        _viItemLength.push_back(_vbitItemDatabase[i].SSE4_count());
    }
    //_viItemLength = viItemLenTemp;
    //_vsizetSortedItemLengthIndex = fnSortIndex(_viItemLength,0);

    srand((unsigned) time(NULL));
    //cout<<"NSGAII::NSGAII Finish!"<<endl;
}

vector<IndividualNode> NSGAII::_fnInitialization()
{
    vector<IndividualNode> vnodePopulations;
    vnodePopulations.reserve(_iPopSize);

    for(int i=0; i<_iPopSize; i++)
    {
        IndividualNode nodeTmp = IndividualNode();
        vnodePopulations.push_back(nodeTmp);
    }
    //Random Initialization
    //for(int i=0; i<_iPopSize; i++)
    //{
        //for(int j=0; j<_iPopDims; j++)
        //{
            //if((double) rand() / RAND_MAX > 0.5)
            //{
                //vnodePopulations[i]._bitTransaction.set(j);
            //}
        //}
    //}


    //Initialization
    if((_iPopSize >> 1) >= _iPopDims)
    {
        for(int i=0; i<_iPopDims; i++)
        {
            vnodePopulations[i]._bitTransaction.set(i);
        }
        for(int i=_iPopDims; i<_iPopSize; i++)
        {
            int iRandIndex                      = rand() % _iDataSize;
            vnodePopulations[i]._bitTransaction = _vbitDatabase[iRandIndex];
        }
    }
    else
    {
        int iInterSize = _iPopSize >> 1;

        vector<size_t> vsizetSortedItemNo = fnSortIndex(_viItemLength, 1);

        for(int i=0; i<iInterSize; i++)
        {
            vnodePopulations[i]._bitTransaction = myBitSet<M>();
            vnodePopulations[i]._bitTransaction.set(vsizetSortedItemNo[i]);
        }
        for(int i=iInterSize; i<_iPopSize; i++)
        {
            int iRandIndex                      = rand() % _iDataSize;
            vnodePopulations[i]._bitTransaction = myBitSet<M>();
            vnodePopulations[i]._bitTransaction = _vbitDatabase[iRandIndex];
        }
    }

    //cout<<"NSGAII::_fnInitialization Finish!"<<endl;
    return vnodePopulations;
}

void NSGAII::_fnCalcFiteness( vector<IndividualNode> &vnodePopulations )
{
    double fMaxArea = 0, fMinArea = MAX_VALUE;
    for(int i=0; i<_iPopSize; i++)
    {
        vnodePopulations[i]._vfFitness.clear();
        int iCurLength   = vnodePopulations[i]._bitTransaction.SSE4_count();
        double fOccupancy = 0;
        myBitSet<N> bitItemTemp;
        bitItemTemp.set();

        for(int j=0; j<_iPopDims; j++)
        {
            if(vnodePopulations[i]._bitTransaction.test(j))
            {
                //cout<<j<<",";
                bitItemTemp &= _vbitItemDatabase[j];
            }
        }

        //cout<<"-----";
        int iFreqNum   = bitItemTemp.SSE4_count();
        //cout<<iCurLength <<", "<<iFreqNum<<endl;
        double fSupport = (double) iFreqNum / _iDataSize;
        double fArea = (double) fSupport * iCurLength;

        //Save the max and min Area
        if(fArea > fMaxArea)
        {
            fMaxArea = fArea;
        }
        if(fArea < fMinArea)
        {
            fMinArea = fArea;
        }

        if(fSupport > 0)
        {
            //Calculate Occupancy
            for(int j=0; j<_iDataSize; j++)
            {
                if(bitItemTemp.test(j))
                {
                    fOccupancy += (double) iCurLength / _viTransLength[j];
                }
            }
            fOccupancy /= (double) iFreqNum;

            if(fOccupancy > 0)
            {
                //Calculate Area

                vnodePopulations[i]._vfFitness.push_back(fSupport);
                vnodePopulations[i]._vfFitness.push_back(fOccupancy);
                vnodePopulations[i]._vfFitness.push_back(fArea);
            }
            else // Occupancy = 0 means this pattern is Nonsense;
            {
                fMinArea = 0;
                vnodePopulations[i]._vfFitness.push_back(0);
                vnodePopulations[i]._vfFitness.push_back(0);
                vnodePopulations[i]._vfFitness.push_back(0);
            }
        }
        else // Support = 0 means this pattern is Nonsense
        {
            fMinArea = 0;
            vnodePopulations[i]._vfFitness.push_back(0);
            vnodePopulations[i]._vfFitness.push_back(0);
            vnodePopulations[i]._vfFitness.push_back(0);
        }
    }

    // Normalization to Area;
    //cout<<fMaxArea<<":"<<fMinArea<<endl;
    if(!(fMaxArea-fMinArea >= -EPSINON && fMaxArea-fMaxArea <= EPSINON))
    {
        for( int i=0; i<_iPopSize; i++ )
        {
            double fArea = vnodePopulations[i]._vfFitness[2];
            fArea = (fArea - fMinArea) / (fMaxArea - fMinArea);
            vnodePopulations[i]._vfFitness[2] = fArea;
        }
    }
    //cout<<"NSGAII::_fnCalcFiteness Finish!"<<endl;
}

vector< vector< int > > NSGAII::_fnNonDominateSort(vector<IndividualNode> & vnodePopulations)
{
    int iPopSize = vnodePopulations.size();
    vector< vector<int> > vviIDominate(iPopSize);
    vector<int> viDominateMe(iPopSize,0);
    vector< vector<int> > vviFrontList(iPopSize);

    for( auto node : vnodePopulations ){
        node._bIsDuplicate = false;
    }

    for( int i=0; i<iPopSize-1; i++ )
    {
        for( int j=i+1; j<iPopSize; j++ )
        {
            bool bMask0 = true, bMask1 = true;
            for( int k=0; k<_iObjDims; k++ )
            {
                bMask0 &= (vnodePopulations[i]._vfFitness[k] >= vnodePopulations[j]._vfFitness[k]);
                bMask1 &= (vnodePopulations[i]._vfFitness[k] <= vnodePopulations[j]._vfFitness[k]);
            }

            if(bMask0 && !bMask1)
            {
                vviIDominate[i].push_back(j);
                viDominateMe[j] ++;
            }
            else if(!bMask0 && bMask1)
            {
                vviIDominate[j].push_back(i);
                viDominateMe[i] ++;
            }
            else if(bMask0 && bMask1)
            {
                myBitSet<M> bitTransTemp;
                bitTransTemp = vnodePopulations[i]._bitTransaction;
                bitTransTemp = bitTransTemp ^ vnodePopulations[j]._bitTransaction;
                if(bitTransTemp.SSE4_count() == 0)
                {
                    if( !vnodePopulations[i]._bIsDuplicate && !vnodePopulations[j]._bIsDuplicate ){
                        vnodePopulations[j]._bIsDuplicate = true;
                    }
                    else if( vnodePopulations[i]._bIsDuplicate ){
                        vnodePopulations[j]._bIsDuplicate = true;
                    }
                    else{
                        #ifdef _DEBUG_
                        cout<<"\033[41m"<<"Calculate Solution Duplicate Err!\033[0m"<<endl;
                        #endif
                    }
                    vviIDominate[i].push_back(j);
                    viDominateMe[j]++;
                }
            }
        }
    }

    for( int i=0; i<iPopSize; i++ )
    {
        if(viDominateMe[i]==0)
        {
            vviFrontList[0].push_back(i);
            vnodePopulations[i]._iFrontNo = 0;
        }
    }

    int iCurFrontIndex = 0, iCounter=vviFrontList[0].size();
    while(iCounter < iPopSize)
    {
        iCurFrontIndex++;
        for(auto iCurPointIndex:vviFrontList[iCurFrontIndex-1])
        {
            for(auto i:vviIDominate[iCurPointIndex])
            {
                viDominateMe[i]--;
                if(viDominateMe[i] == 0)
                {
                    vviFrontList[iCurFrontIndex].push_back(i);
                    vnodePopulations[i]._iFrontNo = iCurFrontIndex;
                    iCounter++;
                }
            }
        }
    }

    //cout<<"NSGAII::_fnNonDominateSort Finish!"<<endl;
    return vviFrontList;
}

//Caculate Crowd Distance
void NSGAII::_fnCalcCrowdDistance
(
 vector<IndividualNode> & vnodePopulations,
 const vector<vector<int> > &vviFrontList
 )
{
    // Assign 0 to each crowd distance;
    for(auto &i:vnodePopulations)
    {
        i._fCrowdDistance = 0;
    }
    for(auto &iCurFront:vviFrontList)
    {
        int iCurFrontSize = iCurFront.size();
        if(iCurFrontSize == 0)
        {
            break;
        }
        for(int i=0; i<_iObjDims; i++)
        {
            vector<double> vfSingleFitness(iCurFrontSize);
            for(int j=0; j<iCurFrontSize; j++)
            {
                vfSingleFitness[j] = vnodePopulations[iCurFront[j]]._vfFitness[i];
            }
            vector<size_t> vsizetSortedFitnessIndex = fnSortIndex(vfSingleFitness, 0);

            double fMaxDistance = vfSingleFitness[vsizetSortedFitnessIndex[iCurFrontSize-1]];
            double fMinDistance = vfSingleFitness[vsizetSortedFitnessIndex[0]];
            vnodePopulations[iCurFront[vsizetSortedFitnessIndex[0]]]._fCrowdDistance = MAX_VALUE;
            vnodePopulations[iCurFront[vsizetSortedFitnessIndex[iCurFrontSize-1]]]._fCrowdDistance = MAX_VALUE;

            if( fMaxDistance - fMinDistance <= EPSINON && fMaxDistance - fMinDistance >= -EPSINON ){
            }
            else{
                for(int j=1; j<iCurFrontSize-1; j++)
                {
                    int iIdx = iCurFront[vsizetSortedFitnessIndex[j]];
                    if( vnodePopulations[iIdx]._fCrowdDistance < MAX_VALUE )
                        vnodePopulations[iIdx]._fCrowdDistance +=
                            (vnodePopulations[iCurFront[vsizetSortedFitnessIndex[j+1]]]._vfFitness[i] - vnodePopulations[iCurFront[vsizetSortedFitnessIndex[j-1]]]._vfFitness[i]) / (fMaxDistance-fMinDistance);
                }
            }
        }
    }
    //cout<<"NSGAII::_fnCalcCrowdDistance Finish!"<<endl;
}

vector<IndividualNode> NSGAII::_fnSelectMatingPool(const vector<IndividualNode> &vnodePopulations)
{
    vector<IndividualNode> vnodeMatePool(_iPopSize);
    for(int i=0; i<_iPopSize; i++)
    {
        int iParentIndex0 = i;
        int iParentIndex1;
        while(true)
        {
            iParentIndex1 = rand() % _iPopSize;
            if(iParentIndex0 != iParentIndex1)
            {
                break;
            }
        }
        if(vnodePopulations[iParentIndex0]._iFrontNo < vnodePopulations[iParentIndex1]._iFrontNo)
        {
            vnodeMatePool[i] = vnodePopulations[iParentIndex0];
        }
        else if(vnodePopulations[iParentIndex0]._iFrontNo > vnodePopulations[iParentIndex1]._iFrontNo)
        {
            vnodeMatePool[i] = vnodePopulations[iParentIndex1];
        }
        else if(vnodePopulations[iParentIndex0]._fCrowdDistance > vnodePopulations[iParentIndex1]._fCrowdDistance)
        {
            vnodeMatePool[i] = vnodePopulations[iParentIndex0];
        }
        else if(vnodePopulations[iParentIndex0]._fCrowdDistance <= vnodePopulations[iParentIndex1]._fCrowdDistance)
        {
            vnodeMatePool[i] = vnodePopulations[iParentIndex1];
        }
    }

    //cout<<"NSGAII::_fnSelectMatingPool Finish!"<<endl;
    return vnodeMatePool;
}

void NSGAII::_fnReproduceOff(vector<IndividualNode> &vnodeMatePool)
{
    //cout<<"NSGAII::_fnReproduceOff start!"<<endl;
    //cout<<"Cross Over"<<endl;
    for(int i=0; i<_iPopSize; i=i+2)
    {
        for(int j=0; j<_iPopDims; j++)
        {
            myBitSet<M> bitTransTemp0 = vnodeMatePool[i]._bitTransaction;
            myBitSet<M> bitTransTemp1 = vnodeMatePool[i+1]._bitTransaction;
            if((((double) rand() / RAND_MAX) > 0.5) && ( bitTransTemp0[j] ^ bitTransTemp1[j] ))
            {
                //cout<<j<<",";
                vnodeMatePool[i]._bitTransaction.flip(j);
                vnodeMatePool[i+1]._bitTransaction.flip(j);
            }
        }
        //cout<<endl;
    }
    //cout<<"Mutation"<<endl;
    for(int i=0; i<_iPopSize; i=i+1)
    {
        for(int j=0; j<_iPopDims; j++)
        {
            if(((double) rand() / RAND_MAX) < ((double) 1 / _iPopDims) )
            {
                //cout<<j<<",";
                vnodeMatePool[i]._bitTransaction.flip(j);
            }
        }
        //cout<<endl;
    }
    //cout<<"NSGAII::_fnReproduceOff Finish!"<<endl;
}

vector<IndividualNode>  NSGAII::_fnNatureSelection
(
 const vector<IndividualNode> &vnodeAllPop,
 const vector<vector<int> > & vviFrontList
 )
{
    vector<IndividualNode> vnodeNextPop(_iPopSize);
    int iCurFrontIndex = 0;
    int iCurSelectNum  = 0;
    int iCurFrontSize  = vviFrontList[iCurFrontIndex].size();
    while(iCurFrontSize + iCurSelectNum <= _iPopSize)
    {
        for(auto iSelectedIndex : vviFrontList[iCurFrontIndex])
        {
            vnodeNextPop[iCurSelectNum++] = vnodeAllPop[iSelectedIndex];
        }
        iCurFrontIndex++;
        iCurFrontSize = vviFrontList[iCurFrontIndex].size();
    }

    if(iCurSelectNum < _iPopSize)
    {
        vector<double> vfCrowdDis;
        for(auto iSelectedIndex : vviFrontList[iCurFrontIndex])
        {
            vfCrowdDis.push_back(vnodeAllPop[iSelectedIndex]._fCrowdDistance);
        }
        vector<size_t> vsizetSortedDisIndex = fnSortIndex(vfCrowdDis,1);
        for(auto iSelectedIndex : vsizetSortedDisIndex)
        {
            if(iCurSelectNum < _iPopSize)
            {
                vnodeNextPop[iCurSelectNum++] = vnodeAllPop[vviFrontList[iCurFrontIndex][iSelectedIndex]];
            }
        }
    }

    //cout<<"NSGAII::_fnNatureSelection Finish!"<<endl;
    return vnodeNextPop;
}


vector<IndividualNode>  NSGAII::_fnNatureSelectionNoDuplicate
(
 const vector<IndividualNode> &vnodeAllPop,
 const vector<vector<int> > & vviFrontList
 )
{
    vector<IndividualNode> vnodeNextPop(_iPopSize);
    int iCurFrontIndex = 0;
    int iCurSelectNum  = 0;
    int iDupCnt        = 0;
    for( auto node : vnodeAllPop ){
        if( node._bIsDuplicate ){
            iDupCnt ++;
        }
    }

    if( iDupCnt < _iPopSize ){
        int iCurFrontSize  = vviFrontList[iCurFrontIndex].size();
        for( auto idx : vviFrontList[iCurFrontIndex] ){
            if(vnodeAllPop[idx]._bIsDuplicate){
                iCurFrontSize--;
            }
        }
        while(iCurFrontSize + iCurSelectNum <= _iPopSize)
        {
            for(auto iSelectedIndex : vviFrontList[iCurFrontIndex])
            {
                if( !vnodeAllPop[iSelectedIndex]._bIsDuplicate )
                    vnodeNextPop[iCurSelectNum++] = vnodeAllPop[iSelectedIndex];
            }
            iCurFrontIndex++;
            iCurFrontSize = vviFrontList[iCurFrontIndex].size();
            for( auto idx : vviFrontList[iCurFrontIndex] ){
                if(vnodeAllPop[idx]._bIsDuplicate){
                    iCurFrontSize--;
                }
            }
        }

        if(iCurSelectNum < _iPopSize)
        {
            vector<double> vfCrowdDis;
            for(auto iSelectedIndex : vviFrontList[iCurFrontIndex])
            {
                vfCrowdDis.push_back(vnodeAllPop[iSelectedIndex]._fCrowdDistance);
            }
            vector<size_t> vsizetSortedDisIndex = fnSortIndex(vfCrowdDis,1);
            for(auto iSelectedIndex : vsizetSortedDisIndex)
            {
                if(iCurSelectNum < _iPopSize && !vnodeAllPop[vviFrontList[iCurFrontIndex][iSelectedIndex]]._bIsDuplicate)
                {
                    vnodeNextPop[iCurSelectNum++] = vnodeAllPop[vviFrontList[iCurFrontIndex][iSelectedIndex]];
                }
            }
        }
    }
    else{
        for( auto node : vnodeAllPop ){
            if( !node._bIsDuplicate )
                vnodeNextPop[iCurSelectNum++] = node;
        }

        int iCurFrontSize  = vviFrontList[iCurFrontIndex].size();
        for( auto idx : vviFrontList[iCurFrontIndex] ){
            if(!vnodeAllPop[idx]._bIsDuplicate){
                iCurFrontSize--;
            }
        }
        while(iCurFrontSize + iCurSelectNum <= _iPopSize)
        {
            for(auto iSelectedIndex : vviFrontList[iCurFrontIndex])
            {
                if( vnodeAllPop[iSelectedIndex]._bIsDuplicate )
                    vnodeNextPop[iCurSelectNum++] = vnodeAllPop[iSelectedIndex];
            }
            iCurFrontIndex++;
            iCurFrontSize = vviFrontList[iCurFrontIndex].size();
            for( auto idx : vviFrontList[iCurFrontIndex] ){
                if(!vnodeAllPop[idx]._bIsDuplicate){
                    iCurFrontSize--;
                }
            }
        }

        if(iCurSelectNum < _iPopSize)
        {
            vector<double> vfCrowdDis;
            for(auto iSelectedIndex : vviFrontList[iCurFrontIndex])
            {
                vfCrowdDis.push_back(vnodeAllPop[iSelectedIndex]._fCrowdDistance);
            }
            vector<size_t> vsizetSortedDisIndex = fnSortIndex(vfCrowdDis,1);
            for(auto iSelectedIndex : vsizetSortedDisIndex)
            {
                if(iCurSelectNum < _iPopSize && vnodeAllPop[vviFrontList[iCurFrontIndex][iSelectedIndex]]._bIsDuplicate)
                {
                    vnodeNextPop[iCurSelectNum++] = vnodeAllPop[vviFrontList[iCurFrontIndex][iSelectedIndex]];
                }
            }
        }
    }

    //cout<<"NSGAII::_fnNatureSelection Finish!"<<endl;
    return vnodeNextPop;
}

bool NSGAII::_fnCheckSimilar( vector<IndividualNode> lhs, vector<IndividualNode> rhs ){
    if( rhs.empty() ) return false;
    auto less = []( IndividualNode i1, IndividualNode i2 ){
            double d1 = i1._vfFitness[0], d2 = i2._vfFitness[0];
            bool res;
            if( d1-d2 >= -EPSINON && d1-d2 <= EPSINON )
                res = i1._vfFitness[1] < i2._vfFitness[1];
            else
                res = d1 < d2;
            return res;
    };
    sort( lhs.begin(), lhs.end(), less );
    sort( rhs.begin(), rhs.end(), less );
    const int size = lhs.size();
    int cnt = 0;
    for( int i=0, j=0; i<size&& j<size; ){
        double ld0 = lhs[i]._vfFitness[0], ld1 = lhs[i]._vfFitness[1];
        double rd0 = rhs[j]._vfFitness[0], rd1 = rhs[j]._vfFitness[1];
        if( ld0 - rd0 >= -EPSINON && ld0 - rd0 <= EPSINON ){
            if( ld1-rd1 >= -EPSINON && ld1-rd1 <= EPSINON )
            {
                cnt++;
                i++, j++;
            }
            else if( ld1 > rd1 ){
                j++;
            }
            else if( ld1 < rd1 ){
                i++;
            }
        }
        else if( ld0 > rd0 ){
            j++;
        }
        else if( ld0 < rd0 ){
            i++;
        }
    }
    double rate = (double) cnt / size;
    #ifdef _DUBUG_
    cout<<cnt<<"/"<<size << endl;
    #endif
    if( rate > 0.9 )
    {
        return true;
    }
    return false;
}

void NSGAII::_fnDebugPrintInfo( ofstream &logs, vector<IndividualNode> & vnodePopulations ){
    for( auto & i:vnodePopulations ){
        for( int j=0; j<_iPopDims; ++j ){
            if( i._bitTransaction.test(j) ){
                logs<<j<<",";
            }
        }
        logs<<"\b : sup="<<i._vfFitness[0]<<",occ="<<i._vfFitness[1]<<",are="<<i._vfFitness[2]<<endl;
    }
}

vector<IndividualNode>  NSGAII::_fnMOEC( int & traversNode, double & spendTime )
{
    //for debug
    //ofstream logfile("./logs");

    clock_t start, end, cmpTime, cstart, cend;
    start = clock();

    vector<IndividualNode> vnodePopulations = _fnInitialization();
    _fnCalcFiteness( vnodePopulations );
    vector<vector<int> > vviFrontList = _fnNonDominateSort( vnodePopulations );
    _fnCalcCrowdDistance( vnodePopulations, vviFrontList );

    vector<IndividualNode> lastPop;
    int cnt = 0;
    int iGene;
    cmpTime = 0;
    for( iGene=0; iGene<50; iGene++ )
    {
        //for debug
        #ifdef _DEBUG_
        //logfile<<endl;
        //logfile<<"Info "<<iGene<<" th Populations"<<endl;
        //_fnDebugPrintInfo(logfile, vnodePopulations);
        #endif

        vector<IndividualNode> vnodeChildPop = _fnSelectMatingPool( vnodePopulations );
        _fnReproduceOff( vnodeChildPop );
        _fnCalcFiteness( vnodeChildPop );
        vector<IndividualNode> vnodeMixedPop;
        vnodeMixedPop.reserve( 2*_iPopSize );
        vnodeMixedPop.insert( vnodeMixedPop.end(),vnodePopulations.begin(),vnodePopulations.end() );
        vnodeMixedPop.insert( vnodeMixedPop.end(),vnodeChildPop.begin(),vnodeChildPop.end() );
        vector<vector<int>> vviMixFrontList = _fnNonDominateSort( vnodeMixedPop );

        _fnCalcCrowdDistance( vnodeMixedPop, vviMixFrontList );

        //for debug
        #ifdef _DEBUG_
        int iFrontCnt        = 0;
        int iIndCnt          = 0;
        int iBackgroundColor = 0;
        for( auto vFrontlist : vviMixFrontList ){
            iBackgroundColor = iFrontCnt % 8;
            for( auto iNo : vFrontlist ){
                cout<<"\033[4"<<iBackgroundColor<<";32;1m"<<setw(3)<<iFrontCnt+1<<"th Front "<<setw(3)<<++iIndCnt<<":";
                for( int i=0; i<M; ++i ){
                    if( vnodeMixedPop[iNo]._bitTransaction.test(i) ){
                        cout<<"\033[4"<<iBackgroundColor<<";35;1m"<<setw(2)<<i<<" \033[0m";
                    }
                }
                cout<<"\033[4"<<iBackgroundColor<<";36;1m"<<":"<<vnodeMixedPop[iNo]._vfFitness[0]<<","<<vnodeMixedPop[iNo]._vfFitness[1]<<","<<vnodeMixedPop[iNo]._vfFitness[2]<<"\033[0m";
                cout<<"\033[4"<<iBackgroundColor<<";37;1m"<<":"<<vnodeMixedPop[iNo]._fCrowdDistance<<"\033[0m";
                if( vnodeMixedPop[iNo]._bIsDuplicate )
                    cout<<"\033[41;34;1m"<<":Duplicate\033[0m"<<endl;
                else
                    cout<<"\033[4"<<iBackgroundColor<<";34;1m"<<":Single\033[0m"<<endl;
            }
            iFrontCnt++;
        }
        //end debug
        #endif

        vnodePopulations.clear();
        vnodePopulations = _fnNatureSelectionNoDuplicate(vnodeMixedPop, vviMixFrontList);

        //cstart = clock();
        if(_fnCheckSimilar( vnodePopulations, lastPop)){
            cnt++;
        }
        if( cnt > 5 )
            break;
        //cend = clock();
        //cmpTime += cend - cstart;

        lastPop.clear();
        lastPop = vnodePopulations;
        cout<<iGene<<"th iterators finished ..."<<endl;
    }

    end = clock();
    spendTime = (double)(end-start) / CLOCKS_PER_SEC;

    traversNode = _iPopSize * iGene;
    vector<IndividualNode> vnodeOutput;
    vnodeOutput.reserve(_iPopSize);
    for(const auto &i:vnodePopulations)
    {
        if(i._iFrontNo == 0 && i._vfFitness[0] > (double) 1 /N )
        {
            vnodeOutput.push_back(i);
        }
    }
    return vnodeOutput;
}

