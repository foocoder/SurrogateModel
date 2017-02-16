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
#include "RadialBasisFunction.h"

using namespace std;

#define EPSINON 1e-18

template <typename T>
vector<size_t> fnSortIndex
(
 const vector<T> & v,
 int iType
 )
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

// 初始化工作
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

// 训练代理模型
vector<RadialBasisFunction> NSGAII::_fnBuildModel
(
 const vector<IndividualNode> & vnodeDatabase
 ){
    int iDBSize = vnodeDatabase.size();
    vector< myBitSet<M> > vviSample;
    vviSample.reserve( iDBSize );
    vector<double> vfRealSup( iDBSize );
    vector<double> vfRealOcc( iDBSize );
    vector<double> vfRealAre( iDBSize );

    for( int i=0; i<iDBSize; ++i ){
        /* vector<int> viSolution( _iPopDims ); */
        /* for( int j=0; j<_iPopDims; ++j ){ */
        /*     viSolution[j] = ( vnodeDatabase[i]._bitTransaction.test(j) ? 1 : 0 ); */
        /* } */
        vviSample.push_back( vnodeDatabase[i]._bitTransaction );
        /* vviSample[i] = viSolution; */
        vfRealSup[i] = vnodeDatabase[i]._vfFitness[0];
        vfRealOcc[i] = vnodeDatabase[i]._vfFitness[1];
        vfRealAre[i] = vnodeDatabase[i]._vfFitness[2];
    }

    RadialBasisFunction rbfSup( iDBSize, iDBSize, _iPopDims, vviSample, vfRealSup );
    RadialBasisFunction rbfOcc( iDBSize, iDBSize, _iPopDims, vviSample, vfRealOcc );
    RadialBasisFunction rbfAre( iDBSize, iDBSize, _iPopDims, vviSample, vfRealAre );

    rbfSup.runRBF();
    rbfOcc.runRBF();
    rbfAre.runRBF();
    return { rbfSup, rbfOcc, rbfAre };
}

//利用代理模型评估个体适应度
void NSGAII::_fnCalcEstimation
(
 const vector<RadialBasisFunction> & vmRBFModels,
 vector<IndividualNode> & vnodePopulations
 ){
    for( int i=0; i<_iPopSize; ++i ){
        /* vector<int> viSolution( _iPopDims ); */
        /* for( int j=0; j<_iPopDims; ++j ){ */
        /*     viSolution[j] = ( vnodePopulations[i]._bitTransaction.test(j) ? 1 : 0 ); */
        /* } */
        vnodePopulations[i]._vfFitness[0]
            = vmRBFModels[0].getEstimation( vnodePopulations[i]._bitTransaction, RadialBasisFunction::GetKernalType("Gaussian") );
        vnodePopulations[i]._vfFitness[1]
            = vmRBFModels[1].getEstimation( vnodePopulations[i]._bitTransaction, RadialBasisFunction::GetKernalType("Gaussian") );
        vnodePopulations[i]._vfFitness[2]
            = vmRBFModels[2].getEstimation( vnodePopulations[i]._bitTransaction, RadialBasisFunction::GetKernalType("Gaussian") );
    }
}

// Local Search Phase
/* void NSGAII::_fnLocalSearchPhase */
/* ( */
/*   vector<IndividualNode> & vnodePopulations */
/*   /1* const vector<IndividualNode> & vnodeGlobalDB *1/ */
/*   ){ */

/*     for( int i=0; i<_iPopSize; ++i ){ */
/*         //生成随机权值 */
/*         vector<double> vdAggrWeight( 3 ); */
/*         vdAggrWeight[0] = 1 - sqrt( static_cast<double>(rand())/RAND_MAX ); */
/*         vdAggrWeight[1] = ( 1 - vdAggrWeight[0] ) * ( 1 - static_cast<double>(rand())/RAND_MAX ); */
/*         vdAggrWeight[2] = 1 - vdAggrWeight[0] - vdAggrWeight[1]; */

/*         // 对当前个体局部搜索 */
/*         vnodePopulations[i] = _fnFindLocalOptima( vnodePopulations[i], vdAggrWeight ); */

/*     } */
/* } */


// Local Search Phase
/* void NSGAII::_fnLocalSearchPhase */
/* ( */
/*   vector<IndividualNode> & vnodePopulations, */
/*   const vector<IndividualNode> & vnodeGlobalDB */
/*   ){ */

/*     // 生成 vector<int> 类型 GlobalDB */
/*     int iDBSize = vnodeGlobalDB.size(); */
/*     vector<vector<int>> vviGlobalDB( iDBSize ); */
/*     for( int i=0; i<iDBSize; ++i ){ */
/*         vector<int> viSolution( _iPopDims ); */
/*         for( int j=0; j<_iPopDims; ++j ){ */
/*             viSolution[j] = ( vnodeGlobalDB[i]._bitTransaction.test(j) ? 1 : 0 ); */
/*         } */
/*         vviGlobalDB[i] = viSolution; */
/*     } */

/*     for( int i=0; i<_iPopSize; ++i ){ */
/*         //生成随机权值 */
/*         vector<double> vdAggrWeight( 3 ); */
/*         vdAggrWeight[0] = 1 - sqrt( static_cast<double>(rand())/RAND_MAX ); */
/*         vdAggrWeight[1] = ( 1 - vdAggrWeight[0] ) * ( 1 - static_cast<double>(rand())/RAND_MAX ); */
/*         vdAggrWeight[2] = 1 - vdAggrWeight[0] - vdAggrWeight[1]; */

/*         #ifdef _DEBUG1_ */
/*         cout<<setw(10)<<vdAggrWeight[0] */
/*             <<setw(10)<<vdAggrWeight[1] */
/*             <<setw(10)<<vdAggrWeight[2]<<endl; */
/*         #endif */

/*         //计算加权聚合函数值 */
/*         vector<double> vdAggrValue( iDBSize ); */
/*         for( int j=0; j<iDBSize; ++j ){ */
/*             vdAggrValue[j] = vdAggrWeight[0] * vnodeGlobalDB[j]._vfFitness[0] */
/*                 + vdAggrWeight[1] * vnodeGlobalDB[j]._vfFitness[1] */
/*                 + vdAggrWeight[2] * vnodeGlobalDB[j]._vfFitness[2]; */
/*         } */

/*         // 利用GlobalDB对聚合函数建模 */
/*         RadialBasisFunction rbfAggr ( iDBSize, 10, _iPopDims, vviGlobalDB, vdAggrValue ); */
/*         rbfAggr.runRBF(); */

/*         // 计算Ensemble 权值 */
/*         vector<double> vdLambda( 3 ); */
/*         double fErrGauss = rbfAggr.getRMSE( RadialBasisFunction::GetKernalType("Gaussian") ); */
/*         double fErrMulti = rbfAggr.getRMSE( RadialBasisFunction::GetKernalType("MultiQuadratic") ); */
/*         double fErrInvrs = rbfAggr.getRMSE( RadialBasisFunction::GetKernalType("InverseMultiQuadratic") ); */
/*         vdLambda[0] = ( fErrMulti+fErrInvrs ) / ( 2 * (fErrGauss+fErrMulti+fErrInvrs) ); */
/*         vdLambda[1] = ( fErrGauss+fErrInvrs ) / ( 2 * (fErrGauss+fErrMulti+fErrInvrs) ); */
/*         vdLambda[2] = ( fErrMulti+fErrGauss ) / ( 2 * (fErrGauss+fErrMulti+fErrInvrs) ); */

/*         // 计算个体0-1编码 */
/*         vector<int> viInd( _iPopDims ); */
/*         for( int j=0; j<_iPopDims; ++j ){ */
/*             viInd[j] = (vnodePopulations[i]._bitTransaction.test(j) ? 1 : 0); */
/*         } */

/*         // 对当前个体局部搜索 */
/*         vector<int> viOptima = _fnFindLocalOptima( viInd, rbfAggr, vdLambda ); */

/*         // 对局部最优个体构建Node */
/*         IndividualNode nodeOptima; */
/*         for( int j=0; j<_iPopDims; ++j ){ */
/*             if( viOptima[j] ) */
/*                 nodeOptima._bitTransaction.set(j); */
/*         } */

/*         // 计算Opt 以及 Origin 适应度 */
/*         _fnCalcFiteness( nodeOptima ); */
/*         _fnCalcFiteness( vnodePopulations[i] ); */

/*         // Strategy 1: MOGLS */
/*         double fOriginAggr = vdAggrWeight[0] * vnodePopulations[i]._vfFitness[0] */
/*             + vdAggrWeight[1] * vnodePopulations[i]._vfFitness[1] */
/*             + vdAggrWeight[2] * vnodePopulations[i]._vfFitness[2]; */
/*         double fOptAggr = vdAggrWeight[0] * nodeOptima._vfFitness[0] */
/*             + vdAggrWeight[1] * nodeOptima._vfFitness[1] */
/*             + vdAggrWeight[2] * nodeOptima._vfFitness[2]; */

/*         if( fOriginAggr < fOptAggr ) */
/*             vnodePopulations[i] = nodeOptima; */

/*         // 计算Opt 与Origin 支配关系 */
/*         /1* bool bMask0 = true, bMask1 = true; *1/ */
/*         /1* for( int j=0; j<_iObjDims; ++j ){ *1/ */
/*         /1*     bMask0 &= (nodeOptima._vfFitness[j] >= vnodePopulations[i]._vfFitness[j]); *1/ */
/*         /1*     bMask1 &= (nodeOptima._vfFitness[j] <= vnodePopulations[i]._vfFitness[j]); *1/ */
/*         /1* } *1/ */

/*         /1* if( bMask0 && !bMask1 ){ // Opt 支配 Origin *1/ */
/*         /1*     vnodePopulations[i] = nodeOptima; *1/ */
/*         /1* } *1/ */
/*         /1* else if( !bMask0 && bMask1 ){ // Origin 支配 Opt *1/ */

/*         /1* } *1/ */
/*         /1* else if( bMask0 && bMask1 ){ // Origin 与 Opt 相等 *1/ */

/*         /1* } *1/ */
/*         /1* else{ // Opt 与 Origin 非支配 *1/ */
/*         /1* } *1/ */

/*     } */

/* } */

// Find Local Optima Using Local Search
/* IndividualNode NSGAII::_fnFindLocalOptima */
/* ( */
/*   IndividualNode nodeInd, */
/*   const vector<double> & vdAggrWeight */
/*   ){ */
/*     double fOptimaValue = vdAggrWeight[0] * nodeInd._vfFitness[0] */
/*         + vdAggrWeight[1] * nodeInd._vfFitness[1] */
/*         + vdAggrWeight[2] * nodeInd._vfFitness[2]; */

/*     int iOptimaIdx = -1; */

/*     for( int i=0; i<_iPopDims; ++i ){ */
/*         nodeInd._bitTransaction.flip(i); */
/*         _fnCalcFiteness( nodeInd ); */

/*         double fNeighborOptimaValue = vdAggrWeight[0] * nodeInd._vfFitness[0] */
/*         + vdAggrWeight[1] * nodeInd._vfFitness[1] */
/*         + vdAggrWeight[2] * nodeInd._vfFitness[2]; */

/*         if( fNeighborOptimaValue > fOptimaValue ) */
/*             iOptimaIdx = i; */
/*         nodeInd._bitTransaction.flip(i); */
/*     } */
/*     if( iOptimaIdx != -1 ) */
/*     { */
/*         nodeInd._bitTransaction.flip(iOptimaIdx); */
/*         _fnCalcFiteness( nodeInd ); */
/*     } */
/*     return nodeInd; */
/* } */

/* // Find Local Optima Using Local Search */
/* vector<int> NSGAII::_fnFindLocalOptima */
/* ( */
/*   vector<int> viOriginInd, */
/*   const RadialBasisFunction & rbfModel, */
/*   const vector<double> & vdLambda */
/*   ){ */
/*     double fOptimaValue = _fnGetEnsembleEstimaion( viOriginInd, rbfModel, vdLambda ); */
/*     int iOptimaIdx = -1; */
/*     for( int i=0; i<_iPopDims; ++i ){ */
/*         viOriginInd[i] ^= 1; */
/*         double fNeighborOptimaValue */
/*             = _fnGetEnsembleEstimaion( viOriginInd, rbfModel, vdLambda ); */
/*         if( fNeighborOptimaValue > fOptimaValue ) */
/*             iOptimaIdx = i; */
/*         viOriginInd[i] ^= 1; */
/*     } */
/*     if( iOptimaIdx != -1 ) viOriginInd[iOptimaIdx] ^= 1; */
/*     return viOriginInd; */
/* } */

// Get Ensemble RBF Model Estimations
/* inline double NSGAII::_fnGetEnsembleEstimaion */
/* ( */
/*   const vector<int> & viInd, */
/*   const RadialBasisFunction & rbfModel, */
/*   const vector<double> & vdLambda */
/*   ){ */
/*     return vdLambda[0] * rbfModel.getEstimation( viInd, RadialBasisFunction::GetKernalType("Gaussian") ) */
/*          + vdLambda[1] * rbfModel.getEstimation( viInd, RadialBasisFunction::GetKernalType("MultiQuadratic") ) */
/*          + vdLambda[2] * rbfModel.getEstimation( viInd, RadialBasisFunction::GetKernalType("InverseMultiQuadratic") ); */
/* } */

// 计算单个个体真实适应度值
void NSGAII::_fnCalcFiteness
(
 IndividualNode & nodeInd
 ){

    myBitSet<N> bitItemTemp;
    bitItemTemp.set();
    nodeInd._vfFitness.clear();

    for( int i=0; i<_iPopDims; ++i ){
        if( nodeInd._bitTransaction.test(i) ){
            bitItemTemp &= _vbitItemDatabase[i];
        }
    }
    int iFreqNum      = bitItemTemp.SSE4_count();
    int iCurLength    = nodeInd._bitTransaction.SSE4_count();
    double fSupport   = static_cast<double>( iFreqNum ) / _iDataSize;
    double fArea      = fSupport * iCurLength;
    double fOccupancy = 0;
    if( fSupport > 0 ){
        for( int i=0; i<_iDataSize; ++i ){
            if( bitItemTemp.test(i) )
                fOccupancy += static_cast<double>(iCurLength) / _viTransLength[i];
        }
        fOccupancy /= static_cast<double>(iFreqNum);
        if( fOccupancy > 0 ){
            nodeInd._vfFitness.push_back( fSupport );
            nodeInd._vfFitness.push_back( fOccupancy );
            nodeInd._vfFitness.push_back( fArea );
        }
        else{
            nodeInd._vfFitness.push_back( 0 );
            nodeInd._vfFitness.push_back( 0 );
            nodeInd._vfFitness.push_back( 0 );
        }
    }
    else{
        nodeInd._vfFitness.push_back( 0 );
        nodeInd._vfFitness.push_back( 0 );
        nodeInd._vfFitness.push_back( 0 );
    }
}

// 计算个体真实适应度
void NSGAII::_fnCalcFiteness
(
 vector<IndividualNode> &vnodePopulations
 )
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
    /* if(!(fMaxArea-fMinArea >= -EPSINON && fMaxArea-fMinArea <= EPSINON)) */
    /* { */
    /*     for( int i=0; i<_iPopSize; i++ ) */
    /*     { */
    /*         double fArea = vnodePopulations[i]._vfFitness[2]; */
    /*         fArea = (fArea - fMinArea) / (fMaxArea - fMinArea); */
    /*         vnodePopulations[i]._vfFitness[2] = fArea; */
    /*     } */
    /* } */
    //cout<<"NSGAII::_fnCalcFiteness Finish!"<<endl;
}

vector< vector< int > > NSGAII::_fnNonDominateSort
(
 vector<IndividualNode> & vnodePopulations
 )
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
            vnodePopulations[iCurFront[vsizetSortedFitnessIndex[0]]]._fCrowdDistance               = MAX_VALUE;
            vnodePopulations[iCurFront[vsizetSortedFitnessIndex[iCurFrontSize-1]]]._fCrowdDistance = MAX_VALUE;

            if( fMaxDistance - fMinDistance <= EPSINON && fMaxDistance - fMinDistance >= -EPSINON ){
            }
            else{
                for(int j=1; j<iCurFrontSize-1; j++)
                {
                    int iIdx = iCurFront[vsizetSortedFitnessIndex[j]];
                    if( vnodePopulations[iIdx]._fCrowdDistance < MAX_VALUE )
                        vnodePopulations[iIdx]._fCrowdDistance +=
                            (vnodePopulations[iCurFront[vsizetSortedFitnessIndex[j+1]]]._vfFitness[i]
                             - vnodePopulations[iCurFront[vsizetSortedFitnessIndex[j-1]]]._vfFitness[i])
                             / (fMaxDistance-fMinDistance);
                }
            }
        }
    }
    //cout<<"NSGAII::_fnCalcCrowdDistance Finish!"<<endl;
}

vector<IndividualNode> NSGAII::_fnSelectMatingPool
(
 const vector<IndividualNode> &vnodePopulations
 )
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

void NSGAII::_fnReproduceOff
(
 vector<IndividualNode> &vnodeMatePool
 )
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

bool NSGAII::_fnCheckSimilar
(
 vector<IndividualNode> lhs,
 vector<IndividualNode> rhs
 ){
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

void NSGAII::_fnDebugPrintInfo
(
 ofstream &logs, vector<IndividualNode> & vnodePopulations
 ){
    for( auto & i:vnodePopulations ){
        for( int j=0; j<_iPopDims; ++j ){
            if( i._bitTransaction.test(j) ){
                logs<<j<<",";
            }
        }
        logs<<"\b : sup="<<i._vfFitness[0]<<",occ="<<i._vfFitness[1]<<",are="<<i._vfFitness[2]<<endl;
    }
}

//MOEC 利用真实适应度
/* vector<IndividualNode>  NSGAII::_fnMOEC0 */
/* ( */
/*  int & traversNode, */
/*  double & spendTime */
/*  ) */
/* { */
/*     clock_t start, end, cmpTime, cstart, cend; */
/*     start = clock(); */

/*     // 初始化生成GlobalDB */
/*     vector<IndividualNode> vnodePopulations = _fnInitialization(); */
/*     _fnCalcFiteness( vnodePopulations ); */

/*     /1* vector<IndividualNode> lastPop; *1/ */
/*     int cnt = 0; */
/*     int iGene; */
/*     cmpTime = 0; */
/*     for( iGene=0; iGene<5; iGene++ ) */
/*     { */

/*         _fnLocalSearchPhase( vnodePopulations ); */

/*         #ifdef _DEBUG1_ */
/*         cout<<iGene<<"th iterators finished ..."<<endl; */
/*         #endif */
/*     } */

/*     end = clock(); */
/*     spendTime = (double)(end-start) / CLOCKS_PER_SEC; */

/*     vector<vector<int> > vviFrontList = _fnNonDominateSort( vnodePopulations ); */

/*     traversNode = _iPopSize * iGene; */
/*     vector<IndividualNode> vnodeOutput; */
/*     vnodeOutput.reserve(_iPopSize); */
/*     for(const auto &i:vnodePopulations) */
/*     { */
/*         if(i._iFrontNo == 0 && i._vfFitness[0] > (double) 1 /N ) */
/*         { */
/*             vnodeOutput.push_back(i); */
/*         } */
/*     } */
/*     return vnodeOutput; */
/* } */

// MOEC 利用代理模型 GLS
/* vector<IndividualNode>  NSGAII::_fnMOEC1 */
/* ( */
/*  int & traversNode, */
/*  double & spendTime */
/*  ) */
/* { */
/*     //for debug */
/*     //ofstream logfile("./logs"); */

/*     clock_t start, end, cmpTime, cstart, cend; */
/*     start = clock(); */

/*     // 初始化生成GlobalDB */
/*     vector<IndividualNode> vnodeGlobalDB = _fnInitialization(); */
/*     _fnCalcFiteness( vnodeGlobalDB ); */

/*     vector<IndividualNode> vnodePopulations( _iPopSize ); */

/*     /1* vector<IndividualNode> lastPop; *1/ */
/*     int cnt = 0; */
/*     int iGene; */
/*     cmpTime = 0; */
/*     for( iGene=0; iGene<5; iGene++ ) */
/*     { */

/*         vector<vector<int> > vviFrontList = _fnNonDominateSort( vnodeGlobalDB ); */
/*         _fnCalcCrowdDistance( vnodeGlobalDB, vviFrontList ); */

/*         vnodePopulations = _fnSelectMatingPool( vnodeGlobalDB ); */
/*         _fnReproduceOff( vnodePopulations ); */

/*         _fnLocalSearchPhase( vnodePopulations, vnodeGlobalDB ); */
/*         vnodeGlobalDB = vnodePopulations; */

/*         /1* _fnSMEC( vmRBFModels, vnodePopulations ); *1/ */
/*         /1* _fnCalcFiteness( vnodePopulations ); *1/ */
/*         /1* vector<vector<int> > vviFrontList = _fnNonDominateSort( vnodePopulations ); *1/ */
/*         /1* _fnCalcCrowdDistance( vnodePopulations, vviFrontList ); *1/ */
/*         /1* vmRBFModels = _fnBuildModel( vnodePopulations ); *1/ */

/*         /1* cstart = clock(); *1/ */
/*         /1* if(_fnCheckSimilar( vnodePopulations, lastPop)){ *1/ */
/*         /1*     cnt++; *1/ */
/*         /1* } *1/ */
/*         /1* if( cnt > 5 ) *1/ */
/*         /1*     break; *1/ */
/*         /1* cend = clock(); *1/ */
/*         /1* cmpTime += cend - cstart; *1/ */

/*         /1* lastPop.clear(); *1/ */
/*         /1* lastPop = vnodePopulations; *1/ */
/*         #ifdef _DEBUG1_ */
/*         cout<<iGene<<"th iterators finished ..."<<endl; */
/*         #endif */
/*     } */

/*     end = clock(); */
/*     spendTime = (double)(end-start) / CLOCKS_PER_SEC; */

/*     /1* vector<vector<int> > vviFrontList = _fnNonDominateSort( vnodePopulations ); *1/ */

/*     /1* traversNode = _iPopSize * iGene; *1/ */
/*     /1* vector<IndividualNode> vnodeOutput; *1/ */
/*     /1* vnodeOutput.reserve(_iPopSize); *1/ */
/*     /1* for(const auto &i:vnodePopulations) *1/ */
/*     /1* { *1/ */
/*     /1*     if(i._iFrontNo == 0 && i._vfFitness[0] > (double) 1 /N ) *1/ */
/*     /1*     { *1/ */
/*     /1*         vnodeOutput.push_back(i); *1/ */
/*     /1*     } *1/ */
/*     /1* } *1/ */
/*     return vnodePopulations; */
/* } */

// MOEC 利用代理模型 NSGAII
vector<IndividualNode>  NSGAII::_fnMOEC
(
 int & traversNode,
 double & spendTime
 )
{
    //for debug
    //ofstream logfile("./logs");

    clock_t start, end, cmpTime, cstart, cend;
    start = clock();

    vector<IndividualNode> vnodePopulations = _fnInitialization();
    _fnCalcFiteness( vnodePopulations );
    vector<vector<int> > vviFrontList = _fnNonDominateSort( vnodePopulations );
    _fnCalcCrowdDistance( vnodePopulations, vviFrontList );

    vector<RadialBasisFunction> vmRBFModels = _fnBuildModel( vnodePopulations );

    /* vector<IndividualNode> lastPop; */
    int cnt = 0;
    int iGene;
    cmpTime = 0;
    for( iGene=0; iGene<10; iGene++ )
    {

        _fnSMEC( vmRBFModels, vnodePopulations );
        _fnCalcFiteness( vnodePopulations );
        vector<vector<int> > vviFrontList = _fnNonDominateSort( vnodePopulations );
        _fnCalcCrowdDistance( vnodePopulations, vviFrontList );
        vmRBFModels = _fnBuildModel( vnodePopulations );

        /* cstart = clock(); */
        /* if(_fnCheckSimilar( vnodePopulations, lastPop)){ */
        /*     cnt++; */
        /* } */
        /* if( cnt > 5 ) */
        /*     break; */
        /* cend = clock(); */
        /* cmpTime += cend - cstart; */

        /* lastPop.clear(); */
        /* lastPop = vnodePopulations; */
        #ifdef _DEBUG1_
        cout<<iGene<<"th iterators finished ..."<<endl;
        #endif
    }

    end = clock();
    spendTime = (double)(end-start) / CLOCKS_PER_SEC;

    _fnCalcFiteness( vnodePopulations );
    vviFrontList = _fnNonDominateSort( vnodePopulations );

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
    return vnodePopulations;
}

void NSGAII::_fnSMEC
(
 const vector<RadialBasisFunction> & vmRBFModels,
 vector<IndividualNode> & vnodePopulations
 )
{

    /* vector<IndividualNode> lastPop; */
    int cnt = 0;
    int iGene;
    for( iGene=0; iGene<20; iGene++ )
    {

        // 选择交配池
        vector<IndividualNode> vnodeChildPop = _fnSelectMatingPool( vnodePopulations );
        // 产生子代个体
        _fnReproduceOff( vnodeChildPop );
        // 使用代理模型估计个体目标值
        _fnCalcEstimation( vmRBFModels, vnodeChildPop );
        // 将子代与父代混合一起
        vector<IndividualNode> vnodeMixedPop;
        vnodeMixedPop.reserve( 2*_iPopSize );
        vnodeMixedPop.insert( vnodeMixedPop.end(),vnodePopulations.begin(),vnodePopulations.end() );
        vnodeMixedPop.insert( vnodeMixedPop.end(),vnodeChildPop.begin(),vnodeChildPop.end() );
        // 对混合种群进行非支配排序
        vector<vector<int>> vviMixFrontList = _fnNonDominateSort( vnodeMixedPop );
        // 对排序后的种群计算拥挤距离
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
                cout<<"\033[4"<<iBackgroundColor<<";36;1m"<<":"
                    <<vnodeMixedPop[iNo]._vfFitness[0]<<","
                    <<vnodeMixedPop[iNo]._vfFitness[1]<<","
                    <<vnodeMixedPop[iNo]._vfFitness[2]<<"\033[0m";
                cout<<"\033[4"<<iBackgroundColor<<";37;1m"<<":"
                    <<vnodeMixedPop[iNo]._fCrowdDistance<<"\033[0m";
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
        // 对混合种群进行环境选择
        vnodePopulations = _fnNatureSelectionNoDuplicate(vnodeMixedPop, vviMixFrontList);

        /* if(_fnCheckSimilar( vnodePopulations, lastPop)){ */
        /*     cnt++; */
        /* } */
        /* if( cnt > 5 ) */
        /*     break; */

        /* lastPop.clear(); */
        /* lastPop = vnodePopulations; */
        #ifdef _DEBUG_
        cout<<iGene<<"th surrogate model iterators finished ..."<<endl;
        #endif
    }

}

