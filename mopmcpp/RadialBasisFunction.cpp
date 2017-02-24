// ---- Program Info Start----
//FileName:     RadialBasisFunction.cpp
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-12-03 14:41:39
// ---- Program Info End  ----

#include <limits>
#include <cassert>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "RadialBasisFunction.h"
#include "matrix.h"

using namespace std;

RBF_KERNAL_TYPE RadialBasisFunction::GetKernalType
(
 const string & strKernalFunc
 ){

    if( strKernalFunc == "Cubic" )
        return k_Cubic;
    else if( strKernalFunc == "ThinPlateSpline" )
        return k_ThinPlateSpline;
    else if( strKernalFunc == "Gaussian" )
        return k_Gaussian;
    else if( strKernalFunc == "MultiQuadratic" )
        return k_MultiQuadratic;
    else if( strKernalFunc == "InverseMultiQuadratic" )
        return k_InverseMultiQuadratic;
    else
        exit( -1 );

}

RadialBasisFunction::RadialBasisFunction(){}
RadialBasisFunction::RadialBasisFunction
(
 int iNum,
 int iHide,
 int iDim,
 const std::vector< myBitSet<M> > &vviSample,
 const std::vector<double> &vdReal
 ):
    _iNumSample( iNum ),
    _iHidNode( iHide ),
    _iInDim( iDim ),
    _vviInSample( vviSample ),
    _vviCenter( iHide ),
    _vdDelta( iHide ),
    _mdOutReal( _iNumSample,1 ),
    _mdGreen( _iNumSample,_iHidNode ),
    _mdWeight( _iHidNode,1 ),
    _vdRMSE( _INT_Kernal_Num, -1 )
{
    for( int i=0; i<_iNumSample; ++i ){
        _mdOutReal.put( i, 0, vdReal[i] );
    }
}
RadialBasisFunction::~RadialBasisFunction(){}

/*产生指定区间上均匀分布的随机数*/
inline double RadialBasisFunction::_fnUniform( double floor, double ceil ){
    return floor+1.0*rand()/RAND_MAX*(ceil-floor);
}

/*产生区间[floor,ceil]上服从正态分布_iInDim[mu,sigma]的随机数*/
inline double RadialBasisFunction::_fnRandomNorm( double mu, double sigma, double floor, double ceil ){
    double x,prob,y;
    do{
        x=_fnUniform(floor,ceil);
        prob=1/sqrt(2*M_PI*sigma)*exp(-1*(x-mu)*(x-mu)/(2*sigma*sigma));
        y=1.0*rand()/RAND_MAX;
    }while(y>prob);
    return x;
}

/*根据网络，由输入得到输出*/
double RadialBasisFunction::getEstimation
(
 const myBitSet<M> & inNode,
 RBF_KERNAL_TYPE rbfKernalFunc
 )const
{
    /* return getEstimation( inNode, rbfKernalFunc ); */
    double dResult = 0.0;
    switch( rbfKernalFunc ){
        case k_Cubic:
            //cout<<"k_Cubic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*iDist);
            }
            break;
        case k_ThinPlateSpline:
            //cout<<"k_ThinPlateSpline"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*log(iDist));
            }
            break;
        case k_Gaussian:
            //cout<<"k_Gaussian"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*exp(-1.0*iDist*iDist/(2*_vdDelta[i]*_vdDelta[i]));
            }
            break;
        case k_MultiQuadratic:
            //cout<<"k_MultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*sqrt(1.0*iDist*iDist+_vdDelta[i]*_vdDelta[i]);
            }
            break;
        case k_InverseMultiQuadratic:
            //cout<<"k_InverseMultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*( 1.0 / sqrt(1.0*iDist*iDist+_vdDelta[i]*_vdDelta[i]) );
            }
            break;
        default:
            cout<<"Err"<<endl;
            exit(-1);
            break;
    }
    return dResult;
}

/*根据网络，由输入得到输出*/
double RadialBasisFunction::getEstimation
(
 const myBitSet<M> & inNode,
 RBF_KERNAL_TYPE rbfKernalFunc
 )
{
    double dResult = 0.0;
    switch( rbfKernalFunc ){
        case k_Cubic:
            //cout<<"k_Cubic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*iDist);
            }
            break;
        case k_ThinPlateSpline:
            //cout<<"k_ThinPlateSpline"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*log(iDist));
            }
            break;
        case k_Gaussian:
            //cout<<"k_Gaussian"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*exp(-1.0*iDist*iDist/(2*_vdDelta[i]*_vdDelta[i]));
            }
            break;
        case k_MultiQuadratic:
            //cout<<"k_MultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*sqrt(1.0*iDist*iDist+_vdDelta[i]*_vdDelta[i]);
            }
            break;
        case k_InverseMultiQuadratic:
            //cout<<"k_InverseMultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = GetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*( 1.0 / sqrt(1.0*iDist*iDist+_vdDelta[i]*_vdDelta[i]) );
            }
            break;
        default:
            cout<<"Err"<<endl;
            exit(-1);
            break;
    }
    return dResult;
}

double RadialBasisFunction::getEnsembleEstimation
(
  const vector<double> & vdWeights,
  const vector<RBF_KERNAL_TYPE> & vtKernalFuncs,
  const myBitSet<M> & inNode
)const
{
    assert( vdWeights.size() == vtKernalFuncs.size() );
    double fResult = 0.0;
    int iFuncNum = vdWeights.size();
    for( int i=0; i<iFuncNum; ++i){
        fResult += vdWeights[i] * getEstimation( inNode, vtKernalFuncs[i] );
    }
    return fResult;
}
double RadialBasisFunction::getEnsembleEstimation
(
  const vector<double> & vdWeights,
  const vector<RBF_KERNAL_TYPE> & vtKernalFuncs,
  const myBitSet<M> & inNode
)
{
    assert( vdWeights.size() == vtKernalFuncs.size() );
    double fResult = 0.0;
    int iFuncNum = vdWeights.size();
    for( int i=0; i<iFuncNum; ++i){
        fResult += vdWeights[i] * getEstimation( inNode, vtKernalFuncs[i] );
    }
    return fResult;
}
// 计算模型均方差
void RadialBasisFunction::_fnCalcRMSE
(
  RBF_KERNAL_TYPE rbfKernalFunc
){
    _vdRMSE[rbfKernalFunc] = 0.0;
    for( int i=0; i<_iNumSample; ++i ){
        double dErr = fabs( getEstimation(_vviInSample[i], rbfKernalFunc) - _mdOutReal.get( i, 0 ) );
        _vdRMSE[rbfKernalFunc] += dErr*dErr;
    }
    _vdRMSE[rbfKernalFunc] = sqrt( _vdRMSE[rbfKernalFunc] );
}

double RadialBasisFunction::getRMSE
(
  RBF_KERNAL_TYPE rbfKernalFunc
){
    if( _vdRMSE[rbfKernalFunc] < 0 ){
        _fnCalcRMSE( rbfKernalFunc );
    }
    return _vdRMSE[rbfKernalFunc];
}

// 计算模型在部分解上的RMSE
double RadialBasisFunction::calcRMSE
(
 RBF_KERNAL_TYPE rbfKernalFunc,
 const vector<myBitSet<M>> & vviInSample,
 const vector<double> & vdReal
)const
{
    double fRMSE = 0.0;
    int iNumSample = vviInSample.size();
    for( int i=0; i<iNumSample; ++i ){
        double dErr = fabs( getEstimation( vviInSample[i], rbfKernalFunc) - vdReal[i] );
        fRMSE += dErr * dErr;
    }
    return sqrt( fRMSE );
}

/*计算样本距离*/
int RadialBasisFunction::GetDistance
(
 const myBitSet<M> &lhs, const myBitSet<M> &rhs
 ){
    /* assert( lhs.size() == rhs.size() ); */
    /* int iDim = lhs.size(); */
    /* int iDist = 0.0; */
    /* for( int i=0; i<iDim; ++i ){ */
    /*     iDist += abs(lhs[i] - rhs[i]); */
    /* } */
    myBitSet<M> bitLHS( lhs );
    bitLHS ^= rhs;
    return bitLHS.SSE4_count();
}

/*寻找样本离哪个中心最近*/
int RadialBasisFunction::_fnGetNearestCenter(const myBitSet<M> &inNode){
    int iIdx=-1;
    int iMinDist=numeric_limits<int>::max();
    for(int i=0;i<_vviCenter.size();++i){
        int iDist = GetDistance( inNode, _vviCenter[i] );
        if(iDist<iMinDist){
            iMinDist=iDist;
            iIdx=i;
        }
    }
    return iIdx;
}

/*计算簇的质心*/
myBitSet<M> RadialBasisFunction::_fnRecalcCenter(const vector< myBitSet<M> > &vviGroup){
    int iGroupSize=vviGroup.size();
    myBitSet<M> viCenter;
    for( int j=0; j<_iInDim; ++j ){
        int iCnt = 0;
        for(int i=0; i<iGroupSize; ++i) {
            iCnt += vviGroup[i].test(j);
        }
        if( iCnt > iGroupSize-iCnt )
            viCenter.set(j);
    }
    return viCenter;
}

/*KMeans聚类法产生数据中心*/
void RadialBasisFunction::_fnKPrototype(){
    assert(_iNumSample % _iHidNode == 0);
    vector<vector<myBitSet<M>>> vvviPackage(_iHidNode);          //记录各个聚类中包含哪些样本
    int iCnt = 0;
    int iSep  = _iNumSample / _iHidNode;
    int iStep = rand() % iSep;
    for(int i=0;i<_iHidNode;++i){   //从_iNumSample个输入样本中随机选_iHidNode个作为初始聚类中心
        _vviCenter[i] = _vviInSample[iSep*i+iStep];
    }
    while(1){
        for(int i=0;i<_iHidNode;++i)
            vvviPackage[i].clear();   //先清空聚类信息
        for(int i=0;i<_iNumSample;++i){       //把所有输入样本归到对应的簇
            int c=_fnGetNearestCenter(_vviInSample[i]);
            vvviPackage[c].push_back(_vviInSample[i]);
        }
        vector<myBitSet<M> > vviNextCenter(_iHidNode);       //存储新的簇心
        for(int i=0;i<_iHidNode;++i){
            vector<myBitSet<M> > vviGroup = vvviPackage[i];
            vviNextCenter[i]=_fnRecalcCenter(vviGroup);
        }
        bool bFlag=false;
        for(int i=0;i<_iHidNode;++i){       //检查前后两次质心的改变量是否都小于gap
            int iDist = GetDistance( vviNextCenter[i], _vviCenter[i] );
            if(iDist != 0){
                bFlag=true;
                break;
            }
        }
        _vviCenter=vviNextCenter;
        iCnt++;
        if(!bFlag)
            break;
    }
}

/*生成Green矩阵*/
void RadialBasisFunction::_fnCalcGreen(){
    for(int i=0;i<_iNumSample;++i){
        for(int j=0;j<_iHidNode;++j){
            int iDist = GetDistance(_vviInSample[i], _vviCenter[j]);
            _mdGreen.put(i,j,exp(-1.0*(iDist)*(iDist)/(2*_vdDelta[j]*_vdDelta[j])));
        }
    }
}

/*求一个矩阵的伪逆*/
Matrix<double> RadialBasisFunction::_fnGetGereralizedInverse(const Matrix<double> &matrix){
    return (matrix.getTranspose()*matrix).getInverse()*(matrix.getTranspose());
}

//自组织选择法, 根据各中心之间最小距离确定扩展常数
void RadialBasisFunction::_fnCalcDelta(){
    for( int i=0; i<_iHidNode; ++i ){
        _vdDelta[i] = numeric_limits<double>::max();
        for( int j=0; j<_iHidNode; ++j ){
            int iDist = GetDistance( _vviCenter[i], _vviCenter[j] );
            if( iDist != 0 ){ // 判断iDist是否为0,相同的点距离为0
                _vdDelta[i] = _vdDelta[i] <= iDist ? _vdDelta[i] : iDist;
            }
        }
    }

}
/*运行RBF网络*/
void RadialBasisFunction::runRBF(){
    //K-protytype 选取中心点
    _fnKPrototype();
    //根据center计算delta
    _fnCalcDelta();
    //计算Green矩阵
    _fnCalcGreen();     //计算Green矩阵
    //计算权值weight
    _mdWeight=_fnGetGereralizedInverse(_mdGreen)*_mdOutReal;      //计算权值矩阵

}
