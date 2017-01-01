// ---- Program Info Start----
//FileName:     RadialBasisFunction.cpp
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-12-03 14:41:39
// ---- Program Info End  ----

#include <cassert>
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

RadialBasisFunction::RadialBasisFunction(){}
RadialBasisFunction::RadialBasisFunction
(
 int iNum,
 int iHide,
 int iDim,
 const std::vector<std::vector<int>> &vviSample,
 const std::vector<double> &vdReal,
 const std::string &strKernalFun
 ):
    _iNumSample( iNum ),
    _iHidNode( iHide ),
    _iInDim( iDim ),
    _vviInSample( vviSample ),
    _vviCenter( iHide ),
    _vdDelta( iHide ),
    _mdOutReal( _iNumSample,1 ),
    _mdGreen( _iNumSample,_iHidNode ),
    _mdWeight( _iHidNode,1 )
{
    for( int i=0; i<_iNumSample; ++i ){
        _mdOutReal.put( i, 0, vdReal[i] );
    }
    if( strKernalFun == "Cubic" )
        _rbfKernalFunc = k_Cubic;
    if( strKernalFun == "ThinPlateSpline" )
        _rbfKernalFunc = k_ThinPlateSpline;
    if( strKernalFun == "Gaussian" )
        _rbfKernalFunc = k_Gaussian;
    if( strKernalFun == "MultiQuadratic" )
        _rbfKernalFunc = k_MultiQuadratic;
    if( strKernalFun == "InverseMultiQuadratic" )
        _rbfKernalFunc = k_InverseMultiQuadratic;
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
 const vector<int> & inNode
 ){
    double dResult = 0.0;
    switch( _rbfKernalFunc ){
        case k_Cubic:
            //cout<<"k_Cubic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*iDist);
            }
            break;
        case k_ThinPlateSpline:
            //cout<<"k_ThinPlateSpline"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*log(iDist));
            }
            break;
        case k_Gaussian:
            //cout<<"k_Gaussian"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*exp(-1.0*iDist*iDist/(2*_vdDelta[i]*_vdDelta[i]));
            }
            break;
        case k_MultiQuadratic:
            //cout<<"k_MultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*sqrt(1.0*iDist*iDist+_vdDelta[i]*_vdDelta[i]);
            }
            break;
        case k_InverseMultiQuadratic:
            //cout<<"k_InverseMultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
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
 const vector<int> & inNode
 ) const
{
    double dResult = 0.0;
    switch( _rbfKernalFunc ){
        case k_Cubic:
            //cout<<"k_Cubic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*iDist);
            }
            break;
        case k_ThinPlateSpline:
            //cout<<"k_ThinPlateSpline"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*(iDist*iDist*log(iDist));
            }
            break;
        case k_Gaussian:
            //cout<<"k_Gaussian"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*exp(-1.0*iDist*iDist/(2*_vdDelta[i]*_vdDelta[i]));
            }
            break;
        case k_MultiQuadratic:
            //cout<<"k_MultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
                dResult += _mdWeight.get(i,0)*sqrt(1.0*iDist*iDist+_vdDelta[i]*_vdDelta[i]);
            }
            break;
        case k_InverseMultiQuadratic:
            //cout<<"k_InverseMultiQuadratic"<<endl;
            for( int i=0; i<_iHidNode; ++i ){
                int iDist = _fnGetDistance( inNode, _vviCenter[i] );
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

// 计算模型均方差
void RadialBasisFunction::_fnCalcRMSE(){
    _dRMSE = 0.0;
    for( int i=0; i<_iNumSample; ++i ){
        double dErr = fabs( getEstimation(_vviInSample[i]) - _vdOutReal[i] );
        _dRMSE += dErr*dErr;
    }
    _dRMSE = sqrt( _dRMSE );
}

inline double RadialBasisFunction::getRMSE(){
    if( _dRMSE < 0 ){
        _fnCalcRMSE();
    }
    return _dRMSE;
}

/*计算样本距离*/
int RadialBasisFunction::_fnGetDistance( const vector<int> &lhs, const vector<int> &rhs ){
    int iDist = 0.0;
    for( int i=0; i<_iInDim; ++i ){
        iDist += abs(lhs[i] - rhs[i]);
    }
    return iDist;
}

/*计算样本距离*/
int RadialBasisFunction::_fnGetDistance( const vector<int> &lhs, const vector<int> &rhs )const{
    int iDist = 0.0;
    for( int i=0; i<_iInDim; ++i ){
        iDist += abs(lhs[i] - rhs[i]);
    }
    return iDist;
}
/*寻找样本离哪个中心最近*/
int RadialBasisFunction::_fnGetNearestCenter(const vector<int> &inNode){
    int iIdx=-1;
    int iMinDist=numeric_limits<int>::max();
    for(int i=0;i<_vviCenter.size();++i){
        int iDist = _fnGetDistance( inNode, _vviCenter[i] );
        if(iDist<iMinDist){
            iMinDist=iDist;
            iIdx=i;
        }
    }
    return iIdx;
}

/*计算簇的质心*/
vector<int> RadialBasisFunction::_fnRecalcCenter(const vector<vector<int>> &vviGroup){
    int iGroupSize=vviGroup.size();
    vector<int> viCenter(_iInDim);
    for( int j=0; j<_iInDim; ++j ){
        int iCnt = 0;
        for(int i=0; i<iGroupSize; ++i) {
            iCnt += vviGroup[i][j];
        }
        if( iCnt > iGroupSize-iCnt )
            viCenter[j] = 1;
        else
            viCenter[j] = 0;
    }
    return viCenter;
}

/*KMeans聚类法产生数据中心*/
void RadialBasisFunction::_fnKPrototype(){
    assert(_iNumSample % _iHidNode == 0);
    vector<vector<vector<int>>> vvviPackage(_iHidNode);          //记录各个聚类中包含哪些样本
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
        vector<vector<int> > vviNextCenter(_iHidNode);       //存储新的簇心
        for(int i=0;i<_iHidNode;++i){
            vector<vector<int> > vviGroup=vvviPackage[i];
            vviNextCenter[i]=_fnRecalcCenter(vviGroup);
        }
        bool bFlag=false;
        for(int i=0;i<_iHidNode;++i){       //检查前后两次质心的改变量是否都小于gap
            int iDist = _fnGetDistance( vviNextCenter[i], _vviCenter[i] );
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
            int iDist = _fnGetDistance(_vviInSample[i], _vviCenter[j]);
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
            int iDist = _fnGetDistance( _vviCenter[i], _vviCenter[j] );
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
