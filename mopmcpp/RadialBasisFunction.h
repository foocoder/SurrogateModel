// ---- Program Info Start----
//FileName:     RadialBasisFunction.h
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-12-03 14:38:14
// ---- Program Info End  ----

#ifndef _RBF_H_
#define _RBF_H_

#include <vector>
#include <cstring>
#include "matrix.h"
#include "myBitSet.h"

#ifndef N
#define N 573652
#endif // define N

#ifndef M
#define M 273
#endif //define M

typedef enum KERNAL {
    k_Cubic = 0,
    k_ThinPlateSpline,
    k_Gaussian,
    k_MultiQuadratic,
    k_InverseMultiQuadratic
} RBF_KERNAL_TYPE;

class RadialBasisFunction{
    private:
        int _iNumSample    = 100;     //输入样本的数量
        int _iHidNode      = 10;   //隐层节点数
        int _iInDim        = 10;    //输入样本维数
        static const int _INT_Kernal_Num = 5;

        /* double _dRMSE      = -1.0;  // 均方差 */
        /* RBF_KERNAL_TYPE _rbfKernalFunc = k_Gaussian; // 核函数类型. */

        std::vector<double> _vdRMSE;  // 不同核函数的均方差
        Matrix<double> _mdGreen;         //Green矩阵
        std::vector< myBitSet<M> > _vviInSample;  //输入样本
        std::vector< myBitSet<M> > _vviCenter;  //M个Green函数的数据中心
        //std::vector<double> _vdOutReal;  //输入样本对应真实值
        Matrix<double> _mdOutReal;
        std::vector<double> _vdDelta;   //M个Green函数的扩展常数
        std::vector<double> _vdWeight;  //权值矩阵
        Matrix<double> _mdWeight;
    private:
        // 计算model均方差
        void _fnCalcRMSE(
            RBF_KERNAL_TYPE
                );
        /*产生指定区间上均匀分布的随机数*/
        double _fnUniform(
            double floor,
            double ceil
                );
        /*产生区间[floor, ceil]上服从正态分布N[mu, sigma]的随机数*/
        double _fnRandomNorm(
            double mu,
            double sigma,
            double floor,
            double ceil
                );
        /*寻找样本离哪个中心最近*/
        int _fnGetNearestCenter(
            const myBitSet<M> &inNode
                );
        /*计算簇的质心*/
        myBitSet<M> _fnRecalcCenter(
            const std::vector<myBitSet<M>> &group
                );
        /*KMeans聚类法产生数据中心*/
        void _fnKPrototype(

                );
        /*生成Green矩阵*/
        void _fnCalcGreen(

                );
        /*求一个矩阵的伪逆*/
        Matrix<double> _fnGetGereralizedInverse(
            const Matrix<double> &matrix
                );
        //自组织选择法 根据各中心之间最小距离确定扩展常数
        void _fnCalcDelta(

                );
    public:
        RadialBasisFunction(

                );
        RadialBasisFunction(
            int iNum,
            int iHide,
            int iDim,
            const std::vector< myBitSet<M> > &vviSample,
            const std::vector<double> &vdReal
                );
        ~RadialBasisFunction(

                );

        /*运行RBF网络*/
        void runRBF(

                );
        /*根据网络，由输入得到输出*/
        double getEstimation(
            const myBitSet<M> &inNode,
            RBF_KERNAL_TYPE
                );
        /*根据网络，由输入得到输出*/
        double getEstimation(
            const myBitSet<M> &inNode,
            RBF_KERNAL_TYPE
                ) const;
        double getEnsembleEstimation (
            const std::vector<double> & vdWeights,
            const std::vector<RBF_KERNAL_TYPE> & vtKernalFuncs,
            const myBitSet<M> & inNode
                )const;
        double getEnsembleEstimation (
            const std::vector<double> & vdWeights,
            const std::vector<RBF_KERNAL_TYPE> & vtKernalFuncs,
            const myBitSet<M> & inNode
                );
        /*获得样本数量*/
        int getNumberSample(

                ){ return _iNumSample; }
        /*获得隐层数量*/
        int getNumberHidden(

                ){ return _iHidNode; }
        /*获得样本维数*/
        int getNumberDim(

                ){ return _iInDim; }

        double getRMSE(
            RBF_KERNAL_TYPE rbfKernalFunc
                );
        double calcRMSE(
            RBF_KERNAL_TYPE rbfKernalFunc,
            const std::vector<myBitSet<M>> & vviInSample,
            const std::vector<double> & vdReal
                )const;
        static RBF_KERNAL_TYPE GetKernalType(
            const std::string &
                );
        /*计算样本点距离*/
        static int GetDistance(
            const myBitSet<M> &lhs,
            const myBitSet<M> &rhs
                );
};

#endif // _RBF_H_
