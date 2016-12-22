// ---- Program Info Start----
//FileName:     RadialBasisFunction.h
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-12-03 14:38:14
// ---- Program Info End  ----

#include <vector>
#include "matrix.h"

class RadialBasisFunction{
    private:
        int _iNumSample    = 100;     //输入样本的数量
        int _iHidNode      = 10;   //隐层节点数
        int _iInDim        = 10;    //输入样本维数

        Matrix<double> _mdGreen;         //Green矩阵
        std::vector<std::vector<int>> _vviInSample;  //输入样本
        std::vector<std::vector<int>> _vviCenter;  //M个Green函数的数据中心
        std::vector<double> _vdOutReal;  //输入样本对应真实值
        Matrix<double> _mdOutReal;
        std::vector<double> _vdDelta;   //M个Green函数的扩展常数
        std::vector<double> _vdWeight;  //权值矩阵
        Matrix<double> _mdWeight;
    private:
        /*产生指定区间上均匀分布的随机数*/
        double _fnUniform(double floor,double ceil);
        /*产生区间[floor,ceil]上服从正态分布N[mu,sigma]的随机数*/
        double _fnRandomNorm(double mu,double sigma,double floor,double ceil);
        /*计算样本点距离*/
        int _fnGetDistance( const std::vector<int> &lhs, const std::vector<int> &rhs );
        /*寻找样本离哪个中心最近*/
        int _fnGetNearestCenter(const std::vector<int> &inNode);
        /*计算簇的质心*/
        std::vector<int> _fnRecalcCenter(const std::vector<std::vector<int>> &group);
        /*KMeans聚类法产生数据中心*/
        void _fnKPrototype();
        /*生成Green矩阵*/
        void _fnCalcGreen();
        /*求一个矩阵的伪逆*/
        Matrix<double> _fnGetGereralizedInverse(const Matrix<double> &matrix);
        //自组织选择法, 根据各中心之间最小距离确定扩展常数
        void _fnCalcDelta();
    public:
        RadialBasisFunction();
        RadialBasisFunction(int iNum, int iHide, int iDim, const std::vector<std::vector<int>> &vviSample, const std::vector<double> &vdReal);
        ~RadialBasisFunction();

        /*运行RBF网络*/
        void runRBF();
        /*根据网络，由输入得到输出*/
        double getEstimation( const std::vector<int> &inNode );
        /*获得样本数量*/
        int getNumberSample(){ return _iNumSample; }
        /*获得隐层数量*/
        int getNumberHidden(){ return _iHidNode; }
        /*获得样本维数*/
        int getNumberDim(){ return _iInDim; }
};
