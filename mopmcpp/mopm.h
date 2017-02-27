// ---- Program Info Start----
//FileName:     MOPM.h
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-01-03 10:54:58
// ---- Program Info End  ----

#ifndef _MOPM_H_
#define _MOPM_H_

#include <vector>
#include "myBitSet.h"
#include "RadialBasisFunction.h"

#ifndef N
#define N 573652
#endif // define N

#ifndef M
#define M 273
#endif //define M

#define MAX_VALUE 10000000

template <typename T>
std::vector<size_t> fnSortIndex(const std::vector<T> & v, int iType);

typedef struct IndividualNode
{
    myBitSet<M> _bitTransaction;
    std::vector<double> _vfFitness;
    //std::vector<double> _vfEstimates;
    //std::vector<int> _viTransaction;
    int _iFrontNo;
    double _fCrowdDistance;
    bool _bIsDuplicate;

    IndividualNode()
    {
        _iFrontNo = 0;
        _fCrowdDistance = 0;
        _bIsDuplicate = false;
    };
}IndividualNode;

class NSGAII
{
    int _iPopSize, _iPopDims, _iDataSize, _iDataDims;
    int _iObjDims;
    std::vector<myBitSet<M> > _vbitDatabase;
    std::vector<myBitSet<N> > _vbitItemDatabase;
    std::vector<int> _viTransLength;
    std::vector<int> _viItemLength;
    std::vector<size_t> _vsizetSortedItemLengthIndex;

    public:
    NSGAII();
    ~NSGAII();

    NSGAII(
            int iPopSize,
            int iPopDims,
            int iDataSize,
            int iDataDims,
            int iObjDims,
            const std::vector<myBitSet<M>> &vbitDatabase
          );
    std::vector<IndividualNode> _fnInitialization(

            );
    std::vector<RadialBasisFunction> _fnBuildModel(
            const std::vector<IndividualNode> &vnodeDatabase,
            std::vector<clock_t> &
            );
    void _fnCalcEstimation(
            const std::vector<std::vector<RadialBasisFunction>> & ,
            std::vector<IndividualNode> &,
            int
            );
    void _fnCalcEstimation(
            const std::vector<RadialBasisFunction> & vmRBFModels,
            std::vector<IndividualNode> & vnodePopulations
            );
    void _fnCalcFiteness(
            std::vector<IndividualNode> &vnodePopulations
            );
    void _fnCalcFiteness(
            IndividualNode &vnodePopulations
            );
    std::vector<std::vector<int> > _fnNonDominateSort(
            std::vector<IndividualNode> &vnodePopulations
            );
    void _fnCalcCrowdDistance(
            std::vector<IndividualNode> & vnodePopulations,
            const std::vector<std::vector<int> > &vviFrontList
            );
    std::vector<IndividualNode> _fnSelectMatingPool(
            const std::vector<IndividualNode> &vnodePopulations
            );
    void _fnReproduceOff(
            std::vector<IndividualNode> &vnodeMatePool
            );
    std::vector<IndividualNode> _fnNatureSelection(
            const std::vector<IndividualNode> &vnodeAllPop,
            const std::vector<std::vector<int> > & vviFrontList
            );
    std::vector<IndividualNode> _fnNatureSelectionNoDuplicate(
            const std::vector<IndividualNode> &vnodeAllPop,
            const std::vector<std::vector<int> > & vviFrontList
            );
    std::vector<IndividualNode> _fnMOEC(
            int &,
            double &
            );
    void _fnSMEC(
            const std::vector<std::vector<RadialBasisFunction>> & ,
            std::vector<IndividualNode> & ,
            int
            );
    void _fnSMEC(
            const std::vector<RadialBasisFunction> & vmRBFModels,
            std::vector<IndividualNode> & vnodePopulations
            );
    bool _fnCheckSimilar(
            std::vector<IndividualNode> lhs,
            std::vector<IndividualNode> rhs
            );
    double _fnVerifyAccuracy(
            const std::vector<RadialBasisFunction> & vmRBFModels,
            const std::vector<IndividualNode> & vnodeDatabase
            );
    void _fnDebugPrintInfo(
            std::ofstream & logs,
            std::vector<IndividualNode> & vnodePopulations
            );
    void _fnLocalSearchPhase(
            std::vector<IndividualNode> &
            /* const std::vector<IndividualNode> & */
            );
    void _fnLocalSearchPhase(
            std::vector<IndividualNode> &,
            const std::vector<IndividualNode> &
            );
    IndividualNode _fnFindLocalOptima(
            IndividualNode,
            const std::vector<double> &
            );
    std::vector<int> _fnFindLocalOptima(
            std::vector<int> ,
            const RadialBasisFunction &,
            const std::vector<double> &
            );
    double _fnGetEnsembleEstimaion(
            const std::vector<int> &,
            const RadialBasisFunction &,
            const std::vector<double> &
            );
};

#endif //_MOPM_H_
