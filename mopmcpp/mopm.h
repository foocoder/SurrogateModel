// ---- Program Info Start----
//FileName:     MOPM.h
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-01-03 10:54:58
// ---- Program Info End  ----

#include <vector>
#include "myBitSet.h"
//#include <bitset>

//#define NBitVector bitset<N>
//#define MBitVector bitset<M>
//#define N 573652
//#define M 273
#define MAX_VALUE 10000000
using namespace std;

template <typename T>
vector<size_t> fnSortIndex(const vector<T> & v, int iType);

typedef struct IndividualNode
{
    myBitSet<M> _bitTransaction;
    vector<double> _vfFitness;
    int _iFrontNo;
    double _fCrowdDistance;
    bool _bIsDuplicate;

    IndividualNode()
    {
        _bitTransaction = myBitSet<M>();
        _iFrontNo = 0;
        _fCrowdDistance = 0;
        _bIsDuplicate = false;
    };
}IndividualNode;

class NSGAII
{
    int _iPopSize, _iPopDims, _iDataSize, _iDataDims;
    int _iObjDims;
    vector<myBitSet<M> > _vbitDatabase;
    vector<myBitSet<N> > _vbitItemDatabase;
    vector<int> _viTransLength;
    vector<int> _viItemLength;
    vector<size_t> _vsizetSortedItemLengthIndex;

public:
    NSGAII();
    ~NSGAII();

    NSGAII(int iPopSize, int iPopDims, int iDataSize, int iDataDims, int iObjDims, const vector<myBitSet<M>> &vbitDatabase);
    vector<IndividualNode> _fnInitialization();
    void _fnCalcFiteness(vector<IndividualNode> &vnodePopulations);
    vector<vector<int> > _fnNonDominateSort(vector<IndividualNode> &vnodePopulations);
    void _fnCalcCrowdDistance( vector<IndividualNode> & vnodePopulations, const vector<vector<int> > &vviFrontList );
    vector<IndividualNode> _fnSelectMatingPool(const vector<IndividualNode> &vnodePopulations);
    void _fnReproduceOff(vector<IndividualNode> &vnodeMatePool);
    vector<IndividualNode> _fnNatureSelection(const vector<IndividualNode> &vnodeAllPop, const vector<vector<int> > & vviFrontList);
    vector<IndividualNode> _fnNatureSelectionNoDuplicate(const vector<IndividualNode> &vnodeAllPop, const vector<vector<int> > & vviFrontList);
    vector<IndividualNode> _fnMOEC( int &, double & );
    bool _fnCheckSimilar( vector<IndividualNode> lhs, vector<IndividualNode> rhs );
    void _fnDebugPrintInfo( ofstream & logs, vector<IndividualNode> & vnodePopulations );
};
