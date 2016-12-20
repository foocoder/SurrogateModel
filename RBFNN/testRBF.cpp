// ---- Program Info Start----
//FileName:     testRBF.cpp
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-12-13 16:55:12
// ---- Program Info End  ----

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <stack>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "RadialBasisFunction.h"

using namespace std;

typedef map<vector<int>, vector<double>> Ind;

Ind DB;  //数据库
const int iInDim        = 15;
const int iNumTrains    = 200;
const int iNumHideNodes = 20;
//const int iObjType      = 1;

void readData( string filename ){
    ifstream ifile( filename.c_str( ) );
    if( !ifile ){
        cerr<<"Open file "<<filename<<" failed!"<<endl;
        exit( -1 );
    }
    string strLine;
    while( getline( ifile, strLine ) ){
        istringstream issLine( strLine );
        vector<int> pattern( iInDim );
        for( int i=0; i<iInDim; ++i ) {
            issLine >> pattern[i];
        }
        char delimiter;
        issLine >> delimiter;
        vector<double> obj( 2 );
        issLine>>obj[0]>>delimiter>>obj[1];

        DB.insert( Ind::value_type( pattern, obj ) );
    }
}

vector<int> Binary2Gray( int bin ){
    stack<int> grayStack;
    vector<int> grayCode( iInDim );
    for( int i=0; i<iInDim; ++i ){
        if( bin & 1 )
            grayStack.push( 1 );
        else
            grayStack.push( 0 );
        bin >>= 1;
    }
    for( int i=0; i<iInDim; ++i ){
        grayCode[i] = grayStack.top( );
        grayStack.pop( );
    }
    return grayCode;
}

int main( int argc, char *argv[] )
{

    /*生成随机种子*/
    srand( time( NULL ) );
    string filename="./data/PatternDB.db";

    readData( filename.c_str( ) );
    int DBSize = pow( 2, iInDim );
    int iRunTotal = 5000;

    for( int iObjType = 0; iObjType<2; ++iObjType )
    {
        //多次运行, 寻找合适步长
        for( int iStep = 1; iStep < 20; ++iStep )
        {
            vector<double> vfErrors(iRunTotal);
            for( int iRun = 0; iRun<iRunTotal; ++iRun )
            {
                vector<int> trainIDX( iNumTrains );
                vector<int> testIDX( iNumTrains );
                int seed = rand( ) % DBSize;
                for( int i=0; i<iNumTrains; ++i ){
                    trainIDX[i] = seed;
                    //testIDX[i] = ( seed+1 ) % DBSize;
                    testIDX[i] = rand() % DBSize;
                    seed = ( seed+iStep ) % DBSize;
                }

                vector<vector<int>> trainData( iNumTrains );
                vector<vector<int>> testData( iNumTrains );
                vector<double> trainReal( iNumTrains );
                for( int i=0; i<iNumTrains; ++i ){
                    trainData[i] = Binary2Gray( trainIDX[i] );
                    testData[i]  = Binary2Gray( testIDX[i] );
                    trainReal[i] = DB[trainData[i]][iObjType];
                }

                RadialBasisFunction rbf( iNumTrains, iNumHideNodes, iInDim, trainData, trainReal );
                rbf.runRBF( );

                double fTotalErr = 0.0;
                for( int i=0; i<iNumTrains; ++i ){
                    //cout<<setw( 12 )<<rbf.getEstimation( testData[i] )<<"\t"<<setw( 12 )<<DB[testData[i]][iObjType]<<"\t"<<setw( 12 )<<fabs( DB[testData[i]][iObjType]-rbf.getEstimation( testData[i] ) )<<endl;
                    fTotalErr += fabs( rbf.getEstimation( testData[i] ) - DB[testData[i]][iObjType] );
                }
                //cout<<"Step: "<<setw(3)<<seed<<"\tTotal Error: "<<setw(12)<<fTotalErr<<endl;
                vfErrors[iRun] = fTotalErr;
            }

            double fAvgErr = 0.0;
            double fStdErr = 0.0;
            for( double i : vfErrors ){
                fAvgErr += i;
            }
            fAvgErr /= vfErrors.size();
            for( double i : vfErrors ){
                fStdErr += fabs( i-fAvgErr ) * fabs( i-fAvgErr );
            }
            fStdErr /= vfErrors.size();
            fStdErr = sqrt( fStdErr );
            cout<<setw(10)<<(iObjType == 0 ? "Support" : "Occupancy")<<"\tStep: "<<setw(4)<<iStep<<"\tAvg Err: "<<setw(12)<<fAvgErr<<"\tStd."<<setw(12)<<fStdErr<<endl;
        }

    }
    return 0;
}
