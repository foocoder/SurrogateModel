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
const int N = 15;
const int P = 200;
const int M = 20;

void readData( string filename ){
    ifstream ifile( filename.c_str() );
    if( !ifile ){
        cerr<<"Open file "<<filename<<" failed!"<<endl;
        exit(-1);
    }
    string strLine;
    while( getline( ifile, strLine ) ){
        istringstream issLine( strLine );
        vector<int> pattern(N);
        for( int i=0; i<N; ++i ) {
            issLine >> pattern[i];
        }
        char delimiter;
        issLine >> delimiter;
        vector<double> obj(2);
        issLine>>obj[1]>>delimiter>>obj[0];

        DB.insert( Ind::value_type( pattern, obj ) );
    }
}

vector<int> Binary2Gray( int bin ){
    stack<int> grayStack;
    vector<int> grayCode( N );
    for( int i=0; i<N; ++i ){
        if( bin & 1 )
            grayStack.push( 1 );
        else
            grayStack.push( 0 );
        bin >>= 1;
    }
    for( int i=0; i<N; ++i ){
        grayCode[i] = grayStack.top();
        grayStack.pop();
    }
    return grayCode;
}

int main(int argc, char *argv[])
{

    /*生成随机种子*/
    srand(time(0));
    string filename="./PatternDB";
    //cout<<"Read Data"<<endl;
    readData(filename.c_str());
    int DBSize = pow(2, N);
    vector<int> trainIDX(P);
    vector<int> testIDX(P);
    int seed = rand() % DBSize;
    //cout<<"Generate Traing Data"<<endl;
    for( int i=0; i<P; ++i ){
        trainIDX[i] = seed;
        testIDX[i] = (seed+1) % DBSize;
        seed = (seed+5) % DBSize;
    }

    vector<vector<int>> trainData(P);
    vector<vector<int>> testData(P);
    vector<double> trainReal(P);
    for( int i=0; i<P; ++i ){
        trainData[i] = Binary2Gray( trainIDX[i] );
        testData[i]  = Binary2Gray( testIDX[i] );
        trainReal[i] = DB[trainData[i]][0];
    }

    //cout<<"Training"<<endl;
    RadialBasisFunction rbf( P, M, N, trainData, trainReal );
    rbf.runRBF();
    //cout<<"Testing"<<endl;
    for( int i=0; i<P; ++i ){
        cout<<setw(12)<<rbf.getEstimation( testData[i] )<<"\t"<<setw(12)<<DB[testData[i]][0]<<"\t"<<setw(12)<<fabs(DB[testData[i]][0]-rbf.getEstimation(testData[i]))<<endl;
    }

    return 0;
}
