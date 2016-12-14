// ---- Program Info Start----
//FileName:     CalclulatePattern.cpp
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-06-06 18:44:30
// ---- Program Info End  ----

#include <iostream>
#include <iomanip>
#include <fstream>
#include "CalclulatePattern.h"

using namespace std;

int main(int argc, char *argv[])
{
    const char * filename;
    if(argc <= 1)
        filename = "./data";
    else
        filename = argv[1];
    CReadFiles readData(filename);
    readData.readBitSet();
    readData.displayBitset();
    CCalculatePattern calcPattern(readData.getRowTrans(), readData.getColumnTrans());

    CGrayCode gray;

    vector<myBitSet<N>> grayCode = gray.getGrayCode();
    ofstream ofile( "./Data/PatternDB" );
    if( !ofile ){
        cerr<<"Open file failed!"<<endl;
        exit(-1);
    }
    for(auto pattern : grayCode){
        vector<double> res = calcPattern.calculate(pattern);
        cout<<pattern<<":"<<setw(12)<<pattern.SSE4_count()<<","<<setw(12)<<res[0]<<",\t"<<setw(12)<<res[1]<<endl;
        ofile<<pattern<<":"<<setw(12)<<res[0]<<",\t"<<setw(12)<<res[1]<<endl;
    }

    return 0;
}
