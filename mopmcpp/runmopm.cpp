// ---- Program Info Start----
//FileName:     runMOPM.cpp
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-01-04 20:58:51
// ---- Program Info End  ----

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include "mopm.h"

using namespace std;

void fnReadData(string strFileName, vector<myBitSet<M> > &vbitDatabase, int & iRowCounter, int & iColumnCounter){
    ifstream infile(strFileName.c_str(),ios::in|ios::binary);
    string strLine;
    int iWord;
    if(!infile){
        cerr << "error not opened!" << endl;
        return;
    }
    //cout<<"Start Read Data"<<endl;
    while(getline(infile,strLine))
    {
        istringstream stream(strLine);
        iColumnCounter = 0;
        myBitSet<M> bitTransTemp;
        while(stream >> iWord)
        {
            if(iWord == 1){
                bitTransTemp.set(iColumnCounter);
            }
            else if(iWord == 0){
            }
            iColumnCounter ++;
        }
        vbitDatabase.push_back(bitTransTemp);
        iRowCounter++;
    }
}

void fnRunNSGAII(const string & strDataPath, const string & strResPath, const string& inputFile)
{
    vector<myBitSet<M>> vbitDatabase;
    vbitDatabase.reserve(N);
    int iRowNum=0, iColumnNum=0;
    string input = strDataPath + "/" + inputFile;
    fnReadData(input, vbitDatabase, iRowNum, iColumnNum);

    int iPopSize = (iColumnNum / 50 + 1) * 50;
    //iPopSize = iPopSize > 1000 ? 1000 : iPopSize;
    //iPopSize = 20;
    cout<<"Starting NSGAII..."<<endl;
    NSGAII algorithm = NSGAII(iPopSize, iColumnNum, iRowNum, iColumnNum, 3, vbitDatabase);
    int traversNode;
    double spendTime;
    const vector<IndividualNode> & vnodeResult = algorithm._fnMOEC(traversNode, spendTime);
    //strDataPath += ".results";
    string resFile = strResPath + "/" + inputFile + ".results";
    ofstream resultFiles( resFile.c_str() );

    for(auto nodeIter:vnodeResult)
    {
        for(int i=0; i<iColumnNum; i++)
        {
            if(nodeIter._bitTransaction.test(i))
            {
                resultFiles<<i<<",";
            }
        }
        resultFiles<<"\b;";
        for(auto i : nodeIter._vfFitness)
        {
            resultFiles<<i<<",";
        }
        resultFiles<<"\b "<<endl;
    }
    resultFiles<<spendTime<<","<<traversNode<<endl;
}

int main(int argc, char *argv[])
{
    if( argc < 4 ){
        cerr<<"Usage: ./runmopm datapath resultpath inputfile"<<endl;
        exit(-1);
    }
    fnRunNSGAII(argv[1], argv[2], argv[3]);
    return 0;
}
