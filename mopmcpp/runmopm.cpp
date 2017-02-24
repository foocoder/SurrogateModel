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
#include <iomanip>
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

void fnRunNSGAII
(
 const string & strDataPath,
 const string & strResPath,
 const string & strInDBName,
 const string & strSaveName
 )
{
    vector<myBitSet<M>> vbitDatabase;
    vbitDatabase.reserve(N);
    int iRowNum=0, iColumnNum=0;
    string input = strDataPath + "/" + strInDBName;
    fnReadData(input, vbitDatabase, iRowNum, iColumnNum);

    int iPopSize = (iColumnNum / 50 + 1) * 50;
    iPopSize = iPopSize > 100 ? 100 : iPopSize;

    //cout<<"Starting NSGAII..."<<endl;
    NSGAII algorithm = NSGAII(iPopSize, iColumnNum, iRowNum, iColumnNum, 3, vbitDatabase);
    int traversNode;
    double spendTime;
    const vector<IndividualNode> & vnodeResult = algorithm._fnMOEC(traversNode, spendTime);
    string strResPatName = strResPath + "/" + strSaveName + ".pat";
    string strResFitName = strResPath + "/" + strSaveName + ".fit";
    ofstream ofsResPatFile( strResPatName.c_str() );
    ofstream ofsResFitFile( strResFitName.c_str() );

    ofsResFitFile<<"#"<<endl;
    for(auto nodeIter:vnodeResult)
    {
        for(int i=0; i<iColumnNum; i++)
        {
            if(nodeIter._bitTransaction.test(i))
            {
                ofsResPatFile<<i<<" ";
            }
        }
        ofsResPatFile<<endl;
        for(auto i : nodeIter._vfFitness)
        {
            ofsResFitFile<<setw(12)<<i<<"\t";
        }
        ofsResFitFile<<endl;
    }
    ofsResFitFile<<"#"<<endl;
    cout<<"File:"<<setw(5)<<strSaveName
        <<"\tTime:"<<setw(8)<<setprecision(10)<<spendTime
        <<"\tNode:"<<setw(8)<<traversNode<<endl;
}

int main(int argc, char *argv[])
{
    if( argc < 5 ){
        cerr<<"Usage: ./a.out datapath resultpath strInDBName strSaveFiles"<<endl;
        exit(-1);
    }
    fnRunNSGAII(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}
