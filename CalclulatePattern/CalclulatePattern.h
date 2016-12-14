// ---- Program Info Start----
//FileName:     CalclulatePattern.h
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-06-06 19:19:19
// ---- Program Info End  ----

#ifndef _CALCPATTERN_H_
#define _CALCPATTERN_H_

#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>
#include "myBitSet.h"

//#define N 5  // The rows of data
//#define M 8  // The number of items of pattern

class CReadFiles{
    public:
        CReadFiles(const char * filename):_fBinPattern(filename, std::ios::in|std::ios::binary),
            _rowTrans(M), _columnTrans(N){
            if(!_fBinPattern){
                std::cerr<<"Error: file "<<filename<<" Opened Failed!"<<std::endl;
                return ;
            }
        }
        ~CReadFiles(){
            if(!_fBinPattern)
                return;
            _fBinPattern.close();
        }
        CReadFiles(const CReadFiles &) = delete;
        CReadFiles & operator=(const CReadFiles &) = delete;

        void readBitSet(){
            if(!_fBinPattern){
                return;
            }
            std::string strline;
            int row = 0, column;
            while(std::getline(_fBinPattern, strline))
            {
                std::istringstream streamline( strline );
                int word;
                column = N-1;
                while( streamline>>word ){
                    if( word==1 ){
                        _rowTrans[row].set(column);
                        _columnTrans[N-column-1].set(row);
                        column--;
                    }
                    else if(word == 0){
                        column--;
                    }
                }
                row++;
            }
            //_row    = row;
            //_column = column;
        }
        void displayBitset(){
            if(_rowTrans.empty()){
                std::cerr<<"Error: transaction is empty"<<std::endl;
                return;
            }
            std::cout<<"Display row transactions"<<std::endl;
            for(auto record : _rowTrans){
                std::cout<<record<<std::endl;
            }
            std::cout<<"Display column transactions"<<std::endl;
            for(auto record : _columnTrans){
                std::cout<<record<<std::endl;;
            }
            std::cout<<"Total: "<<M<<" rows, "<<N<<" columns"<<std::endl;
        }
        std::vector<myBitSet<N>> getRowTrans(){
            return _rowTrans;
        }
        std::vector<myBitSet<M>> getColumnTrans(){
            return _columnTrans;
        }
        //int getRow(){
            //return _row;
        //}
        //int getColumn(){
            //return _column;
        //}
    private:
        std::ifstream _fBinPattern;
        std::vector<myBitSet<N>> _rowTrans;
        std::vector<myBitSet<M>> _columnTrans;
        //int _row,_column;
};

class CCalculatePattern{
    public:
        CCalculatePattern( const std::vector<myBitSet<N>> & rowTrans,
                const std::vector<myBitSet<M>> & columnTrans ):
            _rowTrans(rowTrans), _columnTrans(columnTrans), _rowPatternLen(M){
            if(_rowTrans.empty() || _columnTrans.empty())
            {
                std::cerr<<"Error in CCalculatePattern: Input parameters vector is empty!"<<std::endl;
                return;
            }
            for(int i=0; i<M; ++i){
                _rowPatternLen[i] = _rowTrans[i].SSE4_count();
            }
        }
        ~CCalculatePattern(){}
        CCalculatePattern( const CCalculatePattern & ) = delete;
        CCalculatePattern & operator=( const CCalculatePattern & ) = delete;

        //Calculate support and occupancy of given pattern
        std::vector<double> calculate( const myBitSet<N> & pattern ){
            double sup = 0.0, occ = 0.0;
            std::vector<int> supIndex = pattern.getIndices();
            std::vector<double> resVect{sup,occ};
            if(supIndex.empty())
            {
                return resVect;
            }
            myBitSet<M> res;
            res.set();
            for(auto i : supIndex)
                res &= _columnTrans[i];
            int supCount = res.SSE4_count();
            if(!supCount){
                return resVect;
            }
            double thisPatternLength = static_cast<double>(pattern.SSE4_count());
            resVect[0] = static_cast<double>(supCount) / static_cast<double>(M);
            double sum = 0.0;
            std::vector<int> resIndex = res.getIndices();
            for(auto i : resIndex){
                sum += thisPatternLength / static_cast<double>(_rowPatternLen[i]);
            }
            resVect[1] = sum / static_cast<double>(supCount);
            return resVect;
        }
    private:
        std::vector<myBitSet<N>> _rowTrans;
        std::vector<myBitSet<M>> _columnTrans;
        std::vector<int> _rowPatternLen;
};

class CGrayCode{
    public:
        CGrayCode():_nFullNum(1ul<<(N)), _fullGrayCoding(){
            _fullGrayCoding.reserve(_nFullNum);
            generateGrayCode();
        }
        ~CGrayCode(){}
        CGrayCode(const CGrayCode&) = delete;
        CGrayCode & operator=(const CGrayCode &) = delete;

        void displayGrayCode(){
            if(_fullGrayCoding.empty())
            {
                std::cerr<<"Gray Code Vector is Empty!"<<std::endl;
                return;
            }
            for(auto i : _fullGrayCoding){
                std::cout<<i<<std::endl;
            }
        }
        std::vector<myBitSet<N>> getGrayCode(){
            return _fullGrayCoding;
        }
    private:
        const unsigned long _nFullNum;
        std::vector<myBitSet<N>> _fullGrayCoding;

        void generateGrayCode(){
            for(unsigned long  i=0; i<_nFullNum; ++i){
                myBitSet<N> binCode(i);
                myBitSet<N> grayCode(i>>1);
                grayCode ^= binCode;
                _fullGrayCoding.push_back(grayCode);
            }
        }

};
#endif
