// ---- Program Info Start----
//FileName:     RBFKmean.cpp
//
//Author:       Fuchen Duan
//
//Email:        slow295185031@gmail.com
//
//CreatedAt:    2016-12-07 13:39:24
// ---- Program Info End  ----

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <limits>
#include <cassert>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <map>
#include <iomanip>
#include "matrix.h"

using namespace std;

typedef map<vector<int>, vector<double>> Ind;

const int P=100;        //输入样本的数量
const int N=10;
vector<vector<int>> X(P);  //输入样本
Matrix<double> Y(P,1);        //输入样本对应的期望输出
Matrix<double> Z(P,1);
const int M=10;         //隐藏层节点数目
vector<vector<int>> center(M);       //M个Green函数的数据中心
vector<double> delta(M);        //M个Green函数的扩展常数
Matrix<double> Green(P,M);         //Green矩阵
Matrix<double> weight(M,1);       //权值矩阵
Ind DB;  //数据库

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
        issLine>>obj[0]>>delimiter>>obj[1];

        DB.insert( Ind::value_type( pattern, obj ) );
    }
}
/*产生指定区间上均匀分布的随机数*/
inline double uniform(double floor,double ceil){
    return floor+1.0*rand()/RAND_MAX*(ceil-floor);
}

/*产生区间[floor,ceil]上服从正态分布N[mu,sigma]的随机数*/
inline double RandomNorm(double mu,double sigma,double floor,double ceil){
    double x,prob,y;
    do{
        x=uniform(floor,ceil);
        prob=1/sqrt(2*M_PI*sigma)*exp(-1*(x-mu)*(x-mu)/(2*sigma*sigma));
        y=1.0*rand()/RAND_MAX;
    }while(y>prob);
    return x;
}

template<class Type1, class Type2>
int getDistance( const vector<Type1> & lhs, const vector<Type2> & rhs ){
    int res = 0;
    for( int i=0; i<N; ++i ){
        res += fabs(lhs[i] - rhs[i]);
    }
    return res;
}

double getOutput( const vector<int> &x ){
    double y = 0.0;
    for( int i=0; i<M; ++i ){
        int dis = getDistance( x, center[i] );
        y += weight.get(i,0)*exp(-1.0*dis*dis/(2*delta[i]*delta[i]));
    }
    return y;
}

/*产生输入样本*/
void generateSample(){
    for(int i=0;i<P;++i){
        for( int j=0; j<N; ++j ){
            if( static_cast<double>(rand())/RAND_MAX > 0.5 ) {
                X[i].push_back(1);
            }
            else
                X[i].push_back(0);
        }
        Z.put(i,0,DB[X[i]][1]);
    }
}

void test( int num ){
    vector<int> testX;
    for( int i=0; i<num; ++i ){
        testX.clear();
        for( int j=0; j<N; ++j ){
            if( static_cast<double>(rand())/RAND_MAX > 0.5 ){
                testX.push_back(1);
                cout<<"1 ";
            }
            else{
                testX.push_back(0);
                cout<<"0 ";
            }
        }
        int dis = 0;
        for( int j=0; j<M; ++j ){
            dis += getDistance( testX, center[j] );
        }
        cout<<":"<<setw(12)<<dis<<setw(12)<<DB[testX][1]<<"\t"<<setw(12)<<getOutput(testX)<<"\t"<<setw(12)<<fabs(getOutput(testX)-DB[testX][1])<<endl;
    }
}

/*寻找样本离哪个中心最近*/
int nearest(const vector<vector<int> >& center, const vector<int> & sample){
    int idx=-1;
    int dist=numeric_limits<int>::max();
    for(int i=0;i<center.size();++i){
        int tmpdis = getDistance( sample, center[i] );
        if(tmpdis<dist){
            dist=tmpdis;
            idx=i;
        }
    }
    return idx;
}

/*计算簇的质心*/
vector<int> calCenter(const vector<vector<int> > &g){
    int len=g.size();
    vector<int> c(N);
    for( int j=0; j<N; ++j ){
        int cnt = 0;
        for(int i=0;i<len;++i) {
            cnt += g[i][j];
        }
        if( cnt > len-cnt )
            c[j] = 1;
        else
            c[j] = 0;
    }
    return c;
}

/*KMeans聚类法产生数据中心*/
void KMeans(){
    assert(P%M==0);
    vector<vector<vector<int> > > group(M);          //记录各个聚类中包含哪些样本
    int cnt = 0;
    for(int i=0;i<M;++i){   //从P个输入样本中随机选M个作为初始聚类中心
        for( int j=0; j<N; ++j ){
            center[i].push_back( X[10*i+3][j] );
        }
    }
    while(1){
        for(int i=0;i<M;++i)
            group[i].clear();   //先清空聚类信息
        for(int i=0;i<P;++i){       //把所有输入样本归到对应的簇
            int c=nearest(center,X[i]);
            group[c].push_back(X[i]);
        }
        vector<vector<int> > new_center(M);       //存储新的簇心
        for(int i=0;i<M;++i){
            vector<vector<int> > g=group[i];
            new_center[i]=calCenter(g);
        }
        bool flag=false;
        for(int i=0;i<M;++i){       //检查前后两次质心的改变量是否都小于gap
            int dis = getDistance( new_center[i], center[i] );
            if(dis != 0){
                flag=true;
                break;
            }
        }
        center=new_center;
        cnt++;
        if(!flag)
            break;
    }
}


/*生成Green矩阵*/
void calGreen(){
    for(int i=0;i<P;++i){
        for(int j=0;j<M;++j){
            int dis = getDistance(X[i], center[j]);
            Green.put(i,j,exp(-1.0*(dis)*(dis)/(2*delta[j]*delta[j])));
        }
    }
}

/*求一个矩阵的伪逆*/
Matrix<double> getGereralizedInverse(const Matrix<double> &matrix){
    return (matrix.getTranspose()*matrix).getInverse()*(matrix.getTranspose());
}

int main(int argc,char *argv[]){
    /*生成随机种子*/
    srand(time(0));

    string filename= "./PatternDB";

    readData( filename );

    generateSample();

    KMeans();

    //自组织选择法, 根据各中心之间最小距离确定扩展常数
    for( int i=0; i<M; ++i ){
        delta[i] = 1000;
        for( int j=0; j<M; ++j ){
            int dis = getDistance( center[i], center[j] );
            if( dis > -1e-10 && dis < 1e-10 ){ // 判断dis是否为0
            }
            else{
                delta[i] = delta[i] <= dis ? delta[i] : dis;
            }
        }
        if( delta[i] > 999.999 ){
            cerr<<"Error In Calculate Delta!"<<endl;
            exit(-1);
        }
    }

    calGreen();     //计算Green矩阵
    weight=getGereralizedInverse(Green)*Z;      //计算权值矩阵
    test(100);

    return 0;
}
