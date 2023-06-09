#include<iostream>
#include<fstream>
#include "gamma.h"   //声明全局变量gamma as g
#include "S_invS.h"
#include "Difference_Scheme.h"
#include "Matrix.h"  //矩阵运算
#include "FVS.h"  //通量矢量分裂相关函数
#include "FDS.h" //通量差分分裂相关函数 
using namespace std;
double g=1.4;
int main(){
    //开始计算
    int Nx = 101;
    double dx = 1./(Nx-1),dt = 0.001;
    double t_end = 0.14,Nt = round(t_end/dt);
    cout << "Nt=" << Nt << endl;
    double U[3][Nx];
    
    for(int i=0;i<Nx;i++){
        //cout << i*dx << endl;
        U[0][i] = (i>Nx/2 ? 0.125 : 1);  //rho
        U[1][i] = 0.;  // rho*u
        U[2][i] = (i>Nx/2 ? 0.1 : 1)/(g-1); //E = 1/2*rho*u^2 + p/(g-1)
    }
    updateU(U[0],Nx,dx,dt,Nt,WENO5Z,5,2,LFf);  //LF分裂 使用Weno格式，5个基架点，第一个点为j-2

    //输出
    ofstream csvfile;
    csvfile.open("LF-weno5Z_t=0.14.csv", ios::out | ios::trunc);
    csvfile <<"x" << "," << "rho"<<","<< "u" <<","<<"p"<< endl;
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx <<","<< U[0][i]<<","<<U[1][i]/U[0][i]<<","<<(g-1)*(U[2][i]-0.5*U[1][i]*U[1][i]/U[0][i]) << endl;
    }
    csvfile.close();
    //system("pause");
    return 0;
}