#include<iostream>
#include<fstream>
#include<math.h>
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
    double L=5;
    double rho,u,p;
    int Nx = 2001;
    double dx = 2.0*L/(Nx-1),dt = 0.0001;
    double t_end = 0.25,Nt = round(t_end/dt);
    cout << "Nt=" << Nt << endl;
    double U[3][Nx];
    
    for(int i=0;i<Nx;i++){
        if(i*dx-L<-4.0){
            rho=3.857;u=2.629;p=10.333;
        }
        else{
            rho=1+0.3*sin(40.0*(i*dx-L));
            u=0.0;p=1.0;
        }
        U[0][i] = rho;  //rho
        U[1][i] = rho*u;  // rho*u
        U[2][i] = 0.5*rho*u*u+p/(g-1); //E = 1/2*rho*u^2 + p/(g-1)
    }
    updateU(U[0],Nx,dx,dt,Nt,WENO5,5,2,SWf);  //LF分裂 使用Weno格式，5个基架点，第一个点为j-2

    //输出
    ofstream csvfile;
    csvfile.open("Shu_Osher_Nx=2001_SW-weno5_t=0.25.csv", ios::out | ios::trunc);
    csvfile <<"x" << "," << "rho"<<","<< "u" <<","<<"p"<< endl;
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx-L <<","<< U[0][i]<<","<<U[1][i]/U[0][i]<<","<<(g-1)*(U[2][i]-0.5*U[1][i]*U[1][i]/U[0][i]) << endl;
    }
    csvfile.close();
    //system("pause");
    return 0;
}