#include<iostream>
#include<math.h>
typedef void Scheme(double *,double &,int );
void fU(double *U,double *f){
    //计算f(U)
    f[0] = U[1];
    f[1]=(g-1)*U[2]+0.5*(3-g)*U[1]*U[1]/U[0];
    f[2]=g*U[2]*U[1]/U[0]-0.5*(g-1)*U[1]*U[1]*U[1]/(U[0]*U[0]);
}
void Ubar(double *UL,double *UR,double *Ub){
    //用左值右值计算roe点
    double rhoL=UL[0],uL=UL[1]/UL[0],pL=(g-1)*(UL[2]-0.5*rhoL*uL*uL);
    double hL = 0.5*uL*uL + g/(g-1)*pL/rhoL;

    double rhoR=UR[0],uR=UR[1]/UR[0],pR=(g-1)*(UR[2]-0.5*rhoR*uR*uR);
    double hR = 0.5*uR*uR + g/(g-1)*pR/rhoR;

    double rhob = pow((sqrt(rhoL)+sqrt(rhoR))/2,2);
    double ub = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(2*sqrt(rhob));
    double hb = (sqrt(rhoL)*hL+sqrt(rhoR)*hR)/(2*sqrt(rhob));
    double pb = (g-1)/g*rhob*(hb - 0.5*ub*ub);
    Ub[0] = rhob;Ub[1] = rhob*ub; Ub[2] = rhob*hb-pb;
}
inline double abss(double lam) {return (abs(lam)>1e-2 ? abs(lam) : (lam*lam+1e-4)/(2*1e-2));}
void LBD_abs(double *U,double* LBDa){
    //计算|Lambda|
    double rho,u,p,c;
    rho = U[0];u = U[1]/U[0];p =(g-1)*(U[2]-0.5*rho*u*u);
    c = sqrt(g*p/rho);
    LBDa[0] = abss(u);LBDa[4] = abss(u-c);LBDa[8] = abss(u+c);
}
void Flux_e_fds(double *U,int Nx,double *fo,Scheme DS,int bpn,int sn,int bopn){
    //通量差分分裂，Roe格式
    double UL[3],UR[3],VL[3],VR[3],Ub[3],fUL[3],fUR[3];
    double S[9],invS[9],LBDa[9] = {0.0};
    double Ab[9];
    double tmp[9],tmp2[3],tmp3[3];
    int i,j,k;
    double Ut1[3][bpn],Ut2[3][bpn];
    double Vt1[3][bpn],Vt2[3][bpn];
    double thisU[3];
    for(j=bopn;j<Nx-1-bopn;j++){
        for(i=0;i<3;i++) thisU[i] = 0.5*(U[j+i*Nx]+U[j+1+i*Nx]);
        Su(thisU,S); //更新S_j+1/2
        invSu(thisU,invS); //更新S^-1_j+1/2
        for(k=j-sn;k<=j-sn+bpn-1;k++){
            for(i=0;i<3;i++) Ut1[i][k-j+sn]=U[i*Nx+k];
        }
        for(k=j-(bpn-1-sn)+1;k<=j+sn+1;k++){
            for(i=0;i<3;i++) Ut2[i][k-(j-(bpn-1-sn)+1)] = U[i*Nx+k];
        }
        mul_matrix(S,Ut1[0],Vt1[0],3,3,bpn);
        mul_matrix(S,Ut2[0],Vt2[0],3,3,bpn);
        for(i=0;i<3;i++) {DS(Vt1[i],VL[i],1);DS(Vt2[i],VR[i],-1);}
        mul_matrix(invS,VL,UL,3,3,1);
        mul_matrix(invS,VR,UR,3,3,1);
        fU(UL,fUL);
        fU(UR,fUR);
        Ubar(UL,UR,Ub);
        Su(Ub,S);
        invSu(Ub,invS);
        LBD_abs(Ub,LBDa);
        mul_matrix(invS,LBDa,tmp); // tmp = S^-1*Lambda 
        mul_matrix(tmp,S,Ab);  //|A| = S^-1*Lambda*S
        minus_matrix(UR,UL,tmp2);  //tmp2 = UR-UL
        mul_matrix(Ab,tmp2,tmp3,3,3,1);  //tmp3 = |A|*(UR-UL)
        add_matrix(fUL,fUR,tmp2); //tmp2 = f(UR)+f(UL)
        minus_matrix(tmp2,tmp3,tmp2); //tmp2 = f(UR)+f(UL)-|A|*(UR-UL)
        mul_matrix(tmp2,0.5);
        for(i=0;i<3;i++) fo[i*(Nx-2*bopn-1)+j-bopn] = tmp2[i];
    }
}
void Flux_fds(double *U,int Nx,double *fo,Scheme DS,int bpn,int sn,int bopn){
    //通量差分分裂，Roe格式
    double UL[3],UR[3],Ub[3],fUL[3],fUR[3];
    double S[9],invS[9],LBDa[9] = {0.0};
    double Ab[9];
    double tmp[9],tmp2[3],tmp3[3];
    int i,j,k;
    double Ut1[3][bpn],Ut2[3][bpn];
    for(j=bopn;j<Nx-1-bopn;j++){
        for(k=j-sn;k<=j-sn+bpn-1;k++){
            for(i=0;i<3;i++) Ut1[i][k-j+sn]=U[i*Nx+k];
        }
        for(k=j-(bpn-1-sn)+1;k<=j+sn+1;k++){
            for(i=0;i<3;i++) Ut2[i][k-(j-(bpn-1-sn)+1)] = U[i*Nx+k];
        }
        for(i=0;i<3;i++) {DS(Ut1[i],UL[i],1);DS(Ut2[i],UR[i],-1);}
        fU(UL,fUL);
        fU(UR,fUR);
        Ubar(UL,UR,Ub);
        Su(Ub,S);
        invSu(Ub,invS);
        LBD_abs(Ub,LBDa);
        mul_matrix(invS,LBDa,tmp); // tmp = S^-1*Lambda 
        mul_matrix(tmp,S,Ab);  //|A| = S^-1*Lambda*S
        minus_matrix(UR,UL,tmp2);  //tmp2 = UR-UL
        mul_matrix(Ab,tmp2,tmp3,3,3,1);  //tmp3 = |A|*(UR-UL)
        add_matrix(fUL,fUR,tmp2); //tmp2 = f(UR)+f(UL)
        minus_matrix(tmp2,tmp3,tmp2); //tmp2 = f(UR)+f(UL)-|A|*(UR-UL)
        mul_matrix(tmp2,0.5);
        for(i=0;i<3;i++) fo[i*(Nx-2*bopn-1)+j-bopn] = tmp2[i];
    }
}
void updateU(double *U,int Nx,double dx,double dt,int Nt,Scheme DS,int bpn,int sn){
    //bpn ： 基架点数量 ； sn（shift num)：格式向左偏移数
    int bopn = (sn>bpn-1-sn ? sn : bpn-1-sn); //边界点数量
    double fo[3][Nx-2*bopn-1];  
    double U1[3][Nx],U2[3][Nx],Unext[3][Nx];
    double cfl = dt/dx;
    int i,j,N;
    //时间步推进
    for(N=1;N<=Nt;N++){
        Flux_e_fds(U,Nx,fo[0],DS,bpn,sn,bopn); //计算f_j+1/2
        for(i=0;i<3;i++) {
            for(j=0;j<=bopn;j++) {
                U1[i][j] = U[j+i*Nx];
                U1[i][Nx-j-1] = U[i*Nx+Nx-j-1];
            }
        }
        for(j=bopn+1;j<Nx-bopn-1;j++){
            for(i=0;i<3;i++){
                U1[i][j] = U[j+i*Nx] - cfl*(fo[i][j-bopn]-fo[i][j-bopn-1]);
            }    
        }
        Flux_e_fds(U1[0],Nx,fo[0],DS,bpn,sn,bopn);
        for(i=0;i<3;i++) {
            for(j=0;j<=bopn;j++) {
                U2[i][j] = U1[i][j];
                U2[i][Nx-j-1] = U1[i][Nx-j-1];
            }
        }
        for(j=bopn+1;j<Nx-bopn-1;j++){
            for(i=0;i<3;i++){
                U2[i][j] = 0.75*U[j+i*Nx] + 0.25*(U1[i][j]- cfl*(fo[i][j-bopn]-fo[i][j-bopn-1]));
            }
        }
        Flux_e_fds(U2[0],Nx,fo[0],DS,bpn,sn,bopn);
        for(i=0;i<3;i++) {
            for(j=0;j<=bopn;j++) {
                Unext[i][j] = U2[i][j];
                Unext[i][Nx-j-1] = U2[i][Nx-j-1];
            }
        }
        for(j=bopn+1;j<Nx-bopn-1;j++){
            for(i=0;i<3;i++){
                Unext[i][j] = 1.0*U[j+i*Nx]/3.0 + 2.0/3.0*(U2[i][j]- cfl*(fo[i][j-bopn]-fo[i][j-bopn-1]));
            }
        }
        for(i=0;i<3;i++){
            for(j=0;j<Nx;j++){
                U[j+i*Nx] = Unext[i][j];
            }
        }
    }
}