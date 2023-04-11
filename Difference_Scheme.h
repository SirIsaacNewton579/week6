#include<math.h>
using namespace std;
void WENO5JS(double *f,double &fm,int pm){
    int j=2,p=2;
    double eps = 1e-6;
    double fpi[3] = {f[j-2*pm]/3.0 - 7.0*f[j-1*pm]/6.0 + 11.0/6.0*f[j],
                -f[j-1*pm]/6.0 + 5.0*f[j]/6.0 + 1.0/3.0*f[j+1*pm],
                f[j]/3.0 + 5.0*f[j+1*pm]/6.0 - 1.0/6.0*f[j+2*pm]};
    double IS[3] = {pow(f[j-2*pm]-4.0*f[j-1*pm]+3*f[j],2)/4.0 + pow(f[j-2*pm]-2*f[j-1*pm]+f[j],2)*13.0/12.0,
                    pow(f[j-1*pm]-f[j+1*pm],2)/4.0 + pow(f[j-1*pm]-2*f[j]+f[j+1*pm],2)*13.0/12.0,
                    pow(f[j]-4.0*f[j+1*pm]+3*f[j+2*pm],2)/4.0 + pow(f[j]-2*f[j+1*pm]+f[j+2*pm],2)*13.0/12.0};
    double C[3] = {1.0/10.0,6.0/10.0,3.0/10.0};
    double alp[3] = {C[0]/pow(eps+IS[0],p),C[1]/pow(eps+IS[1],p),C[2]/pow(eps+IS[2],p)};
    double alps = alp[0] + alp[1] + alp[2];
    double w[3] = {alp[0]/alps,alp[1]/alps,alp[2]/alps};
    fm = w[0]*fpi[0] + w[1]*fpi[1] + w[2]*fpi[2];     
}
void WENO5M(double *f,double &fm,int pm){
    int j=2,p=2;
    double eps = 1e-6;
    double fpi[3] = {f[j-2*pm]/3.0 - 7.0*f[j-1*pm]/6.0 + 11.0/6.0*f[j],
                -f[j-1*pm]/6.0 + 5.0*f[j]/6.0 + 1.0/3.0*f[j+1*pm],
                f[j]/3.0 + 5.0*f[j+1*pm]/6.0 - 1.0/6.0*f[j+2*pm]};
    double IS[3] = {pow(f[j-2*pm]-4.0*f[j-1*pm]+3*f[j],2)/4.0 + pow(f[j-2*pm]-2*f[j-1*pm]+f[j],2)*13.0/12.0,
                    pow(f[j-1*pm]-f[j+1*pm],2)/4.0 + pow(f[j-1*pm]-2*f[j]+f[j+1*pm],2)*13.0/12.0,
                    pow(f[j]-4.0*f[j+1*pm]+3*f[j+2*pm],2)/4.0 + pow(f[j]-2*f[j+1*pm]+f[j+2*pm],2)*13.0/12.0};
    double C[3] = {1.0/10.0,6.0/10.0,3.0/10.0};
    double alp[3] = {C[0]/pow(eps+IS[0],p),C[1]/pow(eps+IS[1],p),C[2]/pow(eps+IS[2],p)};
    double alps = alp[0] + alp[1] + alp[2];
    double w[3] = {alp[0]/alps,alp[1]/alps,alp[2]/alps};
    for(int k=0;k<3;k++){
        alp[k] = w[k]*(C[k]+C[k]*C[k]-3*C[k]*w[k]+w[k]*w[k])/(C[k]*C[k]+w[k]*(1-2*C[k]));
    }
    alps = alp[0] + alp[1] + alp[2];
    for(int k=0;k<3;k++) w[k] = alp[k]/alps;
    fm = w[0]*fpi[0] + w[1]*fpi[1] + w[2]*fpi[2];     
}
void WENO5Z(double *f,double &fm,int pm){
    int j=2,p=2;
    double eps = 1e-6;
    double fpi[3] = {f[j-2*pm]/3.0 - 7.0*f[j-1*pm]/6.0 + 11.0/6.0*f[j],
                -f[j-1*pm]/6.0 + 5.0*f[j]/6.0 + 1.0/3.0*f[j+1*pm],
                f[j]/3.0 + 5.0*f[j+1*pm]/6.0 - 1.0/6.0*f[j+2*pm]};
    double IS[3] = {pow(f[j-2*pm]-4.0*f[j-1*pm]+3*f[j],2)/4.0 + pow(f[j-2*pm]-2*f[j-1*pm]+f[j],2)*13.0/12.0,
                    pow(f[j-1*pm]-f[j+1*pm],2)/4.0 + pow(f[j-1*pm]-2*f[j]+f[j+1*pm],2)*13.0/12.0,
                    pow(f[j]-4.0*f[j+1*pm]+3*f[j+2*pm],2)/4.0 + pow(f[j]-2*f[j+1*pm]+f[j+2*pm],2)*13.0/12.0};
    double C[3] = {1.0/10.0,6.0/10.0,3.0/10.0};
    double tau5 = abs(IS[0]-IS[2]);
    double alp[3] ;
    for(int k=0;k<3;k++)
        alp[k] = C[k]*(1+pow(tau5/(eps+IS[k]),p));
    double alps = alp[0] + alp[1] + alp[2];
    double w[3] = {alp[0]/alps,alp[1]/alps,alp[2]/alps};
    fm = w[0]*fpi[0] + w[1]*fpi[1] + w[2]*fpi[2];     
}
void first_windward(double *f,double &fm,int pm) {fm = f[0];}
void VanLeerLim(double *ft,double &ftm,int pm){
    //pm = 1左值，pm=-1右值
    int j=1;
    double dL = ft[j]-ft[j-1*pm],dR=ft[j+1*pm]-ft[j];
    double eps = 1e-12;
    double phi = (abs(dL*dR) + dL*dR+eps)/(abs(dL*dR)+dR*dR+eps);
    ftm = ft[j]+0.5*phi*dR;
}
double minmod(double a,double b){
    if(!((a>0.0)^(b>0.0))) return (abs(a)>abs(b)?b:a);
    else return 0.0;
}
void NND(double *ft,double &ftm,int pm){
    //pm = 1左值，pm=-1右值
    int j=1;
    ftm = ft[j] + 0.5*minmod(ft[j]-ft[j-1*pm],ft[j+1*pm]-ft[j]);
}
void GVC2(double *ft,double &ftm,int pm){
    int j=1;
    if(abs(ft[j]-ft[j-1*pm]) < abs(ft[j+1*pm]-ft[j])){
        ftm = 0.5*(3*ft[j]-ft[j-1*pm]);
    }
    else{
        ftm = 0.5*(ft[j] + ft[j+1*pm]);
    }
}
void MUSCL(double *ft,double &ftm,int pm){
    int j=1;
    double dmf=ft[j]-ft[j-1*pm],dpf=ft[j+1*pm]-ft[j];
    double eps = 1e-6;
    double s = (2*dmf*dpf+eps)/(dmf*dmf+dpf*dpf+eps);
    ftm = ft[j]+s/4.0*((1-s/3.0)*dmf+(1+s/3.0)*dpf);
}