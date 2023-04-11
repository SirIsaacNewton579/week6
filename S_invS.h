#include<math.h>
void Su(double *U,double *S){
    double rho=U[0],u=U[1]/U[0],p=(g-1)*(U[2]-0.5*rho*u*u);
    double c = sqrt(g*p/rho);
    double h = 0.5*u*u + g/(g-1)*p/rho;
    S[0] = 0.5*u*u - c*c/(g-1);
    S[1] = -u;
    S[2] = 1.;
    S[3] = -u -(g-1)/c*0.5*u*u;
    S[4] = 1+(g-1)/c*u;
    S[5] = -(g-1)/c;
    S[6] = -u+(g-1)/c*0.5*u*u;
    S[7] = 1-(g-1)/c*u;
    S[8] = (g-1)/c;
}
void invSu(double *U,double *invS){
    double rho=U[0],u=U[1]/U[0],p=(g-1)*(U[2]-0.5*rho*u*u);
    double c = sqrt(g*p/rho);
    double h = 0.5*u*u + g/(g-1)*p/rho;
    invS[0] = -(g-1)/(c*c);
    invS[1] = -1./(2*c);
    invS[2] = 1./(2*c);
    invS[3] = -(g-1)/(c*c)*u;
    invS[4] = -(u-c)/(2*c);
    invS[5] = (u+c)/(2*c);
    invS[6] = -(g-1)/(c*c)*0.5*u*u;
    invS[7] = -1./(2*c)*(h-u*c);
    invS[8] = 1./(2*c)*(h+u*c);
}