#include <iostream>
#include <array>
#include <math.h>

double angle(double Rt[],int16_t n){
    for(int16_t i = 0; i < n; i++){
        Rt[i] = Rt[i] * M_PI / 180;
    }
    return 0;
}
double degree(double Rt[],int16_t n){
    for(int16_t i = 0; i < n; i++){
        Rt[i] = Rt[i] * 180 / M_PI;
    }
    return 0;
}
double IK(double L[], double xt[], double Rt[],double psi, double theta[]){
    
    angle(Rt,3);
    
    double x_sw[3];
    x_sw[0] = xt[0] - cos(Rt[0]) * sin(Rt[1]) * L[3];
    x_sw[1] = xt[1] - sin(Rt[0]) * sin(Rt[1]) * L[3];
    x_sw[2] = xt[2] - L[0] - cos(Rt[1]) * L[3];
    double nor_xsw2;
    nor_xsw2 = x_sw[0]*x_sw[0] + x_sw[1]*x_sw[1] + x_sw[2]*x_sw[2];
    double r[3][3];
    r[0][0] = cos(Rt[1])*cos(Rt[2]) ;  r[0][1] = -1*sin(Rt[1])*sin(Rt[2]) ;  r[0][2] = sin(Rt[1]);
    r[1][0] = sin(Rt[0])*sin(Rt[1])*cos(Rt[2]) + cos(Rt[1])*sin(Rt[2]);  r[1][1] = -1*sin(Rt[0])*sin(Rt[1])*sin(Rt[2]) + cos(Rt[0])*cos(Rt[2]);  r[1][2] = -1*sin(Rt[0])*cos(Rt[1]);
    r[2][0] = -1*cos(Rt[0])*sin(Rt[1])*cos(Rt[2]) + sin(Rt[0])*sin(Rt[2]);  r[2][1] = cos(Rt[0])*sin(Rt[1])*sin(Rt[2]) + sin(Rt[0])*cos(Rt[2]);  r[2][2] = cos(Rt[0])*cos(Rt[1]);
    double u_sw[3];
    u_sw[0] = x_sw[0]/sqrt(nor_xsw2);
    u_sw[1] = x_sw[1]/sqrt(nor_xsw2);
    u_sw[2] = x_sw[2]/sqrt(nor_xsw2);

    theta[3] = acos((nor_xsw2 - L[1]*L[1] - L[2]*L[2])/(2*L[1]*L[2]));

    double theta20,theta10,S2,C2,M,N;
    theta10 = atan2( x_sw[1] , x_sw[0] );
    N = sin(theta[3])*L[2];
    M = cos(theta[3])*L[2] + L[1];
    S2 = -1*( M*sqrt( x_sw[0]*x_sw[0] + x_sw[1]*x_sw[1] ) + N*x_sw[2] )/( N*N + M*M );
    C2 = ( N*sqrt( x_sw[0]*x_sw[0] + x_sw[1]*x_sw[1] ) - M*x_sw[2] )/( N*N + M*M );
    theta20 = atan2( S2 , C2 );

    double As[3][3],Bs[3][3],Cs[3][3];
    As[0][0] = -1*u_sw[2]*r[1][0] + u_sw[1]*r[2][0];
    As[0][1] = -1*u_sw[2]*r[1][1] + u_sw[1]*r[2][1];
    As[1][0] = u_sw[2]*r[0][0] - u_sw[0]*r[2][0];
    As[1][1] = u_sw[2]*r[0][1] - u_sw[0]*r[2][1];
    As[2][0] = u_sw[0]*r[1][0] - u_sw[1]*r[0][0];
    As[2][1] = u_sw[0]*r[1][1] - u_sw[1]*r[0][1];
    As[2][2] = u_sw[0]*r[1][2] - u_sw[1]*r[0][1];

    Bs[0][0] = -1*u_sw[0]*u_sw[1]*r[1][0] - u_sw[2]*u_sw[0]*r[2][0] + r[0][0]*( u_sw[2]*u_sw[2] + u_sw[1]*u_sw[1] );
    Bs[0][1] = -1*u_sw[0]*u_sw[2]*r[2][1] - u_sw[1]*u_sw[0]*r[1][1] + r[0][1]*( u_sw[2]*u_sw[2] + u_sw[1]*u_sw[1] );
    Bs[1][0] = -1*u_sw[2]*u_sw[1]*r[2][0] - u_sw[1]*u_sw[0]*r[0][0] + r[1][0]*( u_sw[2]*u_sw[2] + u_sw[1]*u_sw[1] );
    Bs[1][1] = -1*u_sw[1]*u_sw[2]*r[2][1] - u_sw[1]*u_sw[0]*r[0][1] + r[1][1]*( u_sw[2]*u_sw[2] + u_sw[1]*u_sw[1] );
    Bs[2][0] = -1*u_sw[0]*u_sw[2]*r[0][0] - u_sw[1]*u_sw[2]*r[1][0] + r[2][0]*( u_sw[0]*u_sw[0] + u_sw[1]*u_sw[1] );
    Bs[2][1] = -1*u_sw[0]*u_sw[2]*r[0][1] - u_sw[1]*u_sw[2]*r[1][1] + r[2][1]*( u_sw[0]*u_sw[0] + u_sw[1]*u_sw[1] );
    Bs[2][2] = -1*u_sw[0]*u_sw[2]*r[0][2] - u_sw[1]*u_sw[2]*r[1][2] + r[2][2]*( u_sw[0]*u_sw[0] + u_sw[1]*u_sw[1] );

    Cs[0][0] = u_sw[0]*u_sw[0]*r[0][0] + u_sw[1]*u_sw[0]*r[1][0] + r[2][0]*u_sw[2]*u_sw[0];
    Cs[0][1] = u_sw[0]*u_sw[0]*r[0][1] + u_sw[1]*u_sw[0]*r[1][1] + r[2][1]*u_sw[2]*u_sw[0];
    Cs[1][0] = u_sw[1]*u_sw[0]*r[0][0] + u_sw[1]*u_sw[1]*r[1][0] + r[2][0]*u_sw[2]*u_sw[1];
    Cs[1][1] = u_sw[0]*u_sw[1]*r[0][1] + u_sw[1]*u_sw[1]*r[1][1] + r[2][1]*u_sw[2]*u_sw[1];
    Cs[2][0] = u_sw[0]*u_sw[2]*r[0][0] + u_sw[1]*u_sw[2]*r[1][0] + r[2][0]*u_sw[2]*u_sw[2];
    Cs[2][1] = u_sw[0]*u_sw[2]*r[0][1] + u_sw[1]*u_sw[2]*r[1][1] + r[2][1]*u_sw[2]*u_sw[2];
    Cs[2][2] = u_sw[0]*u_sw[2]*r[0][2] + u_sw[1]*u_sw[2]*r[1][2] + r[2][2]*u_sw[2]*u_sw[2];

    double theta01,theta02,theta1,theta21,theta22;
    theta01 = -1*( As[1][1]*sin(psi) + Bs[1][1]*cos(psi) + Cs[1][1] );
    theta02 = -1*( As[0][1]*sin(psi) + Bs[0][1]*cos(psi) + Cs[0][1] );
    theta1  = -1*( As[2][1]*sin(psi) + Bs[2][1]*cos(psi) + Cs[2][1] );
    theta21 =      As[2][2]*sin(psi) + Bs[2][2]*cos(psi) + Cs[2][2];
    theta22 = -1*( As[2][0]*sin(psi) + Bs[2][0]*cos(psi) + Cs[2][0] );
    std::cout << "theta01 " << theta01 << " theta02 " << theta02 << " theta1 " << theta1 << " theta21 " << theta21 << " theta22 " << theta22 << std::endl;
    theta[0] = atan2( theta01 , theta02);
    theta[1] = acos( theta1 );
    theta[2] = atan2( theta21 , theta22 );

    double Aw[3][3],Bw[3][3],Cw[3][3];
    Aw[0][2] = r[0][2]*( As[0][0]*cos(L[3]) + As[0][1]*sin(L[3]) ) + r[1][2]*( As[1][0]*cos(L[3]) + As[1][1]*sin(L[3]) ) + r[2][2]*( As[2][0]*cos(L[3]) + As[2][1]*sin(L[3]) );
    Aw[1][2] = r[0][2]*As[0][1] + r[1][2]*As[1][1] + r[2][2]*As[2][1];
    Aw[2][0] = r[0][0]*( As[0][0]*sin(L[3]) - As[0][1]*cos(L[3]) ) + r[1][0]*( As[1][0]*sin(L[3]) - As[1][1]*cos(L[3]) ) + r[2][0]*( As[2][0]*sin(L[3]) + As[2][1]*cos(L[3]) );
    Aw[2][1] = r[0][1]*( As[0][0]*sin(L[3]) - As[0][1]*cos(L[3]) ) + r[1][1]*( As[1][0]*sin(L[3]) - As[1][1]*cos(L[3]) ) + r[2][1]*( As[2][0]*sin(L[3]) + As[2][1]*cos(L[3]) );
    Aw[2][2] = r[0][2]*( As[0][0]*sin(L[3]) - As[0][1]*cos(L[3]) ) + r[1][2]*( As[1][0]*sin(L[3]) - As[1][1]*cos(L[3]) ) + r[2][2]*( As[2][0]*sin(L[3]) + As[2][1]*cos(L[3]) );

    Bw[0][2] = r[0][2]*( Bs[0][0]*cos(L[3]) + Bs[0][1]*sin(L[3]) ) + r[1][2]*( Bs[1][0]*cos(L[3]) + Bs[1][1]*sin(L[3]) ) + r[2][2]*( Bs[2][0]*cos(L[3]) + Bs[2][1]*sin(L[3]) );
    Bw[1][2] = r[0][2]*Bs[0][1] + r[1][2]*Bs[1][1] + r[2][2]*Bs[2][1];
    Bw[2][0] = r[0][0]*( Bs[0][0]*sin(L[3]) - Bs[0][1]*cos(L[3]) ) + r[1][0]*( Bs[1][0]*sin(L[3]) - Bs[1][1]*cos(L[3]) ) + r[2][0]*( Bs[2][0]*sin(L[3]) + Bs[2][1]*cos(L[3]) );
    Bw[2][1] = r[0][1]*( Bs[0][0]*sin(L[3]) - Bs[0][1]*cos(L[3]) ) + r[1][1]*( Bs[1][0]*sin(L[3]) - Bs[1][1]*cos(L[3]) ) + r[2][1]*( Bs[2][0]*sin(L[3]) + Bs[2][1]*cos(L[3]) );
    Bw[2][2] = r[0][2]*( Bs[0][0]*sin(L[3]) - Bs[0][1]*cos(L[3]) ) + r[1][2]*( Bs[1][0]*sin(L[3]) - Bs[1][1]*cos(L[3]) ) + r[2][2]*( Bs[2][0]*sin(L[3]) + Bs[2][1]*cos(L[3]) );

    Cw[0][2] = r[0][2]*( Cs[0][0]*cos(L[3]) + Cs[0][1]*sin(L[3]) ) + r[1][2]*( Cs[1][0]*cos(L[3]) + Cs[1][1]*sin(L[3]) ) + r[2][2]*( Cs[2][0]*cos(L[3]) + Cs[2][1]*sin(L[3]) );
    Cw[1][2] = r[0][2]*Cs[0][1] + r[1][2]*Cs[1][1] + r[2][2]*Cs[2][1];
    Cw[2][0] = r[0][0]*( Cs[0][0]*sin(L[3]) - Cs[0][1]*cos(L[3]) ) + r[1][0]*( Cs[1][0]*sin(L[3]) - Cs[1][1]*cos(L[3]) ) + r[2][0]*( Cs[2][0]*sin(L[3]) + Cs[2][1]*cos(L[3]) );
    Cw[2][1] = r[0][1]*( Cs[0][0]*sin(L[3]) - Cs[0][1]*cos(L[3]) ) + r[1][1]*( Cs[1][0]*sin(L[3]) - Cs[1][1]*cos(L[3]) ) + r[2][1]*( Cs[2][0]*sin(L[3]) + Cs[2][1]*cos(L[3]) );
    Cw[2][2] = r[0][2]*( Cs[0][0]*sin(L[3]) - Cs[0][1]*cos(L[3]) ) + r[1][2]*( Cs[1][0]*sin(L[3]) - Cs[1][1]*cos(L[3]) ) + r[2][2]*( Cs[2][0]*sin(L[3]) + Cs[2][1]*cos(L[3]) );

    theta[4] = atan2( Aw[1][2]*sin(psi) + Bw[1][2]*cos(psi) + Cw[1][2] , Aw[0][2]*sin(psi) + Bw[0][2]*cos(psi) + Cw[0][2] );
    theta[5] = acos( Aw[2][2]*sin(psi) + Bw[2][2]*cos(psi) + Cw[2][2] );
    theta[6] = atan2( Aw[2][1]*sin(psi) + Bw[2][1]*cos(psi) + Cw[2][1] , -1*( Aw[2][0]*sin(psi) + Bw[2][0]*cos(psi) + Cw[2][0] ) );

    degree(theta,7);
    return 0;
}

int main(void){
    double L[4] = { 0.105,0.250,0.250,0.060 };
    double xt[3] = { 0.353,0,0.105 };
    double Rt[3] = { 0,0,90 };
    double theta[7] = {0};
    double psi;psi=30*M_PI/180;

    IK(L,xt,Rt,psi,theta);
    for(auto ar : theta){
        std::cout << ar << std::endl;
    }
    return 0;
}
