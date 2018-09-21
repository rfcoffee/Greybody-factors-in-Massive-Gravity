#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<stdio.h>

using namespace std;
const double pi = 3.1415926;
double RH = 0.1;                           // M is roughly 10
double h0 = 1.e-3; 
double hh = h0*1.e-5;

double r0 = RH*(1.+1.e-1);
double rend = 1.*(1.-1.e-1);
double omega0 = 10;
const int nw = 4;
double dw     = 10;
const int nL = 1;

// First define the preliminary functions.
double f(double r)
{
    return -(((-1 + r)*(r - RH)*(1 + r + RH))/(r*(1 + RH*(1 + RH))));
}

double df(double r)
{
    return -(((-1 + r)*(r - RH))/(r*(1 + RH*(1 + RH)))) - 
   ((-1 + r)*(1 + r + RH))/(r*(1 + RH*(1 + RH))) + 
   ((-1 + r)*(r - RH)*(1 + r + RH))/
    (pow(r,2)*(1 + RH*(1 + RH))) - 
   ((r - RH)*(1 + r + RH))/(r*(1 + RH*(1 + RH)));
}

double dr2f(double r)
{
    return (-2*(-1 + r)*(r - RH)*(1 + r + RH))/(1 + RH*(1 + RH)) + 
    pow(r,2)*(-(((-1 + r)*(r - RH))/(r*(1 + RH*(1 + RH)))) - 
      ((-1 + r)*(1 + r + RH))/(r*(1 + RH*(1 + RH))) + 
      ((-1 + r)*(r - RH)*(1 + r + RH))/
       (pow(r,2)*(1 + RH*(1 + RH))) - 
      ((r - RH)*(1 + r + RH))/(r*(1 + RH*(1 + RH))));
}

// Potential term in Regge-Wheeler equation
double V(double r, int L)
{
    return f(r)*df(r)/r + double(L*(L+1))*f(r)/r/r;	
}

double Function(double r, double yy, double dy, double omega, int L)
{
    return -yy*(pow(omega/f(r),2.0) - double(L*(L+1))/(pow(r,2.0)*f(r))) - dy*(dr2f(r)/(pow(r,2.0)*f(r)));
} 

// define the function to perform the 4th order Runge Kutta algorithm to update the values of phi and dphi
void rk4(double& y1, double& dy1, double& y2, double& dy2, double& r, double h, double omega, int L)
{
    double h2 = h/2.0;
    double k11,k12,k21,k22,k31,k32,k41,k42;
    k11 = h* dy1;
    k12 = h* Function(r,y1,dy1,omega,L);
    k21 = h* (dy1+k12/2.0);
    k22 = h* Function(r+h2, y1+k11/2.0,dy1+k12/2.0, omega,L);
    k31 = h* (dy1+k22/2.0);
    k32 = h* Function(r+h2, y1+k21/2.0, dy1+k22/2.0, omega,L);
    k41 = h* (dy1+k32);
    k42 = h* Function(r+h,  y1+k31, dy1+k32, omega,L);
    y1  = y1  + 1.0/6.0*( k11+2.0*k21+2.0*k31+k41 );
    dy1 = dy1 + 1.0/6.0*( k12+2.0*k22+2.0*k32+k42 );
	
    k11 = h* dy2;
    k12 = h* Function(r,y2,dy2,omega,L);
    k21 = h* (dy2+k12/2.0);
    k22 = h* Function(r+h2, y2+k11/2.0, dy2+k12/2.0, omega,L);
    k31 = h* (dy2+k22/2.0);
    k32 = h* Function(r+h2, y2+k21/2.0, dy2+k22/2.0, omega,L);
    k41 = h* (dy2+k32);
    k42 = h* Function(r+h,  y2+k31,  dy2+k32,  omega,L);
    y2  = y2  + 1.0/6.0*( k11+2.0*k21+2.0*k31+k41 );
    dy2 = dy2 + 1.0/6.0*( k12+2.0*k22+2.0*k32+k42 );
		
    r  = r  + h;
}
								
double greynew(double omega,double C1,double C2)
{
    double diff = RH*RH;
    double sum = 0.5*(C1*C1 + C2*C2);
    if(sum<diff)
    {
        cout<<"ERROR: sum<difference!!"<<endl;
    }
    return 1.0 - (sum-diff)/(sum+diff);
}			
								 
int main()
{
    double r, h = 1.0e-2;                 // stepsize of r
    double x;
    double omega[nw];                     // in unit of k:=sqrt(3.)/l
    int    L=0;                           // angular momentum #
    double Phi1, Phi2, dPhi1, dPhi2;
    double C1=0.,C2=0.;                   // the ratios for the two asymptos, also peak value in the oscillating region
    double greyfactor[nw];
    double evaporation[nw];
    int nwbreak = nw;
    ofstream output("ratios.txt");
    ofstream greybody("greyfactor.txt");
	
	
    // creating the list of values of omega and L	
    for(int i=0;i<nw;i++)
    {
        omega[i] = omega0 + double(i)*dw;
	greyfactor[i]  =0.;
	evaporation[i] = 0.;
    }

    // Find the value of RH, f(RH)=0
    double M    = (RH*(1 + RH))/(1 + RH + pow(RH,2));
    double N    = 1/(1 + RH + pow(RH,2));
    double kH = df(RH)/2.0;
    double beta = 2.*pi/kH; 
    cout<<"RH= "<<RH<<endl;
    cout<<"kH= "<<df(RH)/2.0<<endl;  
    cout<<"beta= "<<beta<<endl;   
    cout<<"M = "<<M<<" N = "<<N<<endl;
    
    // initial value of r
    cout.precision(50);
    
    // Loop over all values of omega    
    for(int i=0;i<nw;i++)
    {
        for(L=0;L<nL;L++)
	{
	    // define the solution near the horizon, real and imaginary parts, separately
            r = r0;
	    h = h0;		
	    Phi1  = 1.;                   // Phi_A1(r,omega[i]);
    	    Phi2  = 0.;                   // Phi_A2(r,omega[i]);
            dPhi1 = 0;                    //-RH/pow((RH+exp(z)),2.0)*exp(z);
            dPhi2 = omega[i]/kH/2./(r-RH);
            C1 = 0.;
	    C2 = 0.;
            while(1)
	    {
	        h=h0;
		rk4(Phi1, dPhi1, Phi2, dPhi2, x, h, omega[i], L);	
                // adjust the step size from that of dx
                
		if(fabs(Phi1)>C1) C1 = fabs(Phi1);
		if(fabs(Phi2)>C2) C2 = fabs(Phi2);
	    	output<<setw(10)<<x<<' '<<setw(20)<<Phi1<<' '<<setw(20)<<Phi2<<' '<<setw(20)<<dPhi1<<' '<<setw(20)<<dPhi2<<'\n'; 
	    	    
                if(r>rend) break;
            } 
            greyfactor[i]  += double(2*L+1) * greynew(omega[i],C1,C2);
	    cout<<"grey = "<<setw(10)<<greynew(omega[i],C1,C2)<<"   L= "<<L<<"   w= "<<omega[i]<<endl;
       }   
   }
	
   for(int i=0; i<nw;i++)
   {
       evaporation[i] = greyfactor[i]/2.0/pi*omega[i]/(exp(beta*omega[i])-1.);
       greybody<<setw(10)<<omega[i]<<' '<<setw(20)<<greyfactor[i]<<' '<<setw(20)<<evaporation[i]<<'\n';	
    }
    return 0;
}
