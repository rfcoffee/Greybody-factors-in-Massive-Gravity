#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<stdio.h>


using namespace std;
const double pi = 3.1415926;
int choose = 1;
double RH = 0.01;       // M is roughly 10
double suppression = 1.; 
double h0 = 1.e-4; 
double r0 = RH*(1.0 + 1.e-3);
double z0 = -1000.;
double zend = 1000.;
double omega0 = 160;
double amplification = 1.; 
const int nw = 10;
double dw     = 0.25;
int    L0 = 3;
const int nL = 1;
double coupling = 1300000.;
char ratioout[] = "ratio.txt";
char greyout[]  = "grey.txt";


// Additional potential term due to coupling with the background Stuckelberg fields
double VMG(double r)
{
   if((-0.013266190475728622 + 0.30995964402877646*r - 2.4687768979404447*pow(r,2) + 11.152088975826482*pow(r,3) - 
   27.202479599173017*pow(r,4) + 37.65296325475077*pow(r,5) - 27.46070216805714*pow(r,6) + 8.030549198146632*pow(r,7))<0.)
   {return 0.;}
   return -0.013266190475728622 + 0.30995964402877646*r - 2.4687768979404447*pow(r,2) + 11.152088975826482*pow(r,3) - 
   27.202479599173017*pow(r,4) + 37.65296325475077*pow(r,5) - 27.46070216805714*pow(r,6) + 8.030549198146632*pow(r,7);	  
}

// First define the preliminary functions.
double f(double r)
{
    return -1/(1 + RH*(RH + 1))/r*(r - RH)*(r - 1)*(r + RH + 1);
}

double dr2f(double r) // calculated using Mathematica

    return (-2*(-1 + r)*(r - RH)*(1 + r + RH))/(1 + RH*(1 + RH)) + 
   pow(r,2)*(-(((-1 + r)*(r - RH))/(r*(1 + RH*(1 + RH)))) - 
      ((-1 + r)*(1 + r + RH))/(r*(1 + RH*(1 + RH))) + 
      ((-1 + r)*(r - RH)*(1 + r + RH))/
       (pow(r,2)*(1 + RH*(1 + RH))) - 
      ((r - RH)*(1 + r + RH))/(r*(1 + RH*(1 + RH))));
}

double df(double r)
{
    return -(((-1 + r)*(r - RH))/(r*(1 + RH*(1 + RH)))) - 
   ((-1 + r)*(1 + r + RH))/(r*(1 + RH*(1 + RH))) + 
   ((-1 + r)*(r - RH)*(1 + r + RH))/
    (pow(r,2)*(1 + RH*(1 + RH))) - 
   ((r - RH)*(1 + r + RH))/(r*(1 + RH*(1 + RH)));
}

double kH()
{
    return df(RH);
}	

double d2f(double r)
{
    return (-2*(-1 + r))/(r*(1 + RH*(1 + RH))) + 
   (2*(-1 + r)*(r - RH))/(pow(r,2)*(1 + RH*(1 + RH))) - 
   (2*(r - RH))/(r*(1 + RH*(1 + RH))) + 
   (2*(-1 + r)*(1 + r + RH))/(pow(r,2)*(1 + RH*(1 + RH))) - 
   (2*(1 + r + RH))/(r*(1 + RH*(1 + RH))) - 
   (2*(-1 + r)*(r - RH)*(1 + r + RH))/
    (pow(r,3)*(1 + RH*(1 + RH))) + 
   (2*(r - RH)*(1 + r + RH))/(pow(r,2)*(1 + RH*(1 + RH)));
}


// Total potential in the Regge-Wheeler equation
double V(double r, int L)
{
    return f(r)*df(r)/r + double(L*(L+1))*f(r)/r/r + coupling*VMG(r);	
}


/* choose the coordinate variable */
double Function(double x, double yy, double dy, double omega, int L)
{
    if(choose == 0)
    {
         double r = x;
         return -yy*(pow(omega/f(r),2.0) - double(L*(L+1))/(pow(r,2.0)*f(r)) - coupling*VMG(r)) - dy*(dr2f(r)/(pow(r,2.0)*f(r))); //massive gravity term
    }
    if(choose == 1)  // near the horizon in z-coord
    {
	     double r = exp(
            -((-1 + RH)*(1 + 2*RH)*
            (x + ((-1 - 2*RH)*(1 + RH + pow(RH,2))*log(1 - RH))/
               ((-1 + RH)*(2 + RH)*(1 + 2*RH)) - 
              ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*
                 log(1 + 2*RH))/((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
         (RH*(1 + RH + pow(RH,2)))) + RH;
         double drdx = 	-(((-1 + RH)*(1 + 2*RH))/
           (exp(((-1 + RH)*(1 + 2*RH)*
           (x + ((-1 - 2*RH)*(1 + RH + pow(RH,2))*log(1 - RH))/
              ((-1 + RH)*(2 + RH)*(1 + 2*RH)) - 
             ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*
                log(1 + 2*RH))/((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
         (RH*(1 + RH + pow(RH,2))))*RH*(1 + RH + pow(RH,2))));
         double d2rdx2 = (pow(-1 + RH,2)*pow(1 + 2*RH,2))/
         (exp(((-1 + RH)*(1 + 2*RH)*
         (x + ((-1 - 2*RH)*(1 + RH + pow(RH,2))*log(1 - RH))/
            ((-1 + RH)*(2 + RH)*(1 + 2*RH)) - 
           ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*
              log(1 + 2*RH))/((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
         (RH*(1 + RH + pow(RH,2))))*pow(RH,2)*
         pow(1 + RH + pow(RH,2),2));

         return -2./r*drdx*dy - yy*(1./r*d2rdx2 + omega*omega - V(r,L));
	}
	if(choose == 2)
	{
	    double r = 1 - 
            exp(-((-1 + RH)*(2 + RH)*(1 + 2*RH)*
            (x + (RH*(1 + RH + pow(RH,2))*log(1 - RH))/
               ((-1 + RH)*(1 + 2*RH)) - 
              ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*
                 log(2 + RH))/((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
            ((-1 - 2*RH)*(1 + RH + pow(RH,2))));
        double drdx = ((-1 + RH)*(2 + RH)*(1 + 2*RH))/
         (exp(((-1 + RH)*(2 + RH)*(1 + 2*RH)*
         (x + (RH*(1 + RH + pow(RH,2))*log(1 - RH))/
            ((-1 + RH)*(1 + 2*RH)) - 
           ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*log(2 + RH))/
            ((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
          ((-1 - 2*RH)*(1 + RH + pow(RH,2))))*(-1 - 2*RH)*
          (1 + RH + pow(RH,2)));
        double d2rdx2 = -((pow(-1 + RH,2)*pow(2 + RH,2)*pow(1 + 2*RH,2))/
           (exp(((-1 + RH)*(2 + RH)*(1 + 2*RH)*
           (x + (RH*(1 + RH + pow(RH,2))*log(1 - RH))/
              ((-1 + RH)*(1 + 2*RH)) - 
             ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*
                log(2 + RH))/((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
           ((-1 - 2*RH)*(1 + RH + pow(RH,2))))*pow(-1 - 2*RH,2)*
           pow(1 + RH + pow(RH,2),2)));
         
         return -2./r*drdx*dy - yy*(1./r*d2rdx2 + omega*omega - V(r,L));
	}
} 


// define the function to perform the 4th order Runge Kutta algorithm to update the values of phi and dphi
void rk4(double& y1, double& dy1, double& y2, double& dy2, double& r, double h, double omega, int L)
{
	double h2 = h/2.0;
	double k11,k12,k21,k22,k31,k32,k41,k42;
	k11 = h* dy1;
	k12 = h* Function(r,y1,dy1,omega,L);
    k21 = h* (dy1+k12/2.0);
	k22 = h* Function(r+h2, y1+k11/2.0, dy1+k12/2.0, omega,L);
	k31 = h* (dy1+k22/2.0);
	k32 = h* Function(r+h2, y1+k21/2.0, dy1+k22/2.0, omega,L);
	k41 = h* (dy1+k32);
	k42 = h* Function(r+h,  y1+k31,     dy1+k32,     omega,L);
	y1  = y1  + 1.0/6.0*( k11+2.0*k21+2.0*k31+k41 );
	dy1 = dy1 + 1.0/6.0*( k12+2.0*k22+2.0*k32+k42 );
	
	k11 = h* dy2;
	k12 = h* Function(r,y2,dy2,omega,L);
    k21 = h* (dy2+k12/2.0);
	k22 = h* Function(r+h2, y2+k11/2.0, dy2+k12/2.0, omega,L);
	k31 = h* (dy2+k22/2.0);
	k32 = h* Function(r+h2, y2+k21/2.0, dy2+k22/2.0, omega,L);
	k41 = h* (dy2+k32);
	k42 = h* Function(r+h,  y2+k31,     dy2+k32,     omega,L);
	y2  = y2  + 1.0/6.0*( k11+2.0*k21+2.0*k31+k41 );
	dy2 = dy2 + 1.0/6.0*( k12+2.0*k22+2.0*k32+k42 );
		
	r  = r  + h;
}


double grey(double C1,double C2)
{
    double diff = RH*RH;
    double sum = 0.5*(C1*C1+C2*C2);
    return 1.- (sum-diff)/(sum+diff);
}
								
															 
int main()
{
	double r, h = 1.0e-2;  // stepsize of r
	double z;
	double hh = h/10., hhh = h/100., hhhh=h/1000.;
	double omega[nw];   // in unit of k:=sqrt(3.)/l
	int    L=0;           // angular momentum #
	double Phi1, Phi2, dPhi1, dPhi2;
	double C1=0.,C2=0.; // the ratios for the two asymptos, also peak value in the oscillating region
	double greyfactor[nw];
	double evaporation[nw];
	int nwbreak = nw;
	ofstream output(ratioout);
	ofstream greybody(greyout);
	
// creating the list of values of omega and L	
	for(int i=0;i<nw;i++){
		omega[i] = omega0 + double(i)*dw;
		greyfactor[i]  =0.;
		evaporation[i] = 0.;
	}

// Find the value of RH, f(RH)=0
    double beta = 2.*pi/kH(); 
    cout<<"RH= "<<RH<<endl;
    cout<<"kH= "<<kH()<<endl;  
    cout<<"beta= "<<beta<<endl;
       
// initial value of r
    cout.precision(50);
    
// Loop over all values of omega    
	for(int i=0;i<nw;i++)
	{
		for(L=L0;L<(nL+L0);L++)
		{
		    // define the solution near the horizon, real and imaginary parts, separately
		    z = z0;
	    	h = h0;		
	    	Phi1  = amplification*1.;
    	    Phi2  = 0.;
            dPhi1 = 0; 
            dPhi2 = amplification*omega[i]; 
            choose = 1;
            C1 = 0;
            C2 = 0;
            h = h0;
            while(1)
	    	{
				if(z > zend) break;
				if(((z/RH) > -10.) && (choose==1))
				{
					double dx1dr = 1./ (-(((-1 + RH)*(1 + 2*RH))/
                    (exp(((-1 + RH)*(1 + 2*RH)*
                    (z + ((-1 - 2*RH)*(1 + RH + pow(RH,2))*log(1 - RH))/
                    ((-1 + RH)*(2 + RH)*(1 + 2*RH)) - 
                    ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*
                    log(1 + 2*RH))/((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
                   (RH*(1 + RH + pow(RH,2))))*RH*(1 + RH + pow(RH,2)))));
         
				    choose = 0;	
				    dPhi1 = dx1dr*dPhi1;
				    dPhi2 = dx1dr*dPhi2;
				    
				    z = exp(-((-1 + RH)*(1 + 2*RH)*
                       (z + ((-1 - 2*RH)*(1 + RH + pow(RH,2))*log(1 - RH))/
                       ((-1 + RH)*(2 + RH)*(1 + 2*RH)) - 
                       ((-1 + pow(RH,2))*(1 + RH + pow(RH,2))*log(1 + 2*RH))/
                       ((-1 + RH)*(2 + RH)*(1 + 2*RH))))/
                       (RH*(1 + RH + pow(RH,2)))) + RH;
     
				    h = h0*1.e-4;
				}
				if((z>1.0-5.e-5) && (choose==0)) 
				{
					double drdx2 = 1./(((-1 - 2*RH)*(1 + RH + pow(RH,2)))/
                           ((1 - z)*(-1 + RH)*(2 + RH)*(1 + 2*RH)));
                           
				    choose = 2;
				    dPhi1 = drdx2*dPhi1;
				    dPhi2 = drdx2*dPhi2;	
				    
				     z = -(((1 + RH + pow(RH,2))*((-1 - 2*RH)*log(1 - z) + 
                        RH*(2 + RH)*log(1 - RH) - (-1 + pow(RH,2))*log(2 + RH)))/
                        ((-1 + RH)*(2 + RH)*(1 + 2*RH)));
				    
				    h = h0;
				}
				
				rk4(Phi1, dPhi1, Phi2, dPhi2, z, h, omega[i], L);	
                // adjust the step size from that of dx
		    	
		    	if((choose == 2) && (fabs(Phi1)>C1)) C1 = Phi1;
		    	if((choose ==2) && (fabs(Phi2)>C2)) C2 = Phi2; /
		     }
	    
            greyfactor[i]  += double(2*L+1) * grey(C1,C2);
		    cout<<"grey = "<<setw(10)<<grey(C1,C2)<<"   L= "<<L<<"   w= "<<omega[i]<<endl;
         }   
	}
	
    for(int i=0; i<nw;i++)
    {
		evaporation[i] = greyfactor[i]/2.0/pi*pow(omega[i],1.0)/(exp(beta*omega[i])-1.);
	    greybody<<setw(10)<<omega[i]<<' '<<setw(20)<<greyfactor[i]<<' '<<setw(20)<<evaporation[i]<<'\n';	
	}
	cout<<"T_H: "<<df(RH)/(4*pi)<<endl;
	return 0;
}
