//g++ -std=c++11 -Wall -O3 structure-factor.cc -o testo
// vim: set ai et cin ts=4 sw=4 tw=80:
//drawing Kitaev honeycomb lattice & calculating energy + magnetization
// in-plane vectors are r1=(0,1,-1), r2= =(-1,1,0), r3 =(0,1,-1) as they lie on 
//the plane formed by cutting the 3 points (1,0,0), (0,1,0) and  (0,0,1)
//the perp vector is (1,1,1)/ \sqrt{3} and we choose the in-plane dirn as r2/\sqrt{2}
// giving H = hmag sin theta (1,1,1)/ \sqrt{3} + hmag cos theta (-1,1,0)/ \sqrt{2}
// tau  is obtained from cross product of h and magnetizaion vector
// unit vector perp to r2 and h is rperp = (2,2,-1)/3.0
// so we need to plot tauperp = tau . rperp



#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <array>
#include <algorithm>
#include <complex>

boost::random::mt19937 gen;
using boost::lexical_cast;
using namespace std;

const double pi = acos(-1.0);

// Define lattice constants
const int axis1 = 10, axis2 = axis1;
const int no_of_sites = 2*axis1*axis2;
const double beta_we_want=1;

typedef boost::multi_array < double, 2 > array_2d_float;
typedef boost::multi_array < int, 2 > array_2d_int;



const double k= -60.0, g = 30.0 ;
const double hmag =50.0;

////////////////////////////////////////////////////////////////////////////////





int main(int argc, char const * argv[])
{     	
	

    array_2d_float sitepos(boost::extents[no_of_sites][2]);
    array_2d_float sitespin(boost::extents[no_of_sites][3]); 

    

    string hmag_str = lexical_cast<string>(hmag);
    string k_str = lexical_cast<string>(-k);
    string g_str = lexical_cast<string>(g);
    string T_str = lexical_cast<string>(1/beta_we_want);
    ofstream f1out(string("Strfac_T_" +T_str+"_h_"+hmag_str + "_minusK_" + k_str + "_G_"
                    +g_str+".dat").c_str()); 

  
    f1out << "qx\t qy \t Sq \t  Sq-squared "<<endl;
    int lab;



   // ifstream gin(string("Spins_T_" +T_str+"_h_"+hmag_str + "_minusK_" + k_str 
              //   + "_G_"+g_str+"_theta_0.031415926535897934.dat").c_str(),ios_base::app); 
ifstream gin("Spins_T_1_h_50_minusK_60_G_30_theta_3.1415926535897931.dat");

 
    string dummyLine;
    getline(gin, dummyLine);//omit reading header line of data file

    for (unsigned int label=0; label< no_of_sites; ++label)
    {
            gin>>lab;
            gin>>sitepos[label][0];
            gin>>sitepos[label][1];
            gin>>sitespin[label][0];
            gin>>sitespin[label][1];
            gin>>sitespin[label][2];
            printf ("%d % f %f %f %f %f\n",lab,sitepos[label][0],sitepos[label][1],
                     sitespin[label][0],sitespin[label][1],sitespin[label][2]);

        
        printf (" \n");
    }


    
/////////////////////////MAIN LOOP STARTS//////////////////////////////////////
     for (unsigned int qxsteps=0; qxsteps<61; ++qxsteps)
    {
        double qx = 0 + qxsteps*2*pi/(10) ;
        

        for (unsigned int qysteps=0; qysteps<61; ++qysteps)
        {
           double sumreal=0, sumimag=0; 
           double qy = 0 + qysteps*2*pi/10 ;
           for (unsigned int i = 0; i < no_of_sites; ++i)
           {      
             sumreal=sumreal+sitespin[i][0]*cos(qx*sitepos[i][0]+qy*sitepos[i][1]);  

             sumimag=sumimag+sitespin[i][0]*sin(qx*sitepos[i][0]+qy*sitepos[i][1]); 


             sumreal=sumreal+sitespin[i][1]*cos(qx*sitepos[i][0]+qy*sitepos[i][1]);  

             sumimag=sumimag+sitespin[i][1]*sin(qx*sitepos[i][0]+qy*sitepos[i][1]); 


             sumreal=sumreal+sitespin[i][2]*cos(qx*sitepos[i][0]+qy*sitepos[i][1]);  

             sumimag=sumimag+sitespin[i][2]*sin(qx*sitepos[i][0]+qy*sitepos[i][1]);                   
             

            }
            f1out.setf( ios_base::fixed, ios_base::floatfield );
            f1out.precision(7);f1out<< setw(12)<< qx
                          << setw(12)<< qy
            << setw(12)<<sqrt(sumreal*sumreal+sumimag*sumimag)/sqrt(no_of_sites)
            << setw(12)<<(sumreal*sumreal+sumimag*sumimag)/double( no_of_sites)
            << endl;
            // printing to file

        }
 
    }
 
	f1out.close();

	
	return 0;
}

