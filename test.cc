//g++ -std=c++11 -Wall -O3 honeycomb_cool-once_spins_recorded.cc -o testo
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

boost::random::mt19937 gen;
using boost::lexical_cast;
using namespace std;

const double pi = acos(-1.0);
const double beta_we_want =1;

// Define lattice constants
const int axis1 = 10, axis2 = axis1;
const int no_of_sites = 2*axis1*axis2;


typedef boost::multi_array < double, 2 > array_2d_float;
typedef boost::multi_array < int, 2 > array_2d_int;


//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);

double energy_tot(array_2d_float sitespin, array_2d_float J1,array_2d_float J2,
       array_2d_float J3,std::array <double, 3> h, array_2d_int A, 
       array_2d_int B, array_2d_int rotateleftA, array_2d_int cornerA);

double nn_energy(array_2d_float sitespin, array_2d_float J1, 
                  array_2d_float J2,array_2d_float J3,
                  std::array <double, 3> h, array_2d_int A, array_2d_int B, 
array_2d_int rotateleftA, array_2d_int cornerA,array_2d_int rotaterightB, 
array_2d_int cornerB, int row, int col, int sublat);


//No.of Monte Carlo updates we want
const unsigned int nmore=1;
const unsigned int N_mc = 1e5;
const unsigned int heating = N_mc;


const double k= -60.0, g = 30.0 ;
const double hmag =50.0;

////////////////////////////////////////////////////////////////////////////////





int main(int argc, char const * argv[])
{     	
	
    ifstream gin("J.dat");
    for (unsigned int comp1=0; comp1<3; ++comp1)
    {
        for (unsigned int comp2=0; comp2<3; ++comp2)
        {

            gin>>J1[comp1][comp2];
            printf (" % f\n",J1[comp1][comp2]);

        }
    }

ofstream fpout("pos.dat");	// Opens a file for output

    array_2d_float sitepos(boost::extents[no_of_sites][2]);
    array_2d_float sitespin(boost::[no_of_sites]extents[3]); 

    

    string hmag_str = lexical_cast<string>(hmag);
    string k_str = lexical_cast<string>(-k);
    string g_str = lexical_cast<string>(g);
    string T_str = lexical_cast<string>(1/beta_we_want);
    ofstream f1out(string("Strfac_T_" +T_str+"_h_"+hmag_str + "_minusK_" + k_str + "_G_"
                    +g_str+".dat").c_str(),ios_base::app); 

  
    f1out << "qx\t qy \t Sq \t  Sq-squared "<<endl;


/////////////////////////MAIN LOOP STARTS//////////////////////////////////////
     for (unsigned int qxsteps=0; qxsteps<101; ++qxsteps)
    {
        double qx = 0 + qxsteps*2*pi/100 ;
        

        for (unsigned int qysteps=0; qysteps<101; ++qysteps)
        {
            

                for (unsigned int l = 0; l < no_of_sites; ++l)
                {
                  mx += sitespin[0][l] ;
                  mx_array[i -1] += sitespin[0][l] ;
 
                  my += sitespin[1][l] ;
                  my_array[i -1] += sitespin[1][l];

                  mz += sitespin[2][l] ;
                  mz_array[i -1] += sitespin[2][l];
                  sout << setw(3) << l+1;
                  sout.setf( ios_base::fixed, ios_base::floatfield );
                  sout.precision(7);sout << setw(12)<< sitepos[l][0]
                  <<setw(12)<< sitepos[l][1]
                  << setw(12)<< sitespin[0][l]
                  << setw(12)<< sitespin[1][l]
                  << setw(12)<< sitespin[2][l]<< endl;
                }
             

       
        fout.setf( ios_base::fixed, ios_base::floatfield );
        fout.precision(2);fout << setw(6) << theta;
        fout.precision(7);fout << setw(25)<< en_sum / (no_of_sites*N_mc) 
                          << setw(15)<< sqrt(sigma_en) / (no_of_sites*N_mc)
                          << setw(12)<< sqrt(sigma_en) / (no_of_sites*N_mc)
                          << setw(12)<< mx / (no_of_sites*N_mc)
                          << setw(12)<< sqrt(sigma_mx)/(no_of_sites*N_mc) 
                          << setw(12)<< my / (no_of_sites*N_mc)
                          << setw(12)<< sqrt(sigma_my)/(no_of_sites*N_mc)
                          << setw(12)<< mz /(no_of_sites*N_mc)   
                          << setw(12)<< sqrt(sigma_mz)/(no_of_sites*N_mc)
                          << endl;
        // printing energy and mag to file

 


 
        mx=0; my=0; mz=0;
    }
 

 
	f1out.close();

	
	return 0;
}

