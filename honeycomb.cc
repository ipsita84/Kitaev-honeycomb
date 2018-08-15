//g++ -std=c++11 -Wall -O3 honeycomb.cc -o testo
// vim: set ai et cin ts=4 sw=4 tw=80:
//drawing Kitaev honeycomb lattice & calculating energy + magnetization

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <math.h>
#include <array>
#include <algorithm>


// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;

const double pi = acos(-1.0);

// Define lattice constants
const int axis1 = 4, axis2 = axis1;
const int no_of_sites = 2*axis1*axis2;


typedef boost::multi_array < double, 2 > array_2d_float;
typedef boost::multi_array < int, 2 > array_2d_int;


//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d_float sitespin, array_2d_float J1,array_2d_float J2,
       array_2d_float J3,std::array <double, 2> h);





int main(int argc, char const * argv[])
{     	
	ofstream fpout("pos.dat");	// Opens a file for output

    array_2d_float sitepos(boost::extents[no_of_sites][2]);
    array_2d_float sitespin(boost::extents[3][no_of_sites]); 

    std::array <double, 2> h = {0,0};
  
	for (unsigned int j = 0; j < no_of_sites ; ++j)
	{	int alpha = (j+1)%(2*axis1);
		int beta = ceil(double(j+1)/ double(2*axis1) );
		int gamma = ceil(double(alpha)/ double(2)) ;
		int delta = alpha % 2;

		sitepos[j][0] = double(gamma)*sqrt(3.0)-double(beta)*sqrt(3.0)*cos(pi/3.0);
		sitepos[j][1] =-double(beta)*sqrt(3.0)*sin(pi/3.0)+ double(delta);

        //printf ("x %f y %f \n", sitepos[j][0], sitepos[j][1]);

        fpout.setf( ios_base::fixed, ios_base::floatfield );
        fpout.precision(7);
        fpout << setw(20) << sitepos[j][0];
        fpout.precision(7);
        fpout << setw(20)<< sitepos[j][1] << endl;
        // writing position coordinates to file "pos.dat"
	}
	fpout.close();

	for (unsigned int j = 0; j < no_of_sites ; ++j)
	{	

            double theta = roll_coin(0,pi);
            double phi = roll_coin(0,2*pi);
            sitespin[0][j] = sin(theta)*cos(phi);
            sitespin[1][j] = sin(theta)*sin(phi);
            sitespin[2][j] = cos(theta); 
            double checksum= pow(sitespin[0][j],2) + pow(sitespin[1][j],2)
                      + pow(sitespin[2][j],2);
            if (checksum > 1.000001) {printf (" initial %f error \n", checksum);}

        }

    array_2d_float J1(boost::extents[3][3]);
    array_2d_float J2(boost::extents[3][3]);
    array_2d_float J3(boost::extents[3][3]);
    double mx=0, my=0, mz=0;

    //Read the random signed bonds for a particular stored realization
    ifstream gin("J1.dat");
    ofstream f1out("mag1_hx.dat",std::fstream::app);	// Opens a file for output
    ofstream fout("Energy1_hx.dat", std::fstream::app);


    double energy(0);
    double en_sum;

    energy = energy_tot(sitespin, J1,J2,J3, h);
    printf("%f\n",energy);
 
	fout.close();
	f1out.close();
	
	return 0;
}



//function to calculate total energy
//for a given spin configuration with pbc
double energy_tot(array_2d_float sitespin, array_2d_float J1, 
                  array_2d_float J2,array_2d_float J3,
                  std::array <double, 2> h)
{
    double energy = 0;

    for (unsigned int i = 0; i < no_of_sites ; ++i)
        {
              energy += -h[0]*sitespin[0][i]-h[1]*sitespin[1][i];
            
        }

    array_2d_int A(boost::extents[axis1][axis2]);//A sublattice has even numbered sites
    array_2d_int B(boost::extents[axis1][axis2]);//B sublattice has even odd sites 


    int counter = 0;
	for (unsigned int i = 0; i < axis1 ; ++i)
	{	
        for (unsigned int j = 0; j < axis2 ; ++j)        
        {   counter = counter +1 ;
            A[i][j] = 2*counter;
            B[i][j] = A[i][j]-1;
            //printf (" % d",A[i][j]);
            //printf (" % d",B[i][j]);
        }
        //printf ("\n");
    }



    // simple rotation to the left
    array_2d_int rotateleftA(boost::extents[axis1][axis2]);
	for (unsigned int i = 0; i < axis1-1 ; ++i)
	{	
        for (unsigned int j = 0; j < axis2 ; ++j)        
        {   rotateleftA[i][j] = A[i+1][j];
            //printf (" % d",rotateleftA[i][j]);
        }
        //printf ("\n");
    }
    for (unsigned int j = 0; j < axis2 ; ++j)        
        {   rotateleftA[axis1-1][j] = A[0][j];
            //printf (" % d",rotateleftA[axis1-1][j]);
        }


   for (unsigned comp1  = 0; comp1 < 3; ++comp1)
    {
        for (unsigned comp2  = 0; comp2 < 3; ++comp2)
        {

            for (unsigned int i = 0; i < axis1 ; ++i)
            {
                for (unsigned int j = 0; j < axis2 ; ++j)
                {
                 energy += J1[comp1][comp2]*sitespin[comp1][A[i][j]-1]
                                           *sitespin[comp2][B[i][j]-1];
                 energy += J1[comp2][comp1]*sitespin[comp2][A[i][j]-1]
                                           *sitespin[comp1][B[i][j]-1];
                 energy += J2[comp1][comp2]*sitespin[comp1][rotateleftA[i][j]-1]
                                           *sitespin[comp2][B[i][j]-1];
                 energy += J2[comp2][comp1]*sitespin[comp2][rotateleftA[i][j]-1]
                                           *sitespin[comp1][B[i][j]-1];
                 
                }
            }
        }
    }
    return energy;
}




//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
    boost::random::uniform_int_distribution <> dist(a, b);
    return dist(gen);
}

//function to generate random real no.
// between 2 integers a & b, including a & excluding b

double random_real(int a, int b)
{
    boost::random::uniform_real_distribution <> dist(a, b);
    // uniform_real_distribution: continuous uniform distribution
    //on some range [min, max) of real number
    return dist(gen);
}


