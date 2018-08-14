//g++ -std=c++11 -Wall -O3 honeycomb.cc -o testo
// vim: set ai et cin ts=4 sw=4 tw=80:
//drawing the honeycomb lattice

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


// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;

const double pi = acos(-1.0);
//const double pi = 3.14;

// Define lattice constants
double J = 1.0;
const int axis1 = 10, axis2 = axis1;
const int no_of_sites = 2*axis1*axis2;


typedef boost::multi_array < double, 2 > array_2d_float;
typedef boost::multi_array < int, 2 > array_2d_int;


//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d_float sitespin, array_2d_float J, std::array <double, 2> h);
//double nn_energy(array_2d_float sitespin,  array_2d J, std::array <double, 2> h,
      // unsigned int row, unsigned int col);




int main(int argc, char const * argv[])
{     	
	ofstream fpout("pos.dat");	// Opens a file for output

    array_2d_float sitepos(boost::extents[no_of_sites][2]);
    array_2d_float sitespin(boost::extents[3][no_of_sites]); 
    array_2d_int A(boost::extents[axis1][axis2]);//A sublattice has even numbered sites
    array_2d_int B(boost::extents[axis1][axis2]);//B sublattice has even odd sites 

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

    int counter = 0;
	for (unsigned int j = 0; j < axis2 ; ++j)
	{	
        for (unsigned int i = 0; i < axis1 ; ++i)        
        {   counter = counter +1 ;
            A[i][j] = 2*counter;
            B[i][j] = A[i][j]-1;
            //printf (" % d",A[i][j]);
            //printf (" % d",B[i][j]);

            double theta = roll_coin(0,pi);
            double phi = roll_coin(0,2*pi);
            sitespin[0][A[i][j]-1] = sin(theta)*cos(phi);
            sitespin[1][A[i][j]-1] = sin(theta)*sin(phi);
            sitespin[2][A[i][j]-1] = cos(theta); 

            theta = roll_coin(0,pi);
            phi = roll_coin(0,2*pi);
            sitespin[0][B[i][j]-1] = sin(theta)*cos(phi);
            sitespin[1][B[i][j]-1] = sin(theta)*sin(phi);
            sitespin[2][B[i][j]-1] = cos(theta); 

        }
        //printf ("\n");
	}

    array_2d_float J(boost::extents[3][3]);
    double mx=0, my=0, mz=0;

    //Read the random signed bonds for a particular stored realization
    ifstream gin("J1.dat");
    ofstream f1out("mag1_hx.dat",std::fstream::app);	// Opens a file for output
    ofstream fout("Energy1_hx.dat", std::fstream::app);



	
	return 0;
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


