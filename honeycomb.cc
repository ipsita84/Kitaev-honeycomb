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

int main()
{     	
	ofstream fout("pos.dat");	// Opens a file for output

    array_2d_float sitepos(boost::extents[no_of_sites][2]); 
    array_2d_int A(boost::extents[axis1][axis2]);//A sublattice has even numbered sites
    array_2d_int B(boost::extents[axis1][axis2]);//B sublattice has even odd sites 
  
	for (unsigned int j = 0; j < no_of_sites ; ++j)
	{	int alpha = (j+1)%(2*axis1);
		int beta = ceil(double(j+1)/ double(2*axis1) );
		int gamma = ceil(double(alpha)/ double(2)) ;
		int delta = alpha % 2;

		sitepos[j][0] = double(gamma)*sqrt(3.0)-double(beta)*sqrt(3.0)*cos(pi/3.0);
		sitepos[j][1] =-double(beta)*sqrt(3.0)*sin(pi/3.0)+ double(delta);

        //printf ("x %f y %f \n", sitepos[j][0], sitepos[j][1]);

        fout.setf( ios_base::fixed, ios_base::floatfield );
        fout.precision(7);
        fout << setw(20) << sitepos[j][0];
        fout.precision(7);
        fout << setw(20)<< sitepos[j][1] << endl;
        // writing position coordinates to file "pos.dat"
	}
    int counter = 0;
	for (unsigned int j = 0; j < axis2 ; ++j)
	{	
        for (unsigned int i = 0; i < axis1 ; ++i)        
        {   counter = counter +1 ;
            A[i][j] = 2*counter;
            printf (" % d",A[i][j]);
        }
        printf ("\n");
	}

	fout.close();
	
	return 0;
}











