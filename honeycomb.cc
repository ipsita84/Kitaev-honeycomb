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

// Define global scopes to use them across all functions
double J = 1.0;
const unsigned int axis1 = 3, axis2 = axis1;
const unsigned int no_of_sites = 2*axis1*axis2;
// above assigns length along each dimension of the 2d configuration

typedef boost::multi_array < double, 2 > array_2d_float;

int main()
{     	
	ofstream fout("pos.dat");	// Opens a file for output

    array_2d_float sitepos(boost::extents[no_of_sites][2]); 
  
	for (unsigned int j = 0; j < no_of_sites ; ++j)
	{	unsigned int alpha = (j+1)%(2*axis1);
		unsigned int beta = ceil(double(j+1)/ double(2*axis1) );
		unsigned int gamma = ceil(double(alpha)/ double(2)) ;
		unsigned int delta = alpha % 2;
        printf (" beta %d\n", beta);

		sitepos[j][0] = gamma*sqrt(3.0) - beta*sqrt(3.0)*cos(pi/3.0);
		sitepos[j][1] =-beta*sqrt(3.0) * sin(pi/3.0) + delta*1.0;
        //printf (" y %f\n", sin(pi/3.0) );
       // printf ("alpha %d beta %d gamma %d\n", alpha, beta, gamma);

        printf ("x %f y %f \n", sitepos[j][0], sitepos[j][1]);

        fout.setf( ios_base::fixed, ios_base::floatfield );
        fout.precision(7);
        fout << setw(20) << sitepos[j][0];
        fout.precision(7);
        fout << setw(20)<< sitepos[j][1] << endl;
        // writing position coordinates to file "pos.dat"
	}
	

	fout.close();
	
	return 0;
}











