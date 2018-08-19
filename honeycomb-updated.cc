//g++ -std=c++11 -Wall -O3 honeycomb-updated.cc -o testo
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

boost::random::mt19937 gen;
using namespace std;

const double pi = acos(-1.0);

// Define lattice constants
const int axis1 = 10, axis2 = axis1;
const int no_of_sites = 2*axis1*axis2;


typedef boost::multi_array < double, 2 > array_2d_float;
typedef boost::multi_array < int, 2 > array_2d_int;


//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);

double energy_tot(array_2d_float sitespin, array_2d_float J1,array_2d_float J2,
       array_2d_float J3,std::array <double, 2> h, array_2d_int A, 
       array_2d_int B, array_2d_int rotateleftA, array_2d_int cornerA);

double nn_energy(array_2d_float sitespin, array_2d_float J1, 
                  array_2d_float J2,array_2d_float J3,
                  std::array <double, 2> h, array_2d_int A, array_2d_int B, 
array_2d_int rotateleftA, array_2d_int cornerA,array_2d_int rotaterightB, 
array_2d_int cornerB, int row, int col, int sublat);


//No.of Monte Carlo updates we want
const unsigned int N_mc = 1e5;

const double beta=0.1;
//const double K= -60, G = 30 ;


////////////////////////////////////////////////////////////////////////////////





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


    array_2d_int A(boost::extents[axis1][axis2]);//A sublattice has even numbered sites
    array_2d_int B(boost::extents[axis1][axis2]);//B sublattice has even odd sites 


    int counter = 0;
	for (unsigned int i = 0; i < axis1 ; ++i)
	{	
        for (unsigned int j = 0; j < axis2 ; ++j)        
        {   counter = counter +1 ;
            B[i][j] = 2*counter;
            A[i][j] = B[i][j]-1;
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

    array_2d_int cornerA(boost::extents[axis1][axis2]);
	for (unsigned int i = 0; i < axis1 ; ++i)
	{	
        for (unsigned int j = 0; j < axis2-1 ; ++j)        
        {   cornerA[i][j] = rotateleftA[i][j+1];
            //printf (" % d",cornerA[i][j]);
        }
        cornerA[i][axis2-1]=rotateleftA[i][0];
        //printf (" % d",cornerA[i][axis2-1]);
        //printf ("\n");
    }

    array_2d_int rotaterightB(boost::extents[axis1][axis2]);
    for (unsigned int j = 0; j < axis2 ; ++j)        
        {   rotaterightB[0][j] = B[axis1-1][j];
            //printf (" % d",rotaterightB[0][j]);
        }
    //printf ("\n");
	for (unsigned int i = 1; i < axis1 ; ++i)
	{	
        for (unsigned int j = 0; j < axis2 ; ++j)        
        {   rotaterightB[i][j] = B[i-1][j];
            //printf (" % d",rotaterightB[i][j]);
        }
        //printf ("\n");
    }


    array_2d_int cornerB(boost::extents[axis1][axis2]);
	for (unsigned int i = 0; i < axis1 ; ++i)
	{	
         cornerB[i][0]=rotaterightB[i][axis2-1];
         //printf (" % d",cornerB[i][0]);
        for (unsigned int j = 1; j < axis2 ; ++j)        
        {   cornerB[i][j] = rotaterightB[i][j-1];
            //printf (" % d",cornerB[i][j]);
        }
        //printf ("\n");
    }
    

    array_2d_float J1(boost::extents[3][3]);
    array_2d_float J2(boost::extents[3][3]);
    array_2d_float J3(boost::extents[3][3]);
    double mplanar=0, mperp=0;

    //Read the random signed bonds for a particular stored realization

    ofstream f1out("mag.dat",std::fstream::app);	// Opens a file for output
    ofstream fout("Energy.dat", std::fstream::app);

    ifstream gin("J1.dat");
    for (unsigned int comp1=0; comp1<3; ++comp1)
    {
        for (unsigned int comp2=0; comp2<3; ++comp2)
        {

            gin>>J1[comp1][comp2];
            printf (" % f\n",J1[comp1][comp2]);

        }
    }
    for (unsigned int comp1=0; comp1<3; ++comp1)
    {
        for (unsigned int comp2=0; comp2<3; ++comp2)
        {

            gin>>J2[comp1][comp2];
            printf (" % f\n",J2[comp1][comp2]);

        }
    }
    for (unsigned int comp1=0; comp1<3; ++comp1)
    {
        for (unsigned int comp2=0; comp2<3; ++comp2)
        {

            gin>>J3[comp1][comp2];
            printf (" % f\n",J2[comp1][comp2]);

        }
    }
    gin.close();


    double energy(0);
    double en_sum;

     for (unsigned int thetasteps=0; thetasteps<101; ++thetasteps)
    {
        double theta = 0 + thetasteps * pi/100 ;
        h[0] = 30.0*cos(theta);
        h[1] = 30.0*sin(theta);
        energy = energy_tot(sitespin,J1,J2,J3,h,A,B,rotateleftA,cornerA);
        en_sum =0;
        std::array <double, N_mc> energy_array =  {0};
        std::array <double, N_mc> m_planar_array={0}, m_perp_array ={0};
        unsigned int heating = 1e5;

        for (unsigned int i = 1; i <=heating + N_mc; ++i)
        {
            for (unsigned int j = 1; j <= no_of_sites; ++j)
            {
            //Now choose a random spin site at (row,col) & sublattice 0 or 1
                
                int row = roll_coin(0, axis1-1);
                int col = roll_coin(0, axis2-1);
                int sublat =roll_coin(0,1);
                  
                int label;
                if (sublat == 0) label = A[row][col]-1;
                else label = B[row][col]-1;

                double s0 = sitespin[0][label];
                double s1 = sitespin[1][label];
                double s2 = sitespin[2][label];
                double energy_old =energy ;
                double energy_minus_rnd_site =energy_old -nn_energy(sitespin,J1,
                J2,J3,h,A,B,rotateleftA,A,rotaterightB,cornerB,row,col,sublat);
  
                double r0 = 0.5*random_real(-1, 1)/beta;
                double r1 = 0.5*random_real(-1, 1)/beta;
                double r2 = 0.5*random_real(-1, 1)/beta;
                double tot = pow( s0+r0, 2)+pow( s1+ r1, 2)+pow(s2 + r2, 2);
                //printf ("tot %f \n",tot);

                sitespin[0][label] = (s0+r0)/sqrt(tot);
                sitespin[1][label] = (s1+r1)/sqrt(tot);
                sitespin[2][label]= (s2+r2)/sqrt(tot);
// algp from pg-24 of
// https://journals-aps-org.proxy.library.cornell.edu/prb/pdf/10.1103/PhysRevB.24.1391

               double checksum= pow(sitespin[0][label],2)
                               +pow(sitespin[1][label],2)
                               +pow(sitespin[2][label],2);
               if (checksum > 1.00001) {printf ("%f error \n", checksum);}

               double energy_new = energy_minus_rnd_site+nn_energy(sitespin,J1,
                J2,J3,h,A,B,rotateleftA,A,rotaterightB,cornerB,row,col,sublat);

               double energy_diff = energy_new - energy_old;
               double acc_ratio = exp(-1.0 * energy_diff* beta/50);
               double r =  random_real(0, 1) ;	//Generate a random no. r such that 0 < r < 1
                //Spin flipped if r <= acceptance ratio
                if (r <= acc_ratio)
                {
                    energy = energy_new ;
    
                }
                if (r > acc_ratio)
                {
                    sitespin[0][label] = s0;
                    sitespin[1][label] = s1;
                    sitespin[2][label] = s2;
                    energy = energy_old ;

                }
 
            }

            if (i >  heating )
            {
                en_sum += energy;
                energy_array[i- heating  -1] = energy;
 
                for (unsigned int l = 0; l < no_of_sites; ++l)
                {
                  mplanar += sitespin[0][l] ;
                  m_planar_array[i- heating  -1] += sitespin[0][l] ;
 
                  mperp += sitespin[1][l] ;
                  m_perp_array[i- heating  -1] += sitespin[1][l];
                }

            }
        }


        double sigma_en = 0, sigma_mplanar = 0, sigma_mperp = 0;
        for (unsigned i=0; i< N_mc; i++)
        {
            sigma_en += (energy_array[i] - en_sum/ N_mc) 
                        * (energy_array[i] - en_sum/ N_mc);
            sigma_mplanar += (m_planar_array[i] - mplanar/ N_mc) 
                               * (m_planar_array[i] - mplanar/ N_mc) ;

            sigma_mperp += (m_perp_array[i] - mperp/ N_mc)
                           * (m_perp_array[i] - mperp/ N_mc) ;
        }

        fout.setf( ios_base::fixed, ios_base::floatfield );
        fout.precision(2);
        fout << setw(6) << theta;
        fout.precision(7);
        fout << setw(15)
             << en_sum / N_mc << setw(15)
             << sqrt(sigma_en) / N_mc << endl;
        // printing energy to file "Energy.dat"

        f1out.setf( ios_base::fixed, ios_base::floatfield );
        f1out.precision(2);
        f1out << setw(6) << theta;
        f1out.precision(7);
        f1out << setw(15)
              <<(mplanar*sin(theta)-mperp*cos(theta))/(no_of_sites*N_mc)
              << setw(15) << mplanar/(no_of_sites*N_mc)
              << setw(15) << sqrt(sigma_mplanar)/(no_of_sites*N_mc)
              << setw(15) << mperp/(no_of_sites*N_mc)
              << setw(15) << sqrt(sigma_mperp)/(no_of_sites*N_mc)  << endl;
// printing magnetization to file "mag.dat"
//  (hz mx -  hz my)/h
        mplanar=0;
        mperp=0;
    }
 
	fout.close();
	f1out.close();
	
	return 0;
}



//function to calculate total energy
//for a given spin configuration with pbc
double energy_tot(array_2d_float sitespin, array_2d_float J1, array_2d_float J2,
array_2d_float J3, std::array <double, 2> h, array_2d_int A, array_2d_int B, 
array_2d_int rotateleftA, array_2d_int cornerA)
{
    double energy = 0;

    for (unsigned int i = 0; i < no_of_sites ; ++i)
        {
              energy += -h[0]*sitespin[0][i]-h[1]*sitespin[1][i];
            
        }



   for (unsigned comp1  = 0; comp1 < 3; ++comp1)
    {
        for (unsigned comp2  = comp1; comp2 < 3; ++comp2)
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
                 energy += J3[comp1][comp2]*sitespin[comp1][cornerA[i][j]-1]
                                           *sitespin[comp2][B[i][j]-1];
                 energy += J3[comp2][comp1]*sitespin[comp2][cornerA[i][j]-1]
                                           *sitespin[comp1][B[i][j]-1];
                 
                }
            }
        }
    }
    return energy;
}



////////////////////////////////////////////////////////////////////
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


double nn_energy(array_2d_float sitespin, array_2d_float J1, 
                  array_2d_float J2,array_2d_float J3,
                  std::array <double, 2> h, array_2d_int A, array_2d_int B, 
array_2d_int rotateleftA, array_2d_int cornerA,array_2d_int rotaterightB, 
array_2d_int cornerB, int row, int col, int sublat)
{
    double nenergy = 0;
    int i=row; int j=col;  
 
    if (sublat == 0)         
    { nenergy += -h[0]*sitespin[0][A[i][j]-1]-h[1]*sitespin[1][A[i][j]-1];
    for (unsigned comp1  = 0; comp1 < 3; ++comp1)
        {
            for (unsigned comp2  = comp1; comp2 < 3; ++comp2)
             {

                nenergy += J1[comp1][comp2]*sitespin[comp1][A[i][j]-1]
                                           *sitespin[comp2][B[i][j]-1];
                nenergy += J1[comp2][comp1]*sitespin[comp2][A[i][j]-1]
                                           *sitespin[comp1][B[i][j]-1];
                nenergy += J2[comp1][comp2]*sitespin[comp1][A[i][j]-1]
                                           *sitespin[comp2][rotaterightB[i][j]-1];
                nenergy += J2[comp2][comp1]*sitespin[comp2][A[i][j]-1]
                                           *sitespin[comp1][rotaterightB[i][j]-1];
                nenergy += J3[comp1][comp2]*sitespin[comp1][A[i][j]-1]
                                           *sitespin[comp2][cornerB[i][j]-1];
                nenergy += J3[comp2][comp1]*sitespin[comp2][A[i][j]-1]
                                           *sitespin[comp1][cornerB[i][j]-1];
                 
              }
        }
    }


    else 
    {nenergy += -h[0]*sitespin[0][B[i][j]-1]-h[1]*sitespin[1][B[i][j]-1];
            

    for (unsigned comp1  = 0; comp1 < 3; ++comp1)
        {
            for (unsigned comp2  = comp1; comp2 < 3; ++comp2)
             {

                nenergy += J1[comp1][comp2]*sitespin[comp1][A[i][j]-1]
                                           *sitespin[comp2][B[i][j]-1];
                nenergy += J1[comp2][comp1]*sitespin[comp2][A[i][j]-1]
                                           *sitespin[comp1][B[i][j]-1];
                nenergy += J2[comp1][comp2]*sitespin[comp1][rotateleftA[i][j]-1]
                                           *sitespin[comp2][B[i][j]-1];
                nenergy += J2[comp2][comp1]*sitespin[comp2][rotateleftA[i][j]-1]
                                           *sitespin[comp1][B[i][j]-1];
                nenergy += J3[comp1][comp2]*sitespin[comp1][cornerA[i][j]-1]
                                           *sitespin[comp2][B[i][j]-1];
                nenergy += J3[comp2][comp1]*sitespin[comp2][cornerA[i][j]-1]
                                           *sitespin[comp1][B[i][j]-1];
                 
              }
        }
    }

    return nenergy;
}




