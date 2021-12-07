/*#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>     // xt::arange
#include <xtensor/xmath.hpp>        // xt::sin, cos
#include <complex>
#include <xtensor/xio.hpp>
#include <array>
#include <vector>
#include <xtensor/xtensor.hpp>
#include "xtensor/xview.hpp"
#include "xtensor/xpad.hpp"
#include "xtensor/xmanipulation.hpp"
using namespace std;
using namespace xt;*/ 

#include <iostream>
#include <istream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include <ctime>
#include <thread>
#include <list>
#include <complex>
#include <omp.h>
#include <mutex>  

//Xtensor preprocessor directives
#include "xtensor/xarray.hpp"
#include <xtensor/xbuilder.hpp>     // xt::arange
#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xcomplex.hpp"
#include "xtensor/xnpy.hpp"
#include "xtensor/xio.hpp"
//#include "xtensor-fftw/basic.hpp"   // rfft, irfft
//#include "xtensor-fftw/helper.hpp"  // rfftscale
//#include "xtensor-io/xhighfive.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xpad.hpp"


using namespace std;
using namespace xt;
using namespace std::this_thread;
using namespace std::chrono;
using namespace std::complex_literals;

// g++ test.cpp -std=c++17 -Wall -Werror -O0 -isystem ../include -isystem ../../../../libs/xtensor/include -isystem ../../../../libs/xtl/include -isystem ../../../../libs/xtensor-io/include -isystem ../../../../libs/xtensor-fftw/include -isystem ../../../../libs/xtensor-blas/include -isystem ../../../../libs/HighFive/include -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -lfftw3 -fopenmp -lblas -llapack -DHAVE_CBLAS=1 -o test.out


bool check_convergence()
{
    // Initialize matrices of overlapping region values 
    xtensor<double,2> A_E = zeros<double>({4-1, 3}); // 4 -> n_subdoms, 3 -> overlap_size
    xtensor<double,2> B_E = zeros<double>({4-1, 3});
    xtensor<double,2> A_H = zeros<double>({4-1, 3}); 
    xtensor<double,2> B_H = zeros<double>({4-1, 3});

    // Initialize error vectors
    xtensor<double,1> E_error = zeros<double>({4-1}); 
    xtensor<double,1> H_error = zeros<double>({4-1});

    // Initialize vectors that will store A - B 
    xtensor<double,1> E_sub = zeros<double>({4-1});
    xtensor<double,1> H_sub = zeros<double>({4-1});


    // Initialize a vector that will store the truth values for each comparison
    xtensor<bool,1> truth_vec;
    truth_vec.resize({4-1}); 

    bool isConverged = true; 

    xtensor<double,1> EleftValues = {1,1,1}; 
    xtensor<double,1> ErightValues = {2,2,2}; 
    xtensor<double,1> HleftValues = {3,3,3}; 
    xtensor<double,1> HrightValues = {4,4,4}; 

    // Loops through each subdomain
    for(int count =  0; count < 4 - 1; count++)
    {
        view(A_E, count, all()) = subdomains[count].getOverlapValues('E',"right");
        view(A_H, count, all()) = subdomains[count].getOverlapValues('H',"right");

        view(B_E, count, all()) = subdomains[count + 1].getOverlapValues('E', "left"); 
        view(B_H, count, all()) = subdomains[count + 1].getOverlapValues('H', "left");
    }

    for(int n =  0; n < 4-1; n++)
    {
        //E_error = linalg::norm(view(A_E, n, all()) - view(B_E, n, all(), 2)); 
        //H_error = linalg::norm(view(A_H, n, all()) - view(B_H, n, all(), 2));

        //view(E_error, n, all()) = linalg::norm(E_sub, n, all()); 
        E_sub = view(A_E, n, all()) - view(B_E, n, all()); 
        H_sub = view(A_H, n, all()) - view(B_H, n, all()); 

        E_error(n) = linalg::norm(E_sub,2); 
        H_error(n) = linalg::norm(H_sub,2); 

        if (E_error(n) <= numeric_limits<double>::epsilon() && H_error(n) <= numeric_limits<double>::epsilon())
        {
            truth_vec(n) = true; 
        } else
        {
            truth_vec(n) = false; 
        }  

        if (truth_vec(n) == false)
        {
            isConverged = false; 
        }
    }
    return isConverged;
}


int main()
{

    xtensor<double,2> A = ones<double>({3,3}); 
    xtensor<double,2> B = ones<double>({3,3})*2; 

    cout << view(A, 0, all()) - view(B, 0, all()) << endl;
    //cout << B << endl; 




    
    return 0;
}