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
    xtensor<double,2> E_sub = zeros<double>({4-1,3});
    xtensor<double,2> H_sub = zeros<double>({4-1,3});


    // Initialize a vector that will store the truth values for each comparison
    xtensor<bool,1> truth_vec;
    truth_vec.resize({4-1}); 

    bool isConverged = true; 

    xtensor<double,1> EleftValues = {1,1,1}; 
    xtensor<double,1> ErightValues = {2,2,2}; 
    xtensor<double,1> HleftValues = {3,3,3}; 
    xtensor<double,1> HrightValues = {4,4,4}; 

    xtensor<double,1> subdomains = ones<double>({4}); 
    // Loops through each subdomain
    for(int count =  0; count < 4 - 1; count++)
    {
        view(A_E, count, all()) = subdomains[count]*ErightValues;
        view(A_H, count, all()) = subdomains[count]*HrightValues;

        view(B_E, count, all()) = subdomains[count + 1]*EleftValues;
        view(B_H, count, all()) = subdomains[count + 1]*HleftValues;
    }

    cout << "A_E\n" << A_E << endl;
    cout << "A_H\n" << A_H << endl;
    cout << "B_E\n" << B_E << endl;
    cout << "B_H\n" << B_H << endl;

    A_E = ones<double>({3,3}); 
    
    B_E = ones<double>({3,3});

    view(B_E, 0, all()) = view(B_E, 0, all()) - numeric_limits<double>::epsilon();

    A_H = ones<double>({3,3});

    B_H = ones<double>({3,3});
    
    view(B_H, 0, all()) = view(B_H, 0, all()) - numeric_limits<double>::epsilon(); 

    cout << "A_E\n" << A_E << endl;
    cout << "A_H\n" << A_H << endl;
    cout << "B_E\n" << B_E << endl;
    cout << "B_H\n" << B_H << endl;

    for(int n =  0; n < 4-1; n++)
    {
        //E_error = linalg::norm(view(A_E, n, all()) - view(B_E, n, all(), 2)); 
        //H_error = linalg::norm(view(A_H, n, all()) - view(B_H, n, all(), 2));

        //view(E_error, n, all()) = linalg::norm(E_sub, n, all()); 
        view(E_sub, n, all()) = view(A_E, n, all()) - view(B_E, n, all()); 
        view(H_sub, n, all()) = view(A_H, n, all()) - view(B_H, n, all()); 

        E_error(n) = linalg::norm(view(E_sub, n, all()),2); 
        H_error(n) = linalg::norm(view(H_sub, n, all()),2); 

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

    cout << "E_sub\n" << E_sub << endl;
    cout << "H_sub\n" << H_sub << endl;
    cout << "E_error\n" << E_error << endl;
    cout << "H_erorr\n" << H_error << endl;
    cout << "truth vec\n" << truth_vec << endl; 
    cout << "Epsilon: " << numeric_limits<double>::epsilon() << endl;
    return isConverged;
}


int main()
{

    xtensor<double,1> A = linspace(0,10,11);
    cout << view(A,all()) << endl;
    cout << view(A,range(0,-2)) << endl;
    cout << view(A,range(1,-1)) << endl;
    cout << view(A,range(1,_)) << endl;
    view(A,range(-1,_)) = 99.9;
    cout << view(A,all()) << endl;
    // cout << view(A, 0, all()) - view(B, 0, all()) << endl;
    // cout << B << endl; 
    cout << view(A,range(0,A.shape(0)));
    cout << A.shape(0) << "  " << A.size() << endl;
    //bool isConverged; 
    //isConverged = check_convergence(); 
    //cout << check_convergence() << endl; 





    
    return 0;
}