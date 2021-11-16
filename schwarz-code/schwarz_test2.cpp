#ifndef FDTD
#define FDTD
 
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

//Xtensor preprocessor directives
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xcomplex.hpp"
#include "xtensor/xnpy.hpp"
#include "xtensor/xrandom.hpp"

using namespace std;
using namespace xt;
using namespace std::this_thread;
using namespace std::chrono;
using namespace std::complex_literals;

//g++ -std=c++17 schwarz_test2.cpp -I ../../../xtensor/include -I ../../../xtl/include -o schwarz_test2.out -o schwarz_test2 && ./schwarz_test2

// Implementing the array of subdomains as a class 
class subdomArray{
    public: 
        int nSubDom; 
        int lSubDom;
        int lOlap;  
        int lNolap; 

        // Function for creating 2D array of subdomains 
        xtensor<double,2> makeSdomArray (){

            xtensor<double,2> omg = zeros<double>({nSubDom,lSubDom});

            int lenArray = omg.shape()[1];

            for (int i = 0; i < nSubDom; i++){
                
                view(omg, i, range(0,lOlap)) = i+1; 
                view(omg, i, range(lenArray-lOlap,lenArray)) = i+1; 

                // Initializing first and last few elements as NaN 
                if (i == 0){
                    view(omg, i, range(0,lOlap)) = NAN;

                } else if (i == nSubDom-1){
                    view(omg, i, range(lenArray-lOlap,lenArray)) = NAN;

                } 
                
            }
            
            return omg; 
        }

    subdomArray(int rows, int cols, int l_olap){
        nSubDom = rows;
        lSubDom = cols; 
        lOlap = l_olap; 
    }
}; 

// Function for transferring data between subdomains 
xtensor<double,2>  transfer(xtensor<double,2> array_in, int lOlap){

    int lenArray = array_in.shape()[1];
    xtensor<double,1> temp = zeros<double>({lOlap}); 
     
     for (int i = 0; i < array_in.shape()[0] - 1; i++){
         temp = view(array_in, i, range(lenArray-lOlap, lenArray));
         view(array_in, i, range(lenArray-lOlap, lenArray)) = view(array_in, i+1, range(0,lOlap));
         view(array_in, i+1, range(0,lOlap)) = temp; 

     }

    return array_in; 
}

int main(){

    int nSubdoms = 4; 
    int lSubdoms = 10;
    int lOlap = 3; 

    subdomArray x = subdomArray(nSubdoms, lSubdoms, lOlap); 
    xtensor<double,2> omg = x.makeSdomArray();  

    cout << "Original matrix: \n" << omg << endl;

    xtensor<double,2> tv = transfer(omg, lOlap);
    cout << "Transfer of data: \n" << tv << endl;

    return 0; 
}

#endif 