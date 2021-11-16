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

// g++ -std=c++17 schwarz_test.cpp -I ../../../xtensor/include -I ../../../xtl/include -o schwarz_test.out -o schwarz_test


/*xtensor<double, 1> createSubdoms(double lDom, double nSubDoms){

    double lSubDom = ceil(lDom/nSubDoms);
    double mod = fmod(lDom, nSubDoms);

    xtensor<double, 1> omg1 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg2 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg3 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg4; 

    if (mod > 0){
        double last = lSubDom - (nSubDoms - mod);
        omg4 = xt::zeros<double>({last});

    } else{
        omg4 = xt::zeros<double>({lSubDom});

    }

    return omg1, omg2, omg3, omg4;

}*/

int main(){

    double lDom = 20; 
    double nSubDoms = 5; 

    xtensor<double,1> compDom;
    // compDom = xt::zeros<double>({20});
    compDom = xt::random::rand<double>({lDom});


    double lSubDom = ceil(lDom/nSubDoms);
    cout << "lSubDom = " << lSubDom << endl;
    double mod = fmod(lDom, nSubDoms);
    cout << "mod = " << mod << endl;

    // Creating Subdomains 
    xtensor<double, 2> **subdomArray;
    //double **subdomArray; 

    subdomArray = new xtensor<double, 2> *[(int)nSubDoms]; 
    //subdomArray = new double *[(int)nSubDoms]; 

    for (int i = 0; i < nSubDoms; ++i){

        subdomArray[i] = new xtensor<double, 2>[(int)lSubDom];
        //subdomArray[i] = new double[(int)lSubDom]; 


        //subdomArray[i][0] = view(compdom, range()); 

        //view(subdomArray, 0, all() ) = view(compDom, range(0,4));

        auto k = view(&subdomArray[0], 0, all()).shape();
        auto l = view(compDom, range(0,4)).shape(); 

        cout << k[0] << " " <<  k[1] << endl; 
        cout << l[0] << " " << l[1] << endl;


        //cast<double>(**subdomArray);

        /*for (int j = 0; j < lDom; j++){
            if (j < lSubDom){
                subdomArray[i][j] = compDom(j);

            } else if (j > lSubDom && j < lSubDom*(i-1) ){
                subdomArray[i][j] = compDom(j);

            } else {
                subdomArray[i][j] = compDom(j);
            }
            
        }*/ 
        //cout << subdomArray[0][1] << endl;
        //cout << sizeof(subdomArray) << endl;

        

        //if(i == 0){
            //view(subdomArray[i], all()) = 1; 
            //view(subdomArray[i],xt::range(0,5)) = xt::view(subdomArray[i],xt::range(0,5)) + 1;
            //view(subdomArray,range(0,5)) = view(subdomArray,range(0,5)) +1;
        //} else if(i < nSubDoms-1){
            
        //} else {
            
        //}

        //cout << subdomArray[i] << endl;

        /* for (int j = 0; i < lDom; j++)
        {
            if (j < lSubDom)
            {
                subdomArray[i][j] ;
            }
            
        }*/ 
        
    }
    

	for( int i = 0; i < nSubDoms; ++i )
    	delete[] subdomArray[i];
 
	delete[] subdomArray;



    //cout << subdomArray[0][0] << endl; 


    //xt::view(subDom[i], xt::all()) = xt::view(compDom, xt::range(0, lSubDom))''
          
    //xt::view(subDom[i], xt::all()) = xt::view(compDom, xt::range((j-1)*lSubDom +1, j*lSubDom));

    //xt::view(subDom[i], xt::all()) = xt::view(compDom, xt::range((j-1)*lSubDom +1, lDom));

    // 4 subdomains     
    /*xtensor<double, 1> omg1 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg2 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg3 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg4; 

    if (mod > 0){
        double last = lSubDom - (nSubDoms - mod);
        omg4 = xt::zeros<double>({last});

    } else{
        omg4 = xt::zeros<double>({lSubDom});

    }*/ 


    /*cout << "omg1: " << omg1 << endl;
    cout << "omg2: " << omg2 << endl; 
    cout << "omg3: " << omg3 << endl; 
    cout << "omg4: " << omg4 << endl;*/ 

    // 8 subdomains     
    /*xtensor<double, 1> omg1 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg2 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg3 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg4 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg5 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg6 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg7 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg8; 

    if (mod > 0){
        double last = lSubDom - (nSubDoms - mod);
        cout << "last: " << last << endl;
        omg8 = xt::zeros<double>({last});

    } else{
        omg8 = xt::zeros<double>({lSubDom});

    }*/ 

    // 16 subdomains
    /*xtensor<double, 1> omg1 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg2 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg3 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg4 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg5 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg6 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg7 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg8 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg9 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg10 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg11 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg12 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg13 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg14 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg15 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg16; 

    if (mod > 0){
        double last = lSubDom - (nSubDoms - mod);
        cout << "last: " << last << endl;
        omg16 = xt::zeros<double>({last});

    } else{
        omg16 = xt::zeros<double>({lSubDom});

    }*/

    // 32 subdomains
    /*xtensor<double, 1> omg1 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg17 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg2 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg18 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg3 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg19 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg4 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg20 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg5 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg21 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg6 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg22 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg7 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg23 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg8 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg24 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg9 = xt::zeros<double>({lSubDom});     xtensor<double, 1> omg25 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg10 = xt::zeros<double>({lSubDom});    xtensor<double, 1> omg26 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg11 = xt::zeros<double>({lSubDom});    xtensor<double, 1> omg27 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg12 = xt::zeros<double>({lSubDom});    xtensor<double, 1> omg28 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg13 = xt::zeros<double>({lSubDom});    xtensor<double, 1> omg29 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg14 = xt::zeros<double>({lSubDom});    xtensor<double, 1> omg30 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg15 = xt::zeros<double>({lSubDom});    xtensor<double, 1> omg31 = xt::zeros<double>({lSubDom});
    xtensor<double, 1> omg16 = xt::zeros<double>({lSubDom});    xtensor<double, 1> omg32;

    if (mod > 0){
        double last = lSubDom - (nSubDoms - mod);
        cout << "last: " << last << endl;
        omg32 = xt::zeros<double>({last});

    } else{
        omg32 = xt::zeros<double>({lSubDom});

    }
  
    cout << "omg1: " << omg29.size() << endl;
    cout << "omg2: " << omg30.size() << endl; 
    cout << "omg3: " << omg31.size() << endl; 
    cout << "omg4: " << omg32.size() << endl;*/

    // 
    /*xtensor<double, 2> x;// = xt::xarray<double>({{3, 4, 1}, {5, 6}}); 
    x.resize({(long unsigned int)nSubDoms,(long unsigned int)lSubDom});
    cout << x << endl;*/ 

 
 


    /*for(size_t i = 0; i < nSubDoms; ++i){
        xtensor<double, 1> subDom[i];
        subDom[i].resize(lSubDom);

        for(size_t j = 0; j < lSubDom; j++){

            if(j == 0){
                xt::view(subDom[i], xt::all()) = xt::view(compDom, xt::range(0, lSubDom));

            } else if(j < lSubDom){
                xt::view(subDom[i], xt::all()) = xt::view(compDom, xt::range((j-1)*lSubDom +1, j*lSubDom));

            } else{
                xt::view(subDom[i], xt::all()) = xt::view(compDom, xt::range((j-1)*lSubDom +1, lDom));
            }
        }
        
    }*/ 

    /*xtensor<double, 2> x = xt::zeros<double>({3,3});
    xtensor<int, 1> y = {1,2,3};
    //view(x,range(0,5)) = view(x,range(0,5)) +1; 
    view(x, 0, all()) = view(y, all());

    view(x, 0, all()) = view(y, all());
    cout << x << endl; */


    
 





    return 0; 
}

#endif