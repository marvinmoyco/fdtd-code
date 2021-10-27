#include <iostream>
#include <istream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <string>
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include <xtensor/xcsv.hpp>
#include <xtensor/xindex_view.hpp>
#include <chrono>
#include <thread>
#include <list>
#include <complex>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xio.hpp>

using namespace std::complex_literals;



using namespace std;
//using namespace xt;


int main(){
    //Shape: {row,column,sheet}
    xt::xtensor<int,1> x;
    x.resize({10});
    x = xt::linspace(1,10,10);
    xt::xtensor<int,1> y = xt::linspace(1,10,10);
    std::cout << x << endl;
    std::cout << y << endl;

    xt::xtensor<int,1> z = xt::zeros<int>({10});

    z = x * y;
    cout << z << endl;
    z = 55*(x / y);
    cout << z << endl;
    //row,col,sheet 1:5
    //range() -> [min,max)
    xt::view(z,xt::range(0,5)) = xt::view(z,xt::range(0,5))/5; //slicing
    cout << z << endl;
    z(5) = 99; //element indexing
    cout << z << endl;
    return 0;
}