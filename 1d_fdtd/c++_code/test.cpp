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
using namespace xt;


int main(){

    xtensor<double,3> x = zeros<double>({3,3,3});
    cout << x << endl << endl;
    view(x,0,all(),all()) = 3;
    view(x,1,all(),all()) = 5;

    cout << x << endl;
    string filename = "out.npy";
    dump_npy(filename,x);
    return 0;
}