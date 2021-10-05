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
#include <xtensor/xarray.hpp>
#include <xtensor/xcomplex.hpp>

using namespace std::complex_literals;



using namespace std;
using namespace xt;
xtensor<double,1> x;

int main(){


    /*xtensor<double,1> a {1.2,23.5,12.56,4100.89};
    xtensor<double,1> b {1,1,1,1};
    xtensor<double,1> c {4,4,4,4};
    xtensor<double,1> f;
    cout << f << endl;
    xtensor<double,1> d = arange(0.0,100.0,0.21254);
    xtensor<double,1> az = ceil(a);
    auto x = amax<double>(sqrt(a*c));
    double y = x(0);
    xshape<7> s;
    double j = 10;
    f.resize({j});
    print_options::set_precision(9);
    cout << a <<endl;
    cout << b << endl;
    cout << a+b << endl;
    cout << a-b << endl;
    cout << a*b << endl;
    cout << 5*b << endl;
    cout << b/a << endl; 
    cout << "f: " <<f << endl;
    cout << x(0)/y << endl;
    cout << y<< endl;
    cout << sum(az) << endl;
    cout << typeid(sum(az)(0)).name() << endl;
    cout <<  d << endl;
    //xtensor<double,1> b {-1,44,55,66,877,888};
    using namespace std::this_thread;
    using namespace std::chrono;
    xtensor<double,1> a {1,1,1,1,1,1,1,1,1,1};
    xtensor<double,1> b {1,1,1,1,1,1,1,1,1,1};
    xtensor<double,1> c {2,4,6,8,10};
    //cout << a << endl;
    //view(a,range(0,a.size()-1)) = view(a,range(0,a.size()-1)) + 5;
    //view(b,range(1,b.size())) = view(b,range(1,b.size())) + 5;


    xtensor<double,2> d {{1,2,3,4,5},{6,7,8,9,10}};
    
    /*for (int j = 1; j <=100; j++)
    {
        cout <<"\rIteration:"<<j<<"/"<<100 << flush;
        sleep_for(milliseconds(100));
    }
    cout << endl;
    //cout << a <<endl;
    cout << a.size() << endl;
    cout << a.size() -1 << endl;*/


/*
    xarray<std::complex<double>> e =
        {{1.0       , 1.0 + 1.0i},
        {1.0 - 1.0i, 1.0       }};

    real(e) = zeros<double>({2, 2});
    // => e = {{0.0, 0.0 + 1.0i}, {0.0 - 1.0i, 0.0}};
    cout << real(e) <<endl;

    cout << imag(e) << endl;
    cout << e << endl;*/
    double df = 1/(3.33564e-12*2100);
    double dt = 3.33564e-12;
    xtensor<double,1> g= {1,2,3,4,5,6};
     x = arange(0,(int) 2.4e9,df);
    cout << df << endl;
    cout << x << endl;
    double h = numeric_constants<double>::PI;
    xtensor<double,1> frange = linspace<double>(0,2.4e9,1000);
    cout << .5/3.33564e-12 << endl;
    cout << exp(-1i*2.0*h*dt*frange)<< endl;
    cout << x.size() << endl;
    cout << "g: " << g << endl;
    view(g,all()) = 0; 
    cout << "g: " << g << endl;
    xtensor<complex<double>,1> j = {1,2,3,4,5,6,7,8,9,10};
    xtensor<complex<double>,1> l = {1,2,3,4,5,6,7,8,9,10};
    cout << j/l << endl; 
    return 0;
}