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

using namespace std;
using namespace xt;

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
    cout <<  d << endl;*/
    //xtensor<double,1> b {-1,44,55,66,877,888};
    xtensor<double,1> a {1,2,3,4,5,6,7,8,9,10};
    auto b = filter(a,a<=5);
    //auto c = view(a,b);
    cout << a << endl << endl;
    cout << b << endl;
    return 0;
}
