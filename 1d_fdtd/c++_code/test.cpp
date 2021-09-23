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

    ifstream input_stream;
    input_stream.open("input.csv");
    auto data = load_csv<double>(input_stream,',',1);

    cout << data << endl;
    for(auto i: col(data,3))
    {
        cout << i << endl;
    }
    return 0;
}
