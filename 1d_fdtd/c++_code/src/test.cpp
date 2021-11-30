#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>     // xt::arange
#include <xtensor/xmath.hpp>        // xt::sin, cos
#include <complex>
#include <xtensor/xio.hpp>
#include <array>
#include <vector>
#include <xtensor/xtensor.hpp>
#include "xtensor/xview.hpp"
#include "xtensor/xpad.hpp"
using namespace std;
using namespace xt;
int main()
{
    xtensor<double,1> x = linspace<double>(1,10,10);
    cout << x << endl; 
    unsigned int a1 = 5;
    vector<size_t> b = {a1};
    xtensor<double,1> pad = zeros<double>(b);
    cout << pad << endl;
    cout << endl << endl;
    //Padding using the pad function
    xtensor<double,1> y = concatenate(xtuple(x,pad),0);
    cout << endl;
    cout << y << endl;
    xtensor<double,2> new_y = y.reshape({3,-1});
    cout << new_y.shape()[0]<< "," << new_y.shape()[1] << endl;
    cout << new_y << endl;
    cout << view(new_y,2,range(3,4)) << endl << transpose(view(new_y,range(0,2),4)) << endl<<endl;
    xtensor<double,1> a = concatenate(
                                       xtuple(view(new_y,2,range(3,4)),view(new_y,range(0,2),4)),
                                       0);
    cout << a << endl;
    xtensor<double,2> new_y1 =  hstack(xtuple(transpose(atleast_2d(a)),new_y));

    int x3 = 6;
    if(x3 > 1 &&
        x3 > 2 &&
        x3 > 5 &&
        x3 > 4)
    {
        cout << x3 <<endl;
    }
    
    cout << new_y1 
         << endl;
    return 0;
}