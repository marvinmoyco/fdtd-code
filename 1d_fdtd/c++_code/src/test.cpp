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
#include "xtensor/xmanipulation.hpp"
using namespace std;
using namespace xt;
int main()
{
    xtensor<double,1> x1 = linspace<double>(1,10,10);
    xtensor<double,1> x2 = linspace<double>(11,20,10);
    xtensor<double,1> x3 = linspace<double>(21,30,10);
    xtensor<double,1> x4 = linspace<double>(31,40,10);
    auto padded_x1 = pad(x1,5,pad_mode::constant,0);
    auto padded_x2 = pad(x2,5,pad_mode::constant,1);
    auto padded_x3 = pad(x3,5,pad_mode::constant,2);
    auto padded_x4 = pad(x4,5,pad_mode::constant,3);
    xtensor<double,2> y;
 
    y = atleast_2d(padded_x1);
    y = vstack(xtuple(y,atleast_2d(padded_x2)));
    y = vstack(xtuple(y,atleast_2d(padded_x3)));
    y = vstack(xtuple(y,atleast_2d(padded_x4)));

    //cout << y << endl;


    xtensor<double,1> y_new;

    unsigned int overlap = 5;
    unsigned int non_overlap = 10;
    unsigned int start = overlap;
    unsigned int stop = non_overlap + overlap;

    for(int i =0; i< 4; i++)
    {
        if(i == 0)
        {
            y_new = view(y,i,range(start,_));
            
        }
        else{
            if(i == 3)
            {
                y_new = hstack(xtuple(y_new,view(y,i,range(start,stop))));
            }
            else{
                y_new = hstack(xtuple(y_new,view(y,i,range(start,_))));
            }
            
        }
        cout << y_new << endl;
    }

    cout << endl << endl << y_new << endl;
    
    return 0;
}