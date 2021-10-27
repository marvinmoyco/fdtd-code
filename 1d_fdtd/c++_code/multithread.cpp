#include <iostream>
#include <omp.h>
#include <mutex>  

//Xtensor preprocessor directives
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xcomplex.hpp"
#include "xtensor/xnpy.hpp"

#define N_THREADS 2

using namespace std;
using namespace xt;
mutex mtx;

int main() 
{
    //Set the number of threads to be created before creating "parallel regions"
    omp_set_num_threads(N_THREADS);
    xtensor<double,2> x = zeros<double>({10000,10000});
    xtensor<double,1> x_1 = zeros<double>({1000});
    cout << "x: " << x << endl;
    cout << "=================================================================" << endl;
    cout << "Starting parallel threads...." << endl;
    double y = 0;
 
    #pragma omp parallel
    {
        
        //mtx.lock();
        cout << "=================================================================" << endl;
        int thread_id = omp_get_thread_num();
        cout << "Number of threads: " << omp_get_num_threads() << endl;
        cout << "Current thread number: " << thread_id << endl;
        if (thread_id == 0) //Identify the main thread
        {
            #pragma omp critical
            {
                cout << "This is the main thread!" << endl;
                cout << "=================================================================" << endl;
                cout << " Thread " << thread_id << " is processing...." << endl;
                y += 5;
                cout << "Thread " << thread_id << " is finished..." << endl;
            }
            
            
        }
        else if(thread_id == 1) //Second thread
        {   
            #pragma omp critical
            {
                cout << "=================================================================" << endl;
                cout << " Thread " << thread_id << " is processing...." << endl;
                y -= 5;
                cout << "Thread " << thread_id << " is finished..." << endl;
            }
            
        }

        
        
        #pragma omp barrier
        cout << "Finished assigning new values" << endl;
        //mtx.unlock();

    }
    cout << "=================================================================" << endl;
    cout << "End of parallel threads...." << endl;

   // cout << "x: " << x << endl;

    //cout << "x_1: " << x_1 << endl;

    
    cout << "y = " << y << endl;
    return 0;
}