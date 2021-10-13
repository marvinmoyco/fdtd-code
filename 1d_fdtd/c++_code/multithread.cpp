#include <iostream>
#include <omp.h>
#define N_THREADS 5

using namespace std;

int main() 
{
    //Set the number of threads to be created before creating "parallel regions"
    omp_set_num_threads(N_THREADS);

    cout << "=================================================================" << endl;
    cout << "Starting parallel threads...." << endl;
    int array[N_THREADS];
    //Initialize the arrays to -1
    cout << endl << endl << "Array: [";
    for(int j = 0; j < N_THREADS; j++)
    {
        array[j] = -1;
        cout << array[j] << ", ";
    }
    cout << "]" << endl;
    #pragma omp parallel
    {
        cout << "=================================================================" << endl;
        int thread_id = omp_get_thread_num();
        cout << "Number of threads: " << omp_get_num_threads() << endl;
        cout << "Current thread number: " << thread_id << endl;
        if (thread_id == 0) //Identify the main thread
        {
            cout << "This is the main thread!" << endl;
        }
        //Access the array variable and edit (1) element
        array[thread_id] = thread_id;
    }
    cout << "=================================================================" << endl;
    cout << "End of parallel threads...." << endl;

    cout << endl << endl << "Array: [";
    for(int j = 0; j < N_THREADS; j++)
    {
        cout << array[j] << ", ";
    }
    cout << "]" << endl;
    return 0;
}