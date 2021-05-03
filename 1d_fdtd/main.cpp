/*
    1D FDTD C++ Code
    
    A simple 1-dimensional FDTD simulator in Free space using Gaussian Pulse as source excitation.

    Based on the lecture videos of Dr. Rumpf (link: https://empossible.net/academics/emp5304/)
    and the book "ELECTROMAGNETIC SIMULATION USING THE FDTD METHOD WITH PYTHON"


*/


#include <iostream>
using namespace std;
#include <vector>
#include <cmath>
//#include "gnuplot-iostream.h"
#include "functions.hpp"
int main()
{
    //User input
    double f_max = 0;
    cout << "Enter the max frequency for the source excitation: ";
    cin >> f_max;
    cout << "The max frequency selected is " << f_max << "\n";

    //Computing grid resolution
    int n_max = 1; //Due to free space
    double lambda_min = c_0/(f_max*n_max);
    int N_lambda = 10;
    double delta_lambda = lambda_min/N_lambda;
    float d_min = 0.01;
    double delta_d = d_min/4;
    double delta_init = min(delta_lambda,delta_d);
    unsigned int spacers = 15;
    float d_critical = 0.7;
    int Nz = ceil(d_critical/delta_init) + spacers;
    double delta_z = d_critical/Nz;
    cout << "==============================================\n";
    cout << "Number of Yee cells: " << Nz << " cells\n";
    cout << "Length of the sides (each cell): " << delta_z << " m\n";

    //Computing material properties
    double mu_r[Nz] = {1};
    double epsilon_r[Nz] = {1};
    double n_r = sqrt(mu_r[0]/epsilon_r[0]);

    //Computing Time step and source excitation
    float n_bc = 1;
    double delta_t = (n_bc*delta_z)/(2*c_0);
    double t_prop = (n_r*Nz*delta_z)/c_0;
    


    return 0;
}