#ifndef COMMON
#define COMMON

#include <iostream>
#include <istream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include <ctime>
#include <thread>
#include <list>
#include <complex>

//Xtensor preprocessor directives
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xcomplex.hpp"
#include "xtensor/xnpy.hpp"
#include "xtensor/xio.hpp"
#include "xtensor-fftw/basic.hpp"   // rfft, irfft
#include "xtensor-fftw/helper.hpp"  // rfftscale
#include "xtensor-io/xhighfive.hpp"

/*
Temporary g++ command

g++ -std=c++14 main.cpp -I ../../../libs/xtensor/include -I ../../../libs/xtl/include -o main.out
*/

using namespace std;
using namespace xt;
using namespace std::this_thread;
using namespace std::chrono;
using namespace std::complex_literals;

//Global variables declaration 
double c_0 = 299792458;
double  epsilon_0 = 8.8541878128E-12;
double  mu_0 = 1.25663706212E-6;

//Struct declarations...
typedef struct Subdomain_data
{
    unsigned int subdomain_id;
    unsigned int subdomain_size;
    unsigned int overlap_size;
    unsigned int num_subdomains = 0;
    bool source_inject = false;
    bool preprocessed = false;
    double dz = 0; //Cell length
    double dt = 0;
    double Nz = 0;
    double Nt = 0;
    int spacers = 0;
    int injection_point = 0;

    string boundary_condition = "";
    string excitation_method = "";

} subdomain_data;

typedef struct Subdomain_fields
{
    xtensor<double,1> E;
    xtensor<double,1> H;
    xtensor<double,1> m_E;
    xtensor<double,1> m_H;
    list<double> E_end {0,0};
    list<double> H_start {0,0};
} subdomain_fields;

typedef struct Simulation_parameters{
    double dz = 0; //Cell length
    double dt = 0;
    double Nz = 0;
    double Nt = 0;
    double fmax = 0;
    double df = 0; // for the Fourier Transform'
    double fmax_fft = 0; //Max freq that can resolved by FDTD in DFT
    double sim_time = 0;
    double n_freq = 0; //For Fourier Transform, indicates how many points
    int spacers = 0;
    int injection_point = 0;
    
    unsigned int subdomain_size = 0;
    unsigned int num_subdomains = 0;
    unsigned int overlap_size = 0;

    string source_type = "";
    string boundary_cond = "";
    string excitation_method = "";

} simulation_parameters;

typedef struct Comp_domain{
    int injection_point = 0;
    int check = 0;
    xtensor<double,1> z; //Computational domain base vector. Basis for other vectors
    xtensor<double,1> mu; //Magnetic permeability vector
    xtensor<double,1> epsilon; // Electric permittivity vector
    xtensor<double,1> n; //Refractive index vector
} computational_domain;

typedef struct Input_data{

    xtensor<double,1> layer_size;
    xtensor<double,1>  magnetic_permeability;
    xtensor<double,1>  electric_permittivity;
    vector<double>  simulation_parameters;

} input_data;

typedef struct Source_parameters{

    double tau = 0;
    double t0 = 0;
    double Nt = 0;
    double Nz = 0;
    double fmax = 0;
    double dt = 0;
    double dz = 0;
    double sim_time = 0;
    double t_prop = 0;
    string source_type = "";
    

} source_parameters;

typedef struct Source_output{

    xtensor<double,1> t;
    xtensor<double,1> Esrc;
    xtensor<double,1> Hsrc;
    xtensor<double,1> t_E;
    xtensor<double,1> t_H;

} source_output_d;

typedef struct Simulation_Fields{
    xtensor<double,1> E;
    xtensor<double,1> H;
    xtensor<complex<double>,1> Reflectance;
    xtensor<complex<double>,1> Transmittance;
    xtensor<complex<double>,1> Con_of_Energy;
    xtensor<complex<double>,1> Kernel_Freq;
    xtensor<complex<double>,1> Source_FFT;
    xtensor<double,1> Freq_range;
    xtensor<double,1> m_E;
    xtensor<double,1> m_H;
    list<double> E_end {0,0};
    list<double> H_start {0,0};

    
} fields;

typedef struct Save_Data{
    xtensor<double,2> E;
    xtensor<double,2> H;
    xtensor<double,2> Reflectance;
    xtensor<double,2> Transmittance;
    xtensor<double,2> Con_of_Energy;
    xtensor<double,2> Source;
    xtensor<double,2> Freq_Axis;
    xtensor<double,2> Source_FFT;

    xtensor<double,2> FFTW_R;
    xtensor<double,2> FFTW_T;
    xtensor<double,2> FFTW_C;
    xtensor<double,2> FFTW_S;
    xtensor<double,2> FFTW_Freq;

    bool simulation_success = false;

} save_data;


#endif