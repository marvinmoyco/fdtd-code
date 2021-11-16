#ifndef FDTD
#define FDTD
 
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
#include <xtensor-fftw/basic.hpp>   // rfft, irfft
#include <xtensor-fftw/helper.hpp>  // rfftscale
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
double pi = numeric_constants<double>::PI;

//Struct declarations...

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

} save_data;


//Class declarations

class Source{
    public:
        source_parameters source_param;
        source_output_d source_output;

        Source(simulation_parameters sim_param)
        {
            source_param.fmax = sim_param.fmax;
            source_param.dz = sim_param.dz;
            source_param.dt = sim_param.dt;
            source_param.Nz = sim_param.Nz;
            source_param.sim_time = sim_param.sim_time;
            source_param.source_type = sim_param.source_type;
        }

        int GaussianSource(double t0_coeff = 6.0,double prop_coeff = 6.0,double tau_coeff = 12.0,double nmax = 1,double nsrc = 1)
        {
            //Calculate the necessary variables
            initialize(t0_coeff,prop_coeff,tau_coeff,nmax);
            cout << "========================================================================" << endl;
            cout << "t0: "<< source_param.t0 << " | tau: " << source_param.tau << " | T: " << source_param.sim_time << endl;

            //Calculate Nt 
            source_param.Nt = ceil(source_param.sim_time/source_param.dt);
            cout << "dt: " << source_param.dt << " seconds" << " | Nt: " << source_param.Nt << " iterations"<< endl;

            //source_output.t = arange<double>(0.0,(source_param.Nt*source_param.dt) ,source_param.dt);
            
            source_output.t = linspace<double>(0,source_param.Nt*source_param.dt,source_param.Nt);

            //Computing the time input for electric field component
            source_output.t_E = (source_output.t - source_param.t0)/source_param.tau;
            
            //Computing the electric field component of the source
            source_output.Esrc = exp(-pow(source_output.t_E,2));
            

            //Computing the time input for the magnetic field component
            double adj_H = (nsrc*source_param.dz/(2*c_0)) + (source_param.dt/2);
            source_output.t_H = ((source_output.t - source_param.t0)+adj_H)/source_param.tau;
            
            //Computing the magnetic field component of the source
            source_output.Hsrc = -sqrt(1)*exp(-pow(source_output.t_H,2));
            
            //cout << "Sizes of the electric and magnetic field component:" << endl;
            //cout << "Esrc: " << source_output.Esrc.size() << " | Hsrc: " << source_output.Hsrc.size() << endl;
            return 0;
        }
        
        int SinusoidalSource(double t0_coeff = 3,double prop_coeff = 3.0,double tau_coeff = 3,double nmax = 1,double nsrc = 1)
        {
            //Calculate the necessary variables
            initialize(t0_coeff,prop_coeff,tau_coeff,nmax);
            cout << "========================================================================" << endl;
            cout << "t0: "<< source_param.t0 << " | tau: " << source_param.tau << " | T: " << source_param.sim_time << endl;

            //Calculate Nt 
            source_param.Nt = ceil(source_param.sim_time/source_param.dt);
            cout << "dt: " << source_param.dt << " seconds" << " | Nt: " << source_param.Nt << " iterations"<< endl;

            
            source_output.t = arange(0.0,source_param.Nt*source_param.dt,source_param.dt);
            
            //source_output.t = linspace<double>(0,source_param.Nt*source_param.dt,source_param.Nt);
            
            //Selecting a sub array inside time-vector to apply exp()
            auto t_condition = filter(source_output.t,source_output.t < source_param.t0);
            long unsigned int t_condition_size = t_condition.size()-1;

            
            //Computing the time input for electric field component
            source_output.t_E = (source_output.t - source_param.t0)/source_param.tau;
            auto t_E_initial = view(source_output.t_E,range(0,t_condition_size+1));
            //Computing the time input for the magnetic field component
            double adj_H = (nsrc*source_param.dz/(2*c_0)) + (source_param.dt/2);
            source_output.t_H = ((source_output.t - source_param.t0)+adj_H)/source_param.tau;
            auto t_H_initial = view(source_output.t_H,range(0,t_condition_size+1));

            cout << "t_E: " << source_output.t_E << endl;
            cout << "t_H: " << source_output.t_H << endl;

            //Resize the source field components 
            source_output.Esrc.resize(source_output.t.shape());
            source_output.Hsrc.resize(source_output.t.shape());
            
            //Computing the electric and magnetic field component of the source before t0
            view(source_output.Esrc,range(0,t_condition_size+1)) = exp(-pow(t_E_initial,2))*(sin(2*numeric_constants<double>::PI*source_param.fmax*t_condition));
            view(source_output.Hsrc,range(0,t_condition_size+1)) = exp(-pow(t_H_initial,2))*(sin(2*numeric_constants<double>::PI*source_param.fmax*t_condition));
            //Computing the electric field and magnetic field component of the source after t0
            view(source_output.Esrc,range(t_condition_size,source_param.Nt)) = (sin(2*numeric_constants<double>::PI*source_param.fmax*view(source_output.t,range(t_condition_size,source_param.Nt))));
            view(source_output.Hsrc,range(t_condition_size,source_param.Nt)) = (sin(2*numeric_constants<double>::PI*source_param.fmax*view(source_output.t,range(t_condition_size,source_param.Nt))));

            return 0;
        }


        int SquareWaveSource(double delay = 0, double width = 0, double nmax = 1, double nsrc = 1)
        {
            initialize(delay,width,12,nmax);
            //Calculate the necessary variables
            cout << "========================================================================" << endl;
            cout << "t0: "<< source_param.t0 << " | tau: " << source_param.tau << " | T: " << source_param.sim_time << endl;

            //Calculate Nt 
            source_param.Nt = ceil(source_param.sim_time/source_param.dt);
            cout << "dt: " << source_param.dt << " seconds" << " | Nt: " << source_param.Nt << " iterations"<< endl;

            
            source_output.t = arange(0.0,source_param.Nt*source_param.dt,source_param.dt);
            

            //Resize the source field components 
            //width = tau
            //delay = t_0
            source_output.Esrc.resize(source_output.t.shape());
            source_output.Hsrc.resize(source_output.t.shape());
            int start_index = floor(source_param.tau/source_param.dt);
            int end_index = ceil((source_param.tau + source_param.t0)/source_param.dt);
            cout << "Start index: " << start_index << " | End index: " << end_index;

            //Computing the electric and magnetic field component of the source before t0
            view(source_output.Esrc,range(start_index,end_index)) = 1;
            source_output.Esrc(start_index-1) = 0.5;
            source_output.Esrc(end_index) = 0.5;
            //Computing the electric field and magnetic field component of the source after t0
            view(source_output.Hsrc,range(start_index,end_index)) = -1;
            source_output.Hsrc(start_index-1) = -0.5;
            source_output.Hsrc(end_index) = -0.5;
            return 0;
        }
        
        int ModulatedSineSource(double t0_coeff = 3,double prop_coeff = 3.0,double tau_coeff = 3,double nmax = 1,double nsrc = 1)
        {
            //Calculate the necessary variables
            initialize(t0_coeff,prop_coeff,tau_coeff,nmax);
            cout << "========================================================================" << endl;
            cout << "t0: "<< source_param.t0 << " | tau: " << source_param.tau << " | T: " << source_param.sim_time << endl;

            //Calculate Nt 
            source_param.Nt = ceil(source_param.sim_time/source_param.dt);
            cout << "dt: " << source_param.dt << " seconds" << " | Nt: " << source_param.Nt << " iterations"<< endl;

            
            source_output.t = arange(0.0,source_param.Nt*source_param.dt,source_param.dt);
            
            //source_output.t = linspace<double>(0,source_param.Nt*source_param.dt,source_param.Nt);
            //Computing the time input for electric field component
            source_output.t_E = (source_outvput.t - source_param.t0)/source_param.tau;
            
            //Computing the electric field component of the source
            source_output.Esrc = (sin(2*numeric_constants<double>::PI*source_param.fmax*source_output.t))*(exp(-pow(source_output.t_E,2)));
            

            //Computing the time input for the magnetic field component
            double adj_H = (nsrc*source_param.dz/(2*c_0)) + (source_param.dt/2);
            source_output.t_H = ((source_output.t - source_param.t0)+adj_H)/source_param.tau;
            
            //Computing the magnetic field component of the source
            source_output.Hsrc = -sqrt(1)*(sin(2*numeric_constants<double>::PI*source_param.fmax*source_output.t))*(exp(-pow(source_output.t_H,2)));
            
           
            return 0;
        }
        
        
        source_output_d get_computed_source()
        {
            return source_output;
        }




    private:

        int initialize(double t0_coeff = 3,double prop_coeff = 1.0,double tau_coeff = 1.0,double nmax = 1)
        {
            if(source_param.source_type == "gaussian")
            {
                //Computing necessary values for the gaussian pulse source
                source_param.tau = 0.5/source_param.fmax;

            }
            else if(source_param.source_type == "sinusoidal")
            {
                source_param.tau = 3/source_param.fmax;
            }
            else if(source_param.source_type == "square")
            {
                source_param.tau = 1/source_param.fmax;
            }
            else{
                cout << "ERROR: Incorrect source type!" << endl;
                return -1;
            }

            source_param.t0 = t0_coeff*source_param.tau;
            source_param.t_prop = (nmax*source_param.Nz*source_param.dz)/c_0;
            double initial_total_time = tau_coeff*source_param.tau +(prop_coeff*source_param.t_prop);
            //cout << "initial_total: " << initial_total_time << " vs. simparam_sim_time: " << source_param.sim_time << endl;
            //Get the smaller total sim time to save memory
            source_param.sim_time = initial_total_time;



            return 0;
        }







};

class Simulation
{

    public:
        //Initialization of the struct data types
        simulation_parameters sim_param;
        computational_domain comp_domain;
        input_data input;
        fields sim_fields;
        int spacer_cells;
        Source* sim_source;
        source_output_d sim_source_fields;
        save_data output;
        //Initializing variables in the constructor
        Simulation(string input_file="")
        {

            if(input_file.empty())  // If there is no input file path, use terminal or stdin/stdout for input/
            {
                /*
                Guide to indices of simulation_parameter vector
                -Index 0: fmax
                -Index 1: source type
                -Index 2: n_model
                -Index 3: t_sim

                */
                double fmax = 0;
                double source_type = 0;
                long unsigned int n_model = 0;
                double t_sim =0;
                cout << "========================================================================" << endl;
                cout << "Initializing simulation object" << endl;
                cout << "Maximum frequency (in Hz): ";
                cin >> fmax;
                cout << "Source Type (0 - Gaussian, 1 - Sinusoidal): ";
                cin >> source_type;
                cout << "Number of layers: ";
                cin >> n_model;
                cout << "Simulation Time (in seconds): ";
                cin >> t_sim;
                input.simulation_parameters.push_back(fmax);
                input.simulation_parameters.push_back(source_type);
                input.simulation_parameters.push_back(n_model);
                input.simulation_parameters.push_back(t_sim);

                //Resize RowVectors to reserve space
                input.layer_size.resize({n_model});
                input.magnetic_permeability.resize({n_model});
                input.electric_permittivity.resize({n_model});

                for(int i =0; i<input.simulation_parameters.at(2);i++)
                {
                    double l_size {0};
                    double l_mu {0};
                    double l_epsilon {0};
                    cout << "========================================================================" << endl;
                    cout << "Layer " << i+1 << " size (in meters): ";
                    cin >> l_size;
                    cout << "Layer " << i+1 << " magnetic permeability 'mu' (relative mu): ";
                    cin >> l_mu;
                    cout << "Layer " << i+1 << " electric permittivity 'epsilon' (relative epsilon): ";
                    cin >> l_epsilon;
                    input.layer_size(i) = l_size;
                    input.magnetic_permeability(i) = l_mu;
                    input.electric_permittivity(i) = l_epsilon;
                }
                cout << "========================================================================" << endl;
                cout << "Input Layer sizes: ";
                for(auto i: input.layer_size){
                    cout << i << " ";
                }
                cout << endl;
                cout << "Input Layer relative magnetic permeability: ";
                for(auto i: input.magnetic_permeability){
                    cout << i << " ";
                }
                cout << endl;
                cout << "Input Layer electric permittivity: ";
                for(auto i: input.electric_permittivity){
                    cout << i << " ";
                }
                cout << endl;
                

            }
            else{ //If there is a filename, process it by load_csv()
                ifstream input_stream;
                input_stream.open(input_file);
                auto input_data = load_csv<double>(input_stream,',',1);
                
                //Store the data into the intialized vectors
                //For simulation parameters
                auto temp_sim_param = col(input_data,3);
                cout << "========================================================================" << endl;
                cout << "INPUT FILE PROPERTIES" << endl;
                cout << "Simulation parameters:" << endl << "fmax,\t source_type,\t n_model,\t t_sim" << endl;
                for(int i=0; i < 4; i++)
                {
                    input.simulation_parameters.push_back(temp_sim_param(i));
                    cout << input.simulation_parameters.at(i) << "\t \t";
                }
                cout << endl;

                unsigned long int n_models = input.simulation_parameters.at(2);

                //Getting the columns for layer size, mu, and epsilon.
                auto temp_l_size = col(input_data,0);
                auto temp_mu = col(input_data,1);
                auto temp_epsilon = col(input_data,2);

                //Resizing the different xtensors
                input.layer_size.resize({n_models});
                input.magnetic_permeability.resize({n_models});
                input.electric_permittivity.resize({n_models});

                //Saving each layer into the vectors
                //For loop is necessary to truncate initial xarray and cutoff excess rows
                cout << endl <<"Device Model Configuration" << endl;
                cout << "Layer # \t Layer size \t Layer mu \t Layer epsilon" << endl;
                for (int i =0; i<(int) n_models;i++)
                {
                    cout << "Layer " << i+1 << ": \t";
                    input.layer_size(i) = temp_l_size(i);
                    input.magnetic_permeability(i)= temp_mu(i);
                    input.electric_permittivity(i) = temp_epsilon(i);
                    cout << input.layer_size(i) << " \t \t" << input.magnetic_permeability(i) << " \t \t" << input.electric_permittivity(i) << endl;
                    
                }

            }   

            
            
        }
        

        //Creating computational domain
        computational_domain create_comp_domain(int spacer_cells = 0,int injection_point = 0, double n_freq = 0)
        {
            /*
                This function is the "pre-processing" stage of the simulation. All of the needed pre-requisite computations are done before the actual
                simulation is done. The simulate() can be used after this function successfully finished.
            */

            
            //Try catch here to make sure that the input struct is not empty
            cout << "========================================================================" << endl;
            try
            {
                bool check_input = input.layer_size.size() != 0 && input.magnetic_permeability.size() != 0  && input.electric_permittivity.size() != 0  && !input.simulation_parameters.empty();
                if(check_input)
                {
                    cout << "Input are properly entered! Continuing the program...." << endl;
                }
                else
                {
                    throw -1;
                }
            }
            catch(...)
            {
                cout << "Error: Input are not yet entered." << endl;
                comp_domain.check = -1;

                exit(EXIT_FAILURE);
            }

            //Computing dz and Nz...
            double n_max =  amax<double>(sqrt(input.magnetic_permeability*input.electric_permittivity))(0);
            
            //Computing cell size based on the smallest wavelength (lambda_min)
            double lambda_min = c_0/(n_max*input.simulation_parameters.at(0));
            double delta_lambda = lambda_min/20; //denominator is  dependent on how many samples you want for a whole wave
            
            //Computing cell size based on the smallest layer size (min. dimension)
            double d_min = amin(input.layer_size)(0);
            double delta_size = d_min/10;  //denominator is the amount of cells that can resolve the smallest dimension
            
            //The final cell size is obtained by getting the smallest of delta_lambda and delta_size 
            //to make sure that the comp domain can resolve all the necessary features (wavelength or dimension)
            
            sim_param.dz = min(delta_lambda,delta_size); //Dividing by 2 further decreases the cell size to make the resolution finer
            
            //Get the total number cells needed for the device model
            xtensor<double,1> model_ncells = ceil((input.layer_size/sim_param.dz));
            sim_param.Nz = sum(model_ncells)(0);
            cout << "Cell amount per layer: " << model_ncells << endl;
            cout << "Number of cells of the device model: " << sim_param.Nz << " cells" << endl;
            cout << "Cell size (dz): "<< sim_param.dz << " m" << endl;
            //Check if there is specified spacer_region
            if(spacer_cells != 0)
            {
                sim_param.Nz += 2*(spacer_cells/sim_param.dz);//spacer_cells are in meters. Convert to cells
                sim_param.spacers = ceil(spacer_cells/sim_param.dz);
            }
            //If there is not, add half of the size of the simulation model to both end of comp domain.
            else
            {
                sim_param.spacers = (int) ceil((sim_param.Nz/2));
                sim_param.Nz += sim_param.Nz;
                
            }
            //Adjust Nz to be divisible by 2
            while(fmod(sim_param.Nz,2.0) != 0)
            {
                //To make Nz divisible by 2, add 1 cell to spacer regions until it becomes divisible by 2
                sim_param.Nz += 1;
            }

            //Check to make sure that the injection point is inside the spacer region. 
            if(injection_point < 0)
            {
                cout << "Error detected: Injection point is invalid" << endl;
                comp_domain.check = -1;
                exit(EXIT_FAILURE);
            }
            else if(injection_point > sim_param.spacers)
            {
                cout << "Error detected: Injection point is inside the device model" <<endl;
                comp_domain.check = -1;
                exit(EXIT_FAILURE);
            }
            else{
                if(injection_point == 0)
                {
                   
                    sim_param.injection_point = (int) ceil(sim_param.spacers/2);
                }
                else
                {
                    sim_param.injection_point = injection_point;
                }
            }
            cout << "Injection point (index/position): " << sim_param.injection_point << "-th cell" <<endl;
            cout << "Spacer cells: " << sim_param.spacers << " cells" <<endl;
            cout << "Total number of cells (model + spacers): " << sim_param.Nz << " cells" << endl;
            cout << "Total length (in m) of the computational domain: " << sim_param.Nz*sim_param.dz << " m" << endl;
            
            
            
            //Creating the computational domain
            comp_domain.z = arange(0.0,sim_param.Nz*sim_param.dz,sim_param.dz);
            //Creating the vectors for mu and epsilon
            //comp_domain.mu = ones<double>(comp_domain.z.shape());
            //comp_domain.epsilon = ones<double>(comp_domain.z.shape());

            comp_domain.mu.resize(comp_domain.z.shape());
            comp_domain.epsilon.resize(comp_domain.z.shape());

            view(comp_domain.mu,all()) = 1.0;
            view(comp_domain.epsilon,all()) = 1.0;

            //Assigning the mu and epsilon values to the proper index in the comp domain
            int start = sim_param.spacers;
            int end = sim_param.spacers;
            for (int i=0;i<input.simulation_parameters.at(2);i++)
            {
                end += model_ncells(i);
                //cout << "range[" << start << ":" << end << "]" << "mu value: " << input.magnetic_permeability(i) << endl;
                view(comp_domain.mu,range(start,end))= input.magnetic_permeability(i);
                view(comp_domain.epsilon,range(start,end)) = input.electric_permittivity(i);
                
                start = end;
            }

 
            /*cout << "+++++++++++++++=+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << endl;
            cout << comp_domain.mu << endl;
            cout << "+++++++++++++++=+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << endl;
            cout << comp_domain.epsilon << endl;*/

            // Computing for the refractive index vector based on the mu and epsilon vectors
            comp_domain.n.resize(comp_domain.z.shape());
            comp_domain.n = sqrt(comp_domain.mu*comp_domain.epsilon);

            //cout << comp_domain.n << endl;
            //Printing the layout of computational domain.
            cout << "========================================================================" << endl;
            cout << "Layout of Computational Domain: " << endl;
            cout << "|----spacer";
            for(int i =0; i < input.simulation_parameters.at(2);i++)
            {
                cout << "----model_layer_" << i+1;
            }
            cout << "----spacer----|" << endl << endl;
            cout << "|---" << sim_param.spacers << " cells";
            for(int i=0;i< input.simulation_parameters.at(2);i++)
            {
                cout << "---" << model_ncells(i) << " cells";
            }
            cout << "---" << sim_param.spacers << " cells---|" << endl;

            //Computing dt or CFL Condition
            sim_param.dt = (1*sim_param.dz)/(2*c_0);
            //Check the sim_time input in the csv file
            if (input.simulation_parameters.at(3) == 0)
            {
                sim_param.sim_time = (comp_domain.n(0)*sim_param.Nz*sim_param.dz)/c_0;
                
            }
            else{
                sim_param.sim_time = input.simulation_parameters.at(3);
            }

            cout << "========================================================================" << endl;
            cout << "Time step (dt): " << sim_param.dt << " seconds | Sim time: " << sim_param.sim_time << " seconds" << endl;

            //store fmax in the proper place
            sim_param.fmax = input.simulation_parameters.at(0);  
            if (input.simulation_parameters.at(1) == 0)
            {
                sim_param.source_type = "gaussian";
            }
            else if(input.simulation_parameters.at(1) == 1)
            {
                sim_param.source_type = "sinusoidal";

            }
            else if(input.simulation_parameters.at(1) == 2)
            {
                sim_param.source_type = "square";
            }
            else{
                cout << "ERROR: Invalid source type!";
                comp_domain.check = -1;
                exit(EXIT_FAILURE);
            }

            //Compute the source that will be used in the project.
            compute_source();
            
            //Resize the fields
            sim_fields.E.resize(comp_domain.z.shape());
            sim_fields.H.resize(comp_domain.z.shape());
            sim_fields.m_E.resize(comp_domain.z.shape());
            sim_fields.m_H.resize(comp_domain.z.shape());

            //Initialize fields to 0
            view(sim_fields.E,range(0,sim_fields.E.size())) = 0;
            view(sim_fields.H,range(0,sim_fields.E.size())) = 0;

            //Compute the update coefficients
            sim_fields.m_E = (c_0*sim_param.dt)/(comp_domain.epsilon*sim_param.dz);
            sim_fields.m_H = (c_0*sim_param.dt)/(comp_domain.mu*sim_param.dz);

            unsigned long int row = sim_param.Nt;
            unsigned long int col = sim_param.Nz;
            //Resize the different save matrices
            output.E.resize({row,col});
            output.H.resize({row,col});

            //Print shape
            //for (auto& el : output.E.shape()) {cout << el << ", ";}
            //cout << output.E.size() <<endl;

            //cout << comp_domain.z.size() << "   " << comp_domain.mu.size() << endl;

            //Computing df - frequency step in Freq. Response 
            
            sim_param.df = 1/(sim_param.dt*sim_param.Nt);
            

            //Initialize the frequency vector for Fourier Transform
            sim_param.n_freq = n_freq;
            sim_param.fmax_fft = 0.5*(1/sim_param.dt); //Upper frequency limit Fs/2
            sim_param.n_freq = sim_param.Nt; //Make sure to match the num of points in other FFT
            sim_fields.Freq_range = linspace<double>(0,sim_param.fmax_fft,sim_param.n_freq);
            cout << "Num. of freq. samples: " << sim_param.n_freq << " | Freq_range shape: (" << sim_fields.Freq_range.shape()[0] << "," << sim_fields.Freq_range.shape()[1] << endl;
            //Initialize the sizes of refl,trans, and kernal vectors
            sim_fields.Kernel_Freq = exp(-1i*2.0*pi*sim_param.dt*sim_fields.Freq_range);
            cout << "Kernel shape: (" << sim_fields.Kernel_Freq.shape()[0] << "," << sim_fields.Kernel_Freq.shape()[1] << ")" << endl; 

            //Initialize all related tensors for FFT
            sim_fields.Reflectance.resize(sim_fields.Kernel_Freq.shape());
            sim_fields.Transmittance.resize(sim_fields.Kernel_Freq.shape());
            sim_fields.Con_of_Energy.resize(sim_fields.Kernel_Freq.shape());
            sim_fields.Source_FFT.resize(sim_fields.Kernel_Freq.shape());

            //Initialize the values to 0
            view(sim_fields.Reflectance,all()) = 0;
            view(sim_fields.Transmittance,all()) = 0;
            view(sim_fields.Con_of_Energy,all()) = 0;
            view(sim_fields.Source_FFT,all()) = 0;


            //Initialize save matrices for FFT;
            unsigned long col_f = sim_param.n_freq;
            output.Reflectance.resize({row,col_f});
            output.Transmittance.resize({row,col_f});
            output.Con_of_Energy.resize({row,col_f});
            output.Source_FFT.resize({row,col_f});

            //Print the information computed
            cout << "========================================================================" << endl;
            cout << "df: " << sim_param.df << " | Frequency range: " << sim_fields.Freq_range.size() << " | Kernel: " << sim_fields.Kernel_Freq.size() << endl;
            cout << "fmax_fft: " << sim_param.fmax_fft << endl;
            cout << "Save Matrices Sizes:" << endl;
            cout << "Reflectance: " << output.Reflectance.size() << " | Transmittance: " << output.Transmittance.size() << " | Conservation of Energy: " << output.Con_of_Energy.size() << endl;
            return comp_domain;
        }
        

        int compute_source()
        {
            //Create a new Source object
            sim_source = new Source(sim_param);
            
            cout << "========================================================================" << endl;
            cout << "Computing source. | Selected source: " << sim_param.source_type << endl;
            if(sim_param.source_type == "gaussian")
            {
                sim_source->GaussianSource(3,12,2,amax(comp_domain.n)(0),comp_domain.n[sim_param.injection_point]);
            }
            else if(sim_param.source_type == "sinusoidal")
            {
                sim_source->SinusoidalSource(3,12,2,amax(comp_domain.n)(0),comp_domain.n[sim_param.injection_point]);
            }
            else if(sim_param.source_type == "square")
            {
                sim_source->SquareWaveSource(6,1,amax(comp_domain.n)(0),comp_domain.n[sim_param.injection_point]);
            }

            //Transfer Nt to sim_param
            sim_param.Nt = sim_source->source_param.Nt;
            sim_source_fields = sim_source->get_computed_source();

            unsigned long int row_s = 3;
            unsigned long int col_s = sim_param.Nt;
            //Resize Source xtensor
            output.Source.resize({row_s,col_s});
    
            //Add source to output_csv
            row(output.Source,0) = sim_source->source_output.t;
            row(output.Source,1) = sim_source->source_output.Esrc;
            row(output.Source,2) = sim_source->source_output.Hsrc;

            //Print the details
            cout << "Sizes of the electric and magnetic field component:" << endl;
            cout << "Esrc size: " << sim_source_fields.Esrc.size() << " | Hsrc size: " << sim_source_fields.Hsrc.size() << endl;
            return 0;
        }

        //FDTD ALgorithm only (serial version)
        save_data simulate(string boundary_condition = "", string excitation = "")
        {
            
            cout << "========================================================================" << endl;
            cout << "Starting 1-D FDTD Simulation..." << endl;
            cout << "Boundary Condition: " << boundary_condition << " | Source Excitation Method: " << excitation << endl;
            unsigned int end_index = sim_param.Nz;
            sim_param.boundary_cond = boundary_condition;
            sim_param.excitation_method = excitation;
            //Initialize variables used for outside boundary terms
            //double E_bounds = 0;
            //double H_bounds = 0;
            
            //Initialize buffer variables in FFT calculation
            xtensor<complex<double>,1> R = zeros<complex<double>>(sim_fields.Kernel_Freq.shape());
            xtensor<complex<double>,1> T = zeros<complex<double>>(sim_fields.Kernel_Freq.shape());
            

            //Remove in production
            comp_domain.injection_point = ceil(sim_param.Nz/2);
            cout << "Start of simulation." << endl;
            //cout << "m_E: " << sim_fields.m_E << endl;
            //cout << "m_H: " << sim_fields.m_H << endl;
            //FDTD Time Loop
            for(int curr_iteration = 0;curr_iteration < sim_param.Nt; curr_iteration++)
            {
                
                //H field computations
                //Store H boundary terms
                if (boundary_condition == "pabc")
                {
                    //cout << " H-pabc ";
                    //Get the front of the list
                    sim_fields.E(0) = sim_fields.H_start.front();
                    //Remove the front element from the list
                    sim_fields.H_start.pop_front();
                    //Add H[0] at the end of the list
                    sim_fields.H_start.push_back(sim_fields.E(1));
                    //Printing the contents of E_end:
                    /*cout << "z_low: [";
                    for (auto iter : sim_fields.H_start)
                    {
                        cout << iter << ",";
                    }
                    cout << "]" << endl;*/

                    //H_bounds = sim_fields.H(0);
                    //sim_fields.H(0) = sim_fields.H(1);
                }
                else if(boundary_condition == "dirichlet")
                {
                    //cout << " H-dirichlet ";
                    //H_bounds = 0;
                    sim_fields.E(0) = 0;
                }

   
                //Update H from E (FDTD Space Loop for H field)
                view(sim_fields.H,range(0,end_index-1)) = view(sim_fields.H,range(0,end_index-1)) + (view(sim_fields.m_H,range(0,end_index-1)))*(view(sim_fields.E,range(1,end_index)) - view(sim_fields.E,range(0,end_index-1)));

               


                //Inject the H source component (only needed when using TFSF)
                /*if(excitation == "hard")
                {
                    //cout << "H-hard ";
                    sim_fields.H(comp_domain.injection_point) = sim_source_fields.Hsrc(curr_iteration);
                }
                else if(excitation == "soft")
                {
                    //cout << "H-soft ";
                    sim_fields.H(comp_domain.injection_point) += sim_source_fields.Hsrc(curr_iteration);
                }*/
                if(excitation == "tfsf")
                {
                    //cout << "H-tfsf ";
                    sim_fields.H(comp_domain.injection_point-1) -= (sim_fields.m_H(comp_domain.injection_point-1)*sim_source_fields.Esrc(curr_iteration));
                }

                //cout << "E_bounds: " << E_bounds << endl;
                //Boundary conditions for H (at the end of the comp. domain)
                //sim_fields.H(end_index-1) = sim_fields.H(end_index-1) + (sim_fields.m_H(end_index-1) * (E_bounds - sim_fields.E(end_index-1)));




                //E field computations
                //Store E boundary terms
                if (boundary_condition == "pabc")
                {
                    //cout << "E-pabc ";
                    //Get the front of the list
                    sim_fields.H(end_index-1) = sim_fields.E_end.front();
                    //Remove the front element from the list
                    sim_fields.E_end.pop_front();
                    //Add E[Nz] at the end of the list
                    sim_fields.E_end.push_back(sim_fields.H(end_index-2)); //In other file H(end_index-2);
                    //Printing the contents of E_end:
                    /*cout << "z_high: [";
                    for (auto iter : sim_fields.E_end)
                    {
                        cout << iter << ",";
                    }
                    cout << "]" << endl;*/
                    //E_bounds = sim_fields.E(end_index -1);
                    //sim_fields.E(end_index-1) = sim_fields.E(end_index-2);
                }
                else if(boundary_condition == "dirichlet")
                {
                    //cout << "E-dirichlet ";
                    //E_bounds = 0;
                    sim_fields.H(end_index-1) = 0;
                }

                //cout << "H_bounds: " << H_bounds << endl;
                //Boundary condition for E (at the start of the comp.domain)
                //sim_fields.E(0) = sim_fields.E(0) + (sim_fields.m_E(0)*(sim_fields.H(0) - H_bounds));


                //Update E from H (FDTD Space Loop for E field)
                view(sim_fields.E,range(1,end_index)) =  view(sim_fields.E,range(1,end_index)) + (view(sim_fields.m_E,range(1,end_index))*(view(sim_fields.H,range(1,end_index))-view(sim_fields.H,range(0,end_index-1))));

                //Inject the E source component
                if(excitation == "hard")
                {
                    //cout << "E-hard ";
                    sim_fields.E(comp_domain.injection_point) = sim_source_fields.Esrc(curr_iteration);
                }
                else if(excitation == "soft")
                {
                    //cout << "E-soft";
                    sim_fields.E(comp_domain.injection_point) += sim_source_fields.Esrc(curr_iteration);
                }
                else if(excitation == "tfsf")
                {
                    //cout << "E-tfsf ";
                    sim_fields.E(comp_domain.injection_point) -= (sim_fields.m_E(comp_domain.injection_point)*sim_source_fields.Hsrc(curr_iteration));
                }
                

                //Compute for the Fourier Transform of the current simulation window
                //More cheaper than saving all the data, it will compute at each iteration
                R = R + (pow(sim_fields.Kernel_Freq,curr_iteration)*sim_fields.E(1));
                T = T + (pow(sim_fields.Kernel_Freq,curr_iteration)*sim_fields.E(end_index-1));
                sim_fields.Source_FFT = sim_fields.Source_FFT  + (pow(sim_fields.Kernel_Freq,curr_iteration)*sim_source_fields.Esrc(curr_iteration));

                //Normalize the computed Fourier Transform
                sim_fields.Reflectance = pow(real(R/sim_fields.Source_FFT ),2);
                sim_fields.Transmittance = pow(real(T/sim_fields.Source_FFT ),2);
                sim_fields.Con_of_Energy = sim_fields.Reflectance + sim_fields.Transmittance;

                //Save the computed fields into the save matrix
                //For the fields...
                row(output.E,curr_iteration) = sim_fields.E;
                row(output.H,curr_iteration) = sim_fields.H;
                row(output.Source_FFT,curr_iteration) = real(sim_fields.Source_FFT);
                //for the Fourier Transform
                row(output.Reflectance,curr_iteration) = real(sim_fields.Reflectance);
                row(output.Transmittance,curr_iteration) = real(sim_fields.Transmittance);
                row(output.Con_of_Energy,curr_iteration) = real(sim_fields.Con_of_Energy);

                //cout << endl;
                cout << "\rCurrent Iteration: "<<curr_iteration + 1<<"/"<<sim_param.Nt ;
            }
            //cout << endl << "Post-processing Fourier Transform" << endl;
            //sim_fields.Reflectance = pow(abs(sim_fields.Reflectance/sim_fields.Source_FFT),2);
            //sim_fields.Transmittance = pow(abs(sim_fields.Transmittance/sim_fields.Source_FFT),2);
            //sim_fields.Con_of_Energy = sim_fields.Reflectance + sim_fields.Transmittance;
            cout << endl << "End of simulation." << endl;

            cout << "Computing FFT using FFTW library..." << endl;

            //Compute FFT using FFTW lib
            //Resize the save data
            unsigned int size_fft = sim_param.Nt/2;
            output.FFTW_R.resize({size_fft,1});
            output.FFTW_T.resize({size_fft,1});
            output.FFTW_C.resize({size_fft,1});
            output.FFTW_S.resize({size_fft,1});
            output.Freq_Axis.resize({sim_param.Nt,1});
            cout << "Transfer xtensor to xarray...." << endl;
            //Transfer xtensor arrays to xarray
            xarray<double> source = zeros<double>({(int) sim_param.Nt,1});
            xarray<double> E_0 = zeros<double>({(int) sim_param.Nt,1});
            xarray<double> E_n = zeros<double>({(int) sim_param.Nt,1});
            view(source,all(),0) = sim_source_fields.Esrc;
            view(E_0,all(),0) = col(output.E,0);
            view(E_n,all(),0) = col(output.E,end_index-1);

            cout << "Compute FFTW..." << source.shape()[0] << ", " << source.shape()[1] << endl;
            cout << sim_source_fields.Esrc.shape()[0] << ", " << sim_source_fields.Esrc.shape()[1] << endl;
            cout << "Esrc (Xtensor): " << sim_source_fields.Esrc << endl;
            cout << "Esrc (Xarray): " << view(source,all(),0).shape()[0] << endl;


            //Compute for the Fourier transform
            output.FFTW_S = abs(real(fftw::rfft(source)));
            cout << "computing source fft.." << endl;
            output.FFTW_R = abs(real(fftw::rfft(E_0)/output.FFTW_S));
            output.FFTW_T = abs(real(fftw::rfft(E_n)/output.FFTW_S));
            output.FFTW_C = abs(real(output.FFTW_R + output.FFTW_T));
            output.FFTW_Freq = fftw::rfftfreq(sim_param.Nt,sim_param.dt);
            cout << "After FFTW_Freq " << output.Source_FFT.shape()[0] << ", " << output.Source_FFT.shape()[1] << endl;
            col(output.Freq_Axis,0) = linspace<double>(0,(sim_param.fmax_fft)*sim_param.df,(output.Source_FFT.shape()[0]));
            cout << "After Freq Axis" << endl;
            return output;
        }

        //FDTD Schwarz Serial Version
        void simulate_fdtd_schwarz()
        {

        }

        int write_to_csv(string output_file = "",xtensor<double,2> data ={{0,0,0},{0,0,0}})
        {
            ofstream out_stream;
            out_stream.open(output_file);
            dump_csv(out_stream,data);
            out_stream.close();
            return 0;
        }
        //Overloading the write_to_csv function to accept 2D complex data type
        int write_to_csv(string output_file = "", xtensor<complex<double>,2> data = {{0,0,0},{0,0,0}})
        {
            ofstream out_stream;
            out_stream.open(output_file);
            dump_csv(out_stream,data);
            out_stream.close();
            return 0;
        }

        int write_to_npy(string output_filename = "", xtensor<double,3> data = {{{0,0,0},{0,0,0}},{{0,0,0},{0,0,0}}})
        {
            
            //Store the time and source fields (column-wise; meaning each vector is stored in 1 column)
            //Resize data matrix
            long unsigned int row = sim_param.Nt;
            long unsigned int col = sim_param.Nz;
            data.resize({row,col,6});
            //Shape: (row,column,sheet) for 3D arrays
            view(data,0,0,0) = sim_param.n_freq;
            view(data,1,0,0) = sim_param.Nz;
            view(data,2,0,0) = sim_param.Nt;
            view(data,3,0,0) = sim_param.fmax_fft;
            view(data,4,0,0) = sim_param.injection_point;
            view(data,5,0,0) = sim_param.dt;
            view(data,6,0,0) = sim_param.dz;
            view(data,all(),1,0) = sim_source_fields.t;
            view(data,all(),2,0) = sim_source_fields.Esrc;
            view(data,all(),3,0) = sim_source_fields.Hsrc;
            //view(data,all(),4,0) = sim_fields.Freq_range;
            cout << view(data,0,all(),all()) << endl;
            //Store E field
            view(data,all(),all(),1) = output.E; 
            //Store H field
            view(data,all(),all(),2) = output.H; 
            //Store Refl
            view(data,all(),sim_param.n_freq,3) = cast<double>(output.Reflectance); 
            //Store Trans
            view(data,all(),sim_param.n_freq,4) = cast<double>(output.Transmittance); 
            //Store Conservation of Power
            cout << "Storing data from Conservation of Power" << endl;
            cout << view(data,all(),all(),5) << endl;
            view(data,all(),sim_param.n_freq,5) = cast<double>(output.Con_of_Energy);
            //write to a npy file
            dump_npy(output_filename,data);

            return 0;
        }

        auto write_to_hdf5(HighFive::File file, string dataset_path, auto data)
        {
            //format dump(filename, key,value)
            return dump(file, dataset_path,data,dump_mode::overwrite);
        }

        int save_to_file(string name = "",string type = "npy",string output_dir = "")
        {
            /*
            
                CSV - produces multiple files in csv format
                NPY - produces single file where the data is stored in 3D Matrix
                HDF5 - produces single file but can have key-value pairs
            */

            //get the current date
            auto now = chrono::system_clock::now();
            auto today = chrono::system_clock::to_time_t(now);
            stringstream string_stream;
            string_stream << put_time(localtime(&today),"%Y-%m-%d"); 
            string date_string = string_stream.str();

            if (type == "csv")
            {
                vector<string> names = {"source.csv","e_field.csv","h_field.csv","refl.csv","trans.csv","refl_trans.csv"};
                
                //cout << date_string + "_" + names[0] << endl;
                cout << "========================================================================" << endl;
                cout << "Saving data to csv files" << endl;
                cout << " Current Date: " << date_string << endl;
                
                



                //Save all data
                for(int i =0; i< (int)names.size();i++) //no reflectance and transmittance at the moment
                {
                    
                    //Concatenate the strings
                    string file_name = output_dir + date_string + "_";
                    if(name.empty())
                    {
                        file_name += names[i];
                    }
                    else{
                        file_name += name + "_" + names[i];
                    }
                    cout << "Filename: " << file_name ;

                    switch(i)
                    {
                        case 0: //for sources
                            //call write_to_csv
                            write_to_csv(file_name,output.Source);
                            break;

                        case 1: //for E-fields
                            //call write_to_csv
                            write_to_csv(file_name,output.E);
                            break;

                        case 2: //for H-fields
                            //call write_to_csv
                            write_to_csv(file_name,output.H);
                            break;

                        case 3: //for the reflectance
                            write_to_csv(file_name,output.Reflectance);
                            break;
                        
                        case 4: //for the transmittance 
                            write_to_csv(file_name,output.Transmittance);
                            break;

                        case 5:
                            write_to_csv(file_name,output.Con_of_Energy);
                            break;
                    }
                    cout << "-----> Successfully saved" << endl;
                    
                }
                
            }
            else if(type == "npy")
            {
                /*
                    3D matrix total size: (face,row,col) = (5,Nz,Nt)
                    npy file will save a 3D matrix containing everything from the simulation parameter to the simulated data.
                    1st 2D matrix (output[0,:,:]) will contain the simulation parameters (Nz, dz, etc.) and the time and source vectors.
                    2nd 2D matrix will contain E field simulation data
                    3rd 2D matrix will contain H field simulation data
                    4th 2D matrix will contain Reflectance simulation data
                    5th 2D matrix will contain Transmittance simulation data.
                */
                string npy_file_name = date_string + "_" + name + ".npy";
                xtensor<double,3> npy_data;
                cout << "Resizing 3D matrix" << endl;
                npy_data.resize({5,(unsigned int) sim_param.Nt,(unsigned int) sim_param.Nz});
                cout << "Initializing 3D Matrix" << endl;
                view(npy_data,all(),all(),all()) = 0;
                cout << "Writing to npy file" << endl;
                string output_filename = output_dir + npy_file_name;
                write_to_npy(output_filename ,npy_data);

            }
            else if(type == "hdf5")
            {
                //The format will be similar to a Python Dictionary, using a key-value pair
                string h5_file_name = output_dir + date_string + "_" + name + ".hdf5";
                cout << "========================================================================" << endl;
                cout << "Creating HDF5 file...." << endl;
                HighFive::File file(h5_file_name, HighFive::File::Overwrite);

                cout << "Saving simulation parameters" << "-----";
                //Store simulation parameters...
                write_to_hdf5(file, string("/Date of Simulation"), date_string);
                write_to_hdf5(file, string("/Cell size (dz)"), sim_param.dz);
                write_to_hdf5(file, string("/Total number of cells (Nz)"), sim_param.Nz);
                write_to_hdf5(file, string("/Time step (dt)"), sim_param.dt);
                write_to_hdf5(file, string("/Total number of time iteration (Nt)"), sim_param.Nt);
                write_to_hdf5(file, string("/Frequency of Interest (fmax)"), sim_param.fmax);
                write_to_hdf5(file, string("/Frequency Resolution (df)"), sim_param.df);
                write_to_hdf5(file, string("/Upper frequency limit (FFT)"), sim_param.fmax_fft);
                write_to_hdf5(file, string("/Total Simulation Time"), sim_param.sim_time);
                write_to_hdf5(file, string("/Number of frequencies (FFT)"), sim_param.n_freq);
                write_to_hdf5(file, string("/Amount of spacing (number of cells)"), sim_param.spacers);
                write_to_hdf5(file, string("/Source injection point (cell index)"), sim_param.injection_point);
                write_to_hdf5(file, string("/Source Type"), sim_param.source_type);
                write_to_hdf5(file, string("/Boundary Condition"), sim_param.boundary_cond);
                write_to_hdf5(file, string("/Source Excitation Method"), sim_param.excitation_method);
                cout << "--->" << "Saved" << endl;


                cout << "Saving Computational Domain parameters -----";
                //Store computational domain parameters
                write_to_hdf5(file, string("/Computational domain z (vector)"), comp_domain.z);
                write_to_hdf5(file, string("/Magnetic permeability mu (vector)"), comp_domain.mu);
                write_to_hdf5(file, string("/Electric permittivity epsilon (vector)"), comp_domain.epsilon);
                write_to_hdf5(file, string("/Refractive Index (vector)"), comp_domain.n);
                cout << "---> Saved!" << endl;

                cout << "Saving Input data -----";
                //Store input data (parsed from the input file)
                write_to_hdf5(file, string("/Input Layer size"), input.layer_size);
                write_to_hdf5(file, string("/Magnetic Permeability"), input.magnetic_permeability);
                write_to_hdf5(file, string("/Electric Permittivity"), input.electric_permittivity);
                cout << "---> Saved!" << endl;

                cout << "Saving Source parameters -----";
                //Store source parameters
                write_to_hdf5(file, string("/Source tau"), sim_source->source_param.tau);
                write_to_hdf5(file, string("/Source t0"), sim_source->source_param.t0);
                write_to_hdf5(file, string("/Source t_prop"), sim_source->source_param.t_prop);
                write_to_hdf5(file, string("/Esrc"), sim_source_fields.Esrc);
                write_to_hdf5(file, string("/Hsrc"), sim_source_fields.Hsrc);
                write_to_hdf5(file, string("/t_E"), sim_source_fields.t_E);
                write_to_hdf5(file, string("/t_H"), sim_source_fields.t_H);
                write_to_hdf5(file, string("/t"), sim_source_fields.t);
                cout << "---> Saved!" << endl;


                cout << "Saving Simulation Data -----";
                //Store the main simulation data
                write_to_hdf5(file, string("/m_E"), sim_fields.m_E);
                write_to_hdf5(file, string("/m_H"), sim_fields.m_H);
                write_to_hdf5(file, string("/FFT Frequency Range"), sim_fields.Freq_range);
                write_to_hdf5(file, string("/FFT Kernel Frequency Vector"), sim_fields.Kernel_Freq);
                write_to_hdf5(file, string("/E"), output.E);
                write_to_hdf5(file, string("/H"), output.H);
                write_to_hdf5(file, string("/Reflectance"), output.Reflectance);
                write_to_hdf5(file, string("/Transmittance"), output.Transmittance);
                write_to_hdf5(file, string("/Conservation of Energy"), output.Con_of_Energy);
                write_to_hdf5(file, string("/Source"), output.Source);
                write_to_hdf5(file, string("/Freq_Axis"), output.Freq_Axis);
                write_to_hdf5(file, string("/Source_FFT"), output.Source_FFT);
                
                //Saving the FFTW-computed data
                write_to_hdf5(file,string("/FFTW_R"),output.FFTW_R);
                write_to_hdf5(file,string("/FFTW_T"),output.FFTW_T);
                write_to_hdf5(file,string("/FFTW_C"),output.FFTW_C);
                write_to_hdf5(file,string("/FFTW_S"),output.FFTW_S);
                write_to_hdf5(file,string("/FFTW_Freq"),output.FFTW_Freq);
                cout << "---> Saved!" << endl;


                cout << "End....";
                cout << "========================================================================" << endl;
            }
            else
            {
                cout << "========================================================================" << endl;
                cout << "Invalid file type. Data will not be saved.";
            }
            
            
            return 0;
        }
};
#endif
