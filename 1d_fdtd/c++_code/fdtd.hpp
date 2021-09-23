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
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include <xtensor/xcsv.hpp>
#include <xtensor/xindex_view.hpp>


/*
Temporary g++ command

g++ -std=c++14 main.cpp -I ../../../libs/xtensor/include -I ../../../libs/xtl/include -o main.out
*/

using namespace std;
using namespace xt;

//Global variables declaration
double c_0 = 299792458;
double  epsilon_0 = 8.8541878128E-12;
double  mu_0 = 1.25663706212E-6;


//Struct declarations...

typedef struct Simulation_parameters{
    double dz = 0; //Cell length
    double dt = 0;
    double Nz = 0;
    double Nt = 0;
    double fmax = 0;
    double sim_time = 0;
    int spacers = 0;
    int injection_point = 0;
    string source_type = "";

} simulation_parameters;

typedef struct Comp_domain{
    int injection_point = 0;
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
    xtensor<double,1> Reflectance;
    xtensor<double,1> Transmittance;
    
} fields;
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

        int GaussianSource(double t0_coeff = 6.0,double prop_coeff = 3.0,double tau_coeff = 12.0,double nmax = 1,double nsrc = 1)
        {
            //Calculate the necessary variables
            initialize(t0_coeff,prop_coeff,tau_coeff,nmax);
            cout << "========================================================================" << endl;
            cout << "t0: "<< source_param.t0 << " | tau: " << source_param.tau << " | T: " << source_param.sim_time << endl;

            //Calculate Nt 
            source_param.Nt = ceil(source_param.sim_time/source_param.dt);
            cout << "dt: " << source_param.dt << " seconds" << " | Nt: " << source_param.Nt << " iterations"<< endl;

            source_output.t = arange(0.0,source_param.Nt*source_param.dt,source_param.dt);

            //Computing the time input for electric field component
            source_output.t_E = (source_output.t - source_param.t0)/source_param.tau;
            

            //Computing the electric field component of the source
            source_output.Esrc = exp(-pow(source_output.t_E,2));
            
            //Computing the time input for the magnetic field component
            double adj_H = (nsrc*source_param.dz/(2*c_0)) + (source_param.dt/2);
            source_output.t_H = ((source_output.t - source_param.t0)+adj_H)/source_param.tau;

            //Computing the magnetic field component of the source
            source_output.Hsrc = -sqrt(1)*exp(-pow(source_output.t_H,2));

            cout << "Sizes of the electric and magnetic field component:" << endl;
            cout << "Esrc: " << source_output.Esrc.size() << " | Hsrc: " << source_output.Hsrc.size() << endl;
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


            //Selecting a sub array inside time-vector to apply exp()
            auto t_condition = filter(source_output.t,source_output.t < source_param.t0);
            long unsigned int t_condition_size = t_condition.size();


            //Computing the time input for electric field component
            source_output.t_E = (source_output.t - source_param.t0)/source_param.tau;
            auto t_E_initial = view(source_output.t_E,range(0,t_condition_size));
            //Computing the time input for the magnetic field component
            double adj_H = (nsrc*source_param.dz/(2*c_0)) + (source_param.dt/2);
            source_output.t_H = ((source_output.t - source_param.t0)+adj_H)/source_param.tau;
            auto t_H_initial = view(source_output.t_H,range(0,t_condition_size));

            //Resize the source field components 
            source_output.Esrc.resize(source_output.t.shape());
            source_output.Hsrc.resize(source_output.t.shape());

            //Computing the electric and magnetic field component of the source before t0
            view(source_output.Esrc,range(0,t_condition_size)) = exp(-pow(t_E_initial,2))*(sin(2*numeric_constants<double>::PI*source_param.fmax*t_condition));
            view(source_output.Hsrc,range(0,t_condition_size)) = exp(-pow(t_H_initial,2))*(sin(2*numeric_constants<double>::PI*source_param.fmax*t_condition));

            //Computing the electric field and magnetic field component of the source after t0
            view(source_output.Esrc,range(t_condition_size,source_param.Nt)) = (sin(2*numeric_constants<double>::PI*source_param.fmax*view(source_output.t,range(t_condition_size,source_param.Nt))));
            view(source_output.Esrc,range(t_condition_size,source_param.Nt)) = (sin(2*numeric_constants<double>::PI*source_param.fmax*view(source_output.t,range(t_condition_size,source_param.Nt))));

            cout << "Sizes of the electric and magnetic field component:" << endl;
            cout << "Esrc: " << source_output.Esrc.size() << " | Hsrc: " << source_output.Hsrc.size() << endl;
            return 0;
        }

        source_output_d get_computed_source()
        {
            return source_output;
        }




    private:

        int initialize(double t0_coeff = 1.0,double prop_coeff = 1.0,double tau_coeff = 1.0,double nmax = 1)
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
            else{
                cout << "ERROR: Incorrect source type!" << endl;
                return -1;
            }

            source_param.t0 = t0_coeff*source_param.tau;
            source_param.t_prop = (nmax*source_param.Nz*source_param.dz)/c_0;
            double initial_total_time = tau_coeff*source_param.tau +(prop_coeff*source_param.t_prop);
            cout << "initial_total: " << initial_total_time << " vs. simparam_sim_time: " << source_param.sim_time << endl;
            //Get the smaller total sim time to save memory
            if (source_param.sim_time > initial_total_time)
            {
                source_param.sim_time = initial_total_time;
            }


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
                for (int i =0; i<n_models;i++)
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
        int create_comp_domain(int spacer_cells = 0,int injection_point = 0)
        {
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
            }

            //Computing dz and Nz...
            double n_max =  amax<double>(sqrt(input.magnetic_permeability*input.electric_permittivity))(0);
            
            //Computing cell size based on the smallest wavelength (lambda_min)
            double lambda_min = c_0/(n_max*input.simulation_parameters.at(0));
            double delta_lambda = lambda_min/25; //denominator is  dependent on how many samples you want for a whole wave
            
            //Computing cell size based on the smallest layer size (min. dimension)
            double d_min = amin(input.layer_size)(0);
            double delta_size = d_min/25;  //denominator is the amount of cells that can resolve the smallest dimension
            
            //The final cell size is obtained by getting the smallest of delta_lambda and delta_size 
            //to make sure that the comp domain can resolve all the necessary features (wavelength or dimension)
            
            sim_param.dz = min(delta_lambda,delta_size)/2; //Dividing by 2 further decreases the cell size to make the resolution finer
            
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
                return -1;
            }
            else if(injection_point > sim_param.spacers)
            {
                cout << "Error detected: Injection point is inside the device model" <<endl;
                return -1;
            }
            else{
                if(injection_point == 0)
                {
                    sim_param.injection_point = (int) ceil(sim_param.spacers/5);
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
            sim_param.dt = (comp_domain.n(0)*sim_param.dz)/(2*c_0);
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
            else{
                cout << "ERROR: Invalid source type!";
                return -1;
            }

            //Compute the source that will be used in the project.
            compute_source();

            return 0;
        }
        

        int compute_source()
        {
            //Create a new Source object
            sim_source = new Source(sim_param);
            if(sim_param.source_type == "gaussian")
            {
                cout << "In Gaussian" << endl;
                sim_source->GaussianSource(6,3,12,amax(comp_domain.n)(0),comp_domain.n[sim_param.injection_point]);
            }
            else if(sim_param.source_type == "sinusoidal")
            {
                cout << "In Sinusoidal" << endl;
                sim_source->SinusoidalSource(3,3,3,amax(comp_domain.n)(0),comp_domain.n[sim_param.injection_point]);
            }

            cout << "Getting the output. FInished computing!" << endl;

            sim_source_fields = sim_source->get_computed_source();

            cout << "Esrc size: " << sim_source_fields.Esrc.size() << " | Hsrc size: " << sim_source_fields.Hsrc.size() << endl;
            return 0;
        }

};
#endif
