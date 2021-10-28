#ifndef FDTD
#define FDTD
 
#include <iostream>
#include <Eigen/Dense>
#include <typeinfo>
#include <stdio.h>
#include <vector>
#include <cmath>


using namespace Eigen;
using namespace std;

//Global variables declaration
double c_0 = 299792458;
double epsilon_0 = 8.8541878128E-12;
double mu_0 = 1.25663706212E-6;


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
} simulation_parameters;

typedef struct Comp_domain{
    int injection_point = 0;
    RowVectorXd z; //Computational domain base vector. Basis for other vectors
    RowVectorXd mu; //Magnetic permeability vector
    RowVectorXd epsilon; // Electric permittivity vector
    RowVectorXd n; //Refractive index vector
} computational_domain;

typedef struct Input_data{

    RowVectorXd layer_size;
    RowVectorXd  magnetic_permeability;
    RowVectorXd  electric_permittivity;
    vector<double>  simulation_parameters;

} input_data;

typedef struct Source_parameters{

    double tau = 0;
    double t0 = 0;
    double Nt = 0;
    double fmax = 0;
    double dt = 0;
    double dz = 0;
    double sim_time = 0;

} source_param;

typedef struct Source_output{

    RowVectorXd t;
    RowVectorXd Esrc;
    RowVectorXd Hrc;

} source_output;

typedef struct Simulation_Fields{
    RowVectorXd E;
    RowVectorXd H;
    RowVectorXd Reflectance;
    RowVectorXd Transmittance;
    
} fields;
//Class declarations

//class Source{
    
//};

class Simulation
{

    public:
        //Initialization of the struct data types
        simulation_parameters sim_param;
        computational_domain comp_domain;
        input_data input;
        fields sim_fields;
        int spacer_cells;
        //Source sim_source;

        //Initializing variables in the constructor
        Simulation(string input_file="")
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
            double n_model = 0;
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
            input.layer_size.resize(n_model);
            input.magnetic_permeability.resize(n_model);
            input.electric_permittivity.resize(n_model);
            

            for(int i =0; i<n_model;i++)
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
            cout << "Product: " << input.magnetic_permeability.cwiseProduct(input.electric_permittivity) << endl;
            
            
            //Computing dz and Nz...
            double n_max = ((input.magnetic_permeability.cwiseProduct(input.electric_permittivity)).cwiseSqrt()).maxCoeff();
            cout << "n_max: "<< n_max << "\tfmax: " << input.simulation_parameters.at(0) <<endl;
            
            //Computing cell size based on the smallest wavelength (lambda_min)
            double lambda_min = c_0/(n_max*input.simulation_parameters.at(0));
            double delta_lambda = lambda_min/25; //denominator is  dependent on how many samples you want for a whole wave
            
            //Computing cell size based on the smallest layer size (min. dimension)
            double d_min = input.layer_size.minCoeff();
            double delta_size = d_min/25;  //denominator is the amount of cells that can resolve the smallest dimension
            
            //The final cell size is obtained by getting the smallest of delta_lambda and delta_size 
            //to make sure that the comp domain can resolve all the necessary features (wavelength or dimension)
            sim_param.dz = min(delta_lambda,delta_size)/2; //Dividing by 2 further decreases the cell size to make the resolution finer
            
            //Get the total number cells needed for the device model
            RowVectorXd model_ncells = (input.layer_size/sim_param.dz).array().ceil();
            sim_param.Nz = model_ncells.sum();
            cout << "Model cells: " << sim_param.Nz << endl;
            cout << "dz: "<< sim_param.dz << endl;
            //Check if there is specified spacer_region
            if(spacer_cells != 0)
            {
                sim_param.Nz += 2*(spacer_cells/sim_param.dz);//spacer_cells are in meters. Convert to cells
                sim_param.spacers = ceil(spacer_cells/sim_param.dz);
            }
            //If there is not, add half of the size of the simulation model to both end of comp domain.
            else
            {
                sim_param.spacers = (int) ceil((sim_param.Nz/2)/sim_param.dz);
                sim_param.Nz += sim_param.Nz;
                
            }
            cout << "Model cells + spacers: " << sim_param.Nz << endl;
            //Adjust Nz to be divisible by 2
            while(fmod(sim_param.Nz,2.0) != 0)
            {
                //To make Nz divisible by 2, add 1 cell to spacer regions until it becomes divisible by 2
                sim_param.Nz += 1;
            }

            //Check to make sure that the injection point is inside the spacer region. 
            try{
                if(injection_point > sim_param.spacers && injection_point < 0)
                {
                    throw injection_point;
                }
            }
            catch(...)
            {
                cout << "Injection point at cell " << injection_point << " is not valid.";
                return -1;
            }
            
            if(injection_point == 0)
            {
                sim_param.injection_point = (int) ceil(sim_param.spacers/5);
            }
            else
            {
                sim_param.injection_point = injection_point;
            }
            
            
            //Creating the computational domain
            comp_domain.z.setLinSpaced(((((sim_param.Nz*sim_param.dz))-0)/sim_param.dz)+1,0,0+sim_param.dz*(sim_param.Nz-1));
            cout << comp_domain.z << endl;
            
            
            
            
            
            return 0;
        }
        

};
#endif
