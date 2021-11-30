#ifndef SIMULATION
#define SIMULATION
 


#include "common.hpp"
#include "source.hpp"
#include "subdomain.hpp"


//Class declarations
class Simulation
{

    public:
        //Initialization of the struct data types
        simulation_parameters sim_param;
        computational_domain comp_domain;
        input_data input;
        fields sim_fields;
        Source* sim_source;
        source_output_d sim_source_fields;
        save_data output;
        vector<Subdomain> subdomains;

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
                    if(temp_l_size(i) == 0){
                        continue;
                    }
                    else{
                        cout << "Layer " << i+1 << ": \t";
                        input.layer_size(i) = temp_l_size(i);
                        input.magnetic_permeability(i)= temp_mu(i);
                        input.electric_permittivity(i) = temp_epsilon(i);
                        cout << input.layer_size(i) << " \t \t" 
                            << input.magnetic_permeability(i) << " \t \t" 
                            << input.electric_permittivity(i) << endl;
                    }
                    
                    
                }

            }   

            
            
        }
        

        //Creating computational domain
        computational_domain create_comp_domain(int spacer_cells = 0,
                                                int injection_point = 0, 
                                                double n_freq = 0, 
                                                unsigned int num_subdomains=1,
                                                double overlap=0, 
                                                bool multithread = false, 
                                                string algorithm = "fdtd-schwarz")
        {
            /*
                This function is the "pre-processing" stage of the simulation. All of the needed pre-requisite computations are done before the actual
                simulation is done. The simulate() can be used after this function successfully finished.
                Input arguments:
                1. spacer_cells - in meters. amount of spacing between end of comp domain and the device model (will be adjusted in the computations)
                2. injection_point - Position (in cell index) of where the source is injected in the comp domain. This will always be within the spacer region.
                3. n_freq - Number of frequency points to be used in the Frequency Response.
                4. num_subdomains - Number of subdomains created in FDTD-Schwarz Algorithm.
                5. multithread - flag to determine whether the algorithm used is serial or parallel.
                6. overlap - amount of overlap used in Schwarz method. Range (0,1) 0% to 100%.
            */

            
            //Try catch here to make sure that the input struct is not empty
            cout << "========================================================================" << endl;
            try
            {
                bool check_input = input.layer_size.size() != 0 && 
                                   input.magnetic_permeability.size() != 0  && 
                                   input.electric_permittivity.size() != 0  && 
                                   !input.simulation_parameters.empty();
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

            //Store the number of subdomains in sim_param struct....
            sim_param.num_subdomains = num_subdomains;

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
                sim_param.left_spacers = ceil(spacer_cells/sim_param.dz);
                sim_param.right_spacers = ceil(spacer_cells/sim_param.dz);
            }
            //If there is not, add half of the size of the simulation model to both end of comp domain.
            else
            {
                sim_param.left_spacers = (int) ceil((sim_param.Nz/2));
                sim_param.right_spacers = (int) ceil((sim_param.Nz/2));
                sim_param.Nz += sim_param.Nz;
                
            }
            //Adjust Nz to be divisible by the numbe of subdomains
            while((sim_param.Nz % sim_param.num_subdomains) != 0)
            {
                
                //To make Nz divisible by number of subdomains, add 1 cell to the left spacer region until it becomes divisible by it
                sim_param.Nz += 1;
                sim_param.left_spacers +=1;
            }

            //Check to make sure that the injection point is inside the spacer region. 
            if(injection_point < 0)
            {
                cout << "Error detected: Injection point is invalid" << endl;
                comp_domain.check = -1;
                exit(EXIT_FAILURE);
            }
            else if(injection_point > sim_param.left_spacers)
            {
                cout << "Error detected: Injection point is inside the device model" <<endl;
                comp_domain.check = -1;
                exit(EXIT_FAILURE);
            }
            else{
                if(injection_point == 0)
                {
                   
                    sim_param.injection_point = (int) ceil(sim_param.left_spacers/2);
                }
                else
                {
                    sim_param.injection_point = injection_point;
                }
            }
            cout << "Injection point (index/position): " << sim_param.injection_point << "-th cell" <<endl;
            cout << "Spacer cells (Left): " << sim_param.left_spacers << " cells" 
                 << "Spacer cells (Right): " << sim_param.right_spacers << " cells" << endl;
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
            int start = sim_param.left_spacers;
            int end = sim_param.left_spacers;
            for (int i=0;i<input.simulation_parameters.at(2);i++)
            {
                end += model_ncells(i);
                //cout << "range[" << start << ":" << end << "]" << "mu value: " << input.magnetic_permeability(i) << endl;
                view(comp_domain.mu,range(start,end))= input.magnetic_permeability(i);
                view(comp_domain.epsilon,range(start,end)) = input.electric_permittivity(i);
                
                start = end;
            }


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
            cout << "|---" << sim_param.left_spacers << " cells";
            for(int i=0;i< input.simulation_parameters.at(2);i++)
            {
                cout << "---" << model_ncells(i) << " cells";
            }
            cout << "---" << sim_param.right_spacers << " cells---|" << endl;

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

            //Get the proper source type based on the input file...
            switch((int) input.simulation_parameters.at(1))
            {
                case 0:
                    sim_param.source_type = "gaussian";
                    break;
                case 1:
                    sim_param.source_type = "sinusoidal";
                    break;
                case 2:
                    sim_param.source_type = "square";
                    break;
                case 3:
                    sim_param.source_type = "modulatedsine";
                    break;
                default:
                    cout << "ERROR: Invalid source type!";
                    comp_domain.check = -1;
                    exit(EXIT_FAILURE);
                    break;
            }
            //Compute the source that will be used in the project.
            compute_source();
            
            /*
            At this point, the code needs to check if it is in serial or parallel mode. 
            Meaning, all of the pre-processing for both serial and parallel versions are done in this method.
            The only data transferred to the subdomain class are the simulation parameters and the
            computational domain vectors
            
            */
            unsigned long int row = sim_param.Nt;
            unsigned long int col = sim_param.Nz;
            sim_param.multithread = multithread;
            sim_param.algorithm = algorithm;
            if(sim_param.algorithm == "fdtd")
            {
                //FDTD Basic Version (Serial)

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

                
                //Resize the different save matrices
                output.E.resize({row,col});
                output.H.resize({row,col});

                //Print shape
                //for (auto& el : output.E.shape()) {cout << el << ", ";}
                //cout << output.E.size() <<endl;

                //cout << comp_domain.z.size() << "   " << comp_domain.mu.size() << endl;

                
            }
            else if(sim_param.algorithm == "fdtd-schwarz")
            {
                /*
                In this part, do the following: (this assumes that the multithreading part is the FDTD-Schwarz Serial Version)
                1. Compute the size of each subdomain.
                2. Compute the size of each overlapping region.
                3. Find out which subdomain to inject the source.
                4. call the pre-process subdomain method.
                */

                //Exception Handling to make sure that the number of subdomain and the overlap size is correct.
                //2^x
                
                try
                {
                    //Check to verify that the number of subdomains is valid...
                    cout << "Number of subdomains: " << sim_param.num_subdomains << endl;
                    if(sim_param.num_subdomains % 2 == 0 && sim_param.num_subdomains <= 64) // Divisible by 2 (2^x number of subdomains)
                    {
                        
                        cout << "Proper number of subdomains are detected..." << endl;
                    }
                    else
                    {
                        throw -1;
                    }

              
                }
                catch(...)
                {
                    cout << "Error: Invalid number of subdomain" << endl;
                    comp_domain.check = -1;
                    exit(EXIT_FAILURE);
                }

                /*
                Layout of a subdomain:
                ||--overlapping region--|--subdomain--|--overlapping region--||
                */

                
                
                //Computing the size of overlapping region
                if(overlap > 0 && overlap < 1) //If the overlap input is a percentage...
                {
                    sim_param.overlap_size = sim_param.Nz*overlap; //25% of the original computational domain...
                }
                else //If it is the amount of overlap (in cells)...
                {
                    sim_param.overlap_size = overlap;
                    //Check if the overlap_size is valid (if Nz + overlap is divisible by the number of subdomain)
                    while((sim_param.Nz + sim_param.overlap_size) % sim_param.num_subdomains != 0)
                    {
                        sim_param.overlap_size += 1;
                    }
                    
                }

                //Computing the sizes of subdomain for pre-processing
                sim_param.non_overlap_size = (sim_param.Nz - ((sim_param.num_subdomains -1)*sim_param.overlap_size))/sim_param.num_subdomains;

                /*
                 * The method used here is similar to how a manual STFT is done with 
                 * frame_size = subdomain size
                 * hop_size = how much the 'frame' moves in the computational domain.
                 */  

                unsigned int frame_size = sim_param.non_overlap_size + 2*sim_param.overlap_size;
                unsigned int start = 0;
                unsigned int stop = frame_size;
                unsigned int hop_size = sim_param.non_overlap_size + sim_param.overlap_size;

                cout << "Non_overlap size: " << sim_param.non_overlap_size << " | Overlap size: " << sim_param.overlap_size << endl;
                cout << "Frame size: " << frame_size << " | Hop size: " << hop_size << endl;

                //Pad the computed mu and epsilon vectors with 0.
                auto padded_mu = pad(comp_domain.mu,sim_param.overlap_size,pad_mode::constant,0);
                auto padded_epsilon = pad(comp_domain.epsilon,sim_param.overlap_size,pad_mode::constant,0);

                
                cout << "Padded mu: " << padded_mu << endl
                     << "Padded epsilon: " << padded_epsilon << endl;
                
         
                //Initialize 2D matrices..
                xtensor<double,2> mu_2D;
                xtensor<double,2> epsilon_2D;

                //Get the subsets for each subdomain.
                for(int i=0;i<sim_param.num_subdomains;i++)
                {   
                    cout << "=====================================" << endl;
                    cout << "Start: " << start << " | Stop: " << stop << endl;

                    if(i == 0)
                    {
                        //In the 1st subdomain, the vector is only needed to be inserted to the 2D matrix..
                        mu_2D = atleast_2d(view(padded_mu,range(start,stop)));
                        epsilon_2D = atleast_2d(view(padded_epsilon,range(start,stop)));
                    }
                    else
                    {
                        /*
                        * After the 1st subdomain, we need to stack (vertically) the values inside the frame to create a 2D matrix
                        * where each row is the subdomain while the columns are the cells inside the subdomains.
                        */
                        mu_2D = vstack(xtuple(
                                                mu_2D,
                                                atleast_2d(view(padded_mu,range(start,stop)))
                        ));

                        epsilon_2D = vstack(xtuple(
                                                epsilon_2D,
                                                atleast_2d(view(padded_epsilon,range(start,stop)))
                        ));
                    }

                    //Adjust the indices by the hop size
                    start += hop_size;
                    stop += hop_size;
                   // cout << "Stacked MU "<<  i << ": " << endl << mu_2D << endl;
                }

                cout << "Stacked MU: " << endl << mu_2D << endl
                     << "Stacked EPSILON: " << epsilon_2D << endl;

 
                //Find out which subdomain to to insert the source.
                //At this point, we can just assume that the source will always be injected into the 1st subdomain (center position)

                //Call the preprocess subdomain to create the subdomain objects.
                preprocess_subdomain(mu_2D,epsilon_2D,sim_param.num_subdomains,sim_param.subdomain_size,sim_param.overlap_size);



            }


            //Computing df - frequency step in Freq. Response 
                
            sim_param.df = 1/(sim_param.dt*sim_param.Nt);
            

            //Initialize the frequency vector for Fourier Transform
            sim_param.n_freq = n_freq;
            sim_param.fmax_fft = 0.5*(1/sim_param.dt); //Upper frequency limit Fs/2
            sim_param.n_freq = sim_param.Nt; //Make sure to match the num of points in other FFT
            sim_fields.Freq_range = linspace<double>(0,sim_param.fmax_fft,sim_param.n_freq);
            cout << "Num. of freq. samples: " << sim_param.n_freq << " | Freq_range shape: (" 
                 << sim_fields.Freq_range.shape()[0] << "," << sim_fields.Freq_range.shape()[1] << endl;
            //Initialize the sizes of refl,trans, and kernal vectors
            sim_fields.Kernel_Freq = exp(-1i*2.0*numeric_constants<double>::PI*sim_param.dt*sim_fields.Freq_range);
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
            cout << "df: " << sim_param.df << " | Frequency range: " << sim_fields.Freq_range.size() 
                 << " | Kernel: " << sim_fields.Kernel_Freq.size() << endl;
            cout << "fmax_fft: " << sim_param.fmax_fft << endl;
            cout << "Save Matrices Sizes:" << endl;
            cout << "Reflectance: " << output.Reflectance.size() << " | Transmittance: " 
                 << output.Transmittance.size() << " | Conservation of Energy: " << output.Con_of_Energy.size() << endl;
        
            return comp_domain;
        }
        

        int compute_source()
        {
            //Create a new Source object
            sim_source = new Source(sim_param,comp_domain);
            
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
            else if(sim_param.source_type == "modulatedsine")
            {
                sim_source->ModulatedSineSource(3,12,2,amax(comp_domain.n)(0),comp_domain.n[sim_param.injection_point]);
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


        int preprocess_subdomain(xtensor<double,2> mu_2D,
                                 xtensor<double,2> epsilon_2D,
                                 unsigned int n_subdomain=0,
                                 unsigned int size=0, 
                                 unsigned int overlap_size=0,
                                 unsigned int non_overlap_size =0)
        {
            //Create new subdomain objects based on the number of subdomains
           
            for(int sdomain_index=0;sdomain_index < n_subdomain; sdomain_index++)
            {
                //Create a new object (creating a vector of Subdomain objects)
                subdomains.push_back(Subdomain(sim_param,sim_source_fields,comp_domain,size,overlap_size,non_overlap_size,sdomain_index));
                switch(sdomain_index)
                {
                    case 0:
                        subdomains.front().subdomain_param.source_inject = true;
                        break;
                }
                //Store the computed mu and epsilon for each subdomain (for computation of update coeff)...
                //TO BE TESTED
                subdomains[sdomain_index].subdomain.mu = row(mu_2D,sdomain_index);
                subdomains[sdomain_index].subdomain.epsilon = row(epsilon_2D,sdomain_index);

                //Modify the preprocessed flag
                subdomains[sdomain_index].subdomain_param.preprocessed = true;
                cout << "============================================" << endl;
                cout << "Subdomain " << sdomain_index + 1 << endl;
                cout << "Mu: " << endl << subdomains[sdomain_index].subdomain.mu << endl;
                cout << "Epsilon: " << endl << subdomains[sdomain_index].subdomain.epsilon << endl;
                

                
            }
            //Flag the first and last subdomain...
            
            

            

            return 0;
        }


        //This is a calling function in general to simulate FDTD. It can distinguish between serial and parallel config.
        save_data simulate(string boundary_condition = "", string excitation_method = "")
        {
            //Check if the configuration is serial or parallel...
            if(sim_param.algorithm == "fdtd")
            {
                //call the simulate_serial function..
                auto output_data = simulate_serial(boundary_condition,excitation_method);
            }
            else if(sim_param.algorithm == "fdtd-schwarz")
            {
                //Call the fdtd-schwarz method if the 
                simulate_fdtd_schwarz();
            }
            return output;
        }

        //FDTD ALgorithm only (serial version)
        save_data simulate_serial(string boundary_condition = "", string excitation = "")
        {
            //Checking if the pre-processing is done successfully...
            cout << "========================================================================" << endl;
            cout << "Checking to make sure that the pre-processing is done...." << endl;
            if(comp_domain.check == -1)
            {
                cout << "Pre-processing has failed. Exiting program" << endl;
                exit(EXIT_FAILURE);
            }

            cout << "Pre-processing is completed!" << endl;
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
                view(sim_fields.H,range(0,end_index-1)) = view(sim_fields.H,range(0,end_index-1)) + 
                                                          (view(sim_fields.m_H,range(0,end_index-1)))*(view(sim_fields.E,range(1,end_index)) - 
                                                          view(sim_fields.E,range(0,end_index-1)));

               


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
                    sim_fields.H(comp_domain.injection_point-1) -= (sim_fields.m_H(comp_domain.injection_point-1)*
                                                                    sim_source_fields.Esrc(curr_iteration));
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
                view(sim_fields.E,range(1,end_index)) =  view(sim_fields.E,range(1,end_index)) + 
                                                         (view(sim_fields.m_E,range(1,end_index))*(view(sim_fields.H,range(1,end_index))-
                                                         view(sim_fields.H,range(0,end_index-1))));
 
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
                    sim_fields.E(comp_domain.injection_point) -= (sim_fields.m_E(comp_domain.injection_point)*
                                                                  sim_source_fields.Hsrc(curr_iteration));
                }
                

                //Compute for the Fourier Transform of the current simulation window
                //More cheaper than saving all the data, it will compute at each iteration
                R = R + (pow(sim_fields.Kernel_Freq,curr_iteration)*sim_fields.E(1));
                T = T + (pow(sim_fields.Kernel_Freq,curr_iteration)*sim_fields.E(end_index-1));
                sim_fields.Source_FFT = sim_fields.Source_FFT  + (pow(sim_fields.Kernel_Freq,curr_iteration)*
                                                                  sim_source_fields.Esrc(curr_iteration));

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

            cout << "Computing Frequency Response using FFTW library..." << endl;

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

            //cout << "Compute Frequency Response using FFTW Library..." << source.shape()[0] << ", " << source.shape()[1] << endl;
            //cout << sim_source_fields.Esrc.shape()[0] << ", " << sim_source_fields.Esrc.shape()[1] << endl;
            //cout << "Esrc (Xtensor): " << sim_source_fields.Esrc << endl;
            //cout << "Esrc (Xarray): " << view(source,all(),0).shape()[0] << endl;


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

            //Set the simulation_success flag to true
            output.simulation_success = true;

            return output;
        }

        //FDTD Schwarz Serial Version
        void simulate_fdtd_schwarz(string direction = "right")
        {
            
            if(sim_param.multithread == false)
            {
                /*
                * Part of the code when OpenMP is not implemented (FDTD-Schwarz in Serial configuration)
                * Iterates through every subdomain so concurrent simulation will not be possible.
                */


                //This is the FDTD Time Loop
                for(int curr_iter=0;curr_iter<sim_param.Nt;curr_iter++)
                {
                    bool isConverged = false;
                    //Loop the FDTD-Schwarz until the data has converged...
                    while(isConverged == false)
                    {
                        //Iterate through the subdomain objects...
                        //While loop and dependent on the return value of the convergence function...
                        for(int subdom_index=0;subdom_index<sim_param.num_subdomains;subdom_index++)
                        {
                            //Call the simulate() method of the Subdomain class to proceed to the FDTD Space loop...
                            subdomains[subdom_index].simulate(curr_iter,sim_param.boundary_cond,sim_param.excitation_method);
                        
                            //Transfer the internal boundary data of the adjacent subdomains here....
                            //Toggle here if you want the direction to be left or right in transferring the boundary data...
                            if (direction == "right")
                            {
                                if(subdom_index == sim_param.num_subdomains -1)
                                {
                                    //Skip the last subdomain (since there is no adjacent subdomain in the right side)
                                    continue;
                                }
                                subdomains[subdom_index].transfer_boundary_data(subdomains[subdom_index + 1],direction);
                            }
                            else if(direction == "left")
                            {
                                if(subdom_index == 0)
                                {   
                                    //Skip the first subdomain (since there is no adjacent subdomain in the left side)
                                    continue;
                                }
                                subdomains[subdom_index].transfer_boundary_data(subdomains[subdom_index - 1],direction);
                            }
                            

                        }

                        //Check for convergence here....
                        //isConverged = convergence_function()

                        


                        //Continue the loop if not converged...

                    }

                    //If converged, reconstruct the whole comp domain here...


                    //Update the FFT here....

                }
                
            }
            else if(sim_param.multithread == true)
            {
                /*
                * Part of the code when OpenMP is utilized. For loop is still used but each iteration will be a separate thread.
                */
               cout << "OPENMP PART OF THE CODE...." << endl;
            }
        }

        bool check_convergence()
        {
            bool computed_l2norm = false;
            int counter = 0;
            for(int s_domain=0;s_domain<sim_param.num_subdomains-1;s_domain++)
            {
                computed_l2norm = subdomains[s_domain].compute_L2norm("right");
                
            }

            //Compare final L2 norm sa machine epsilon
            //  numeric_limits<double>::epsilon()
            return true;
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

            //Make sure that the simulation is done successfully...
            if(output.simulation_success == false)
            {
                cout << "Simulation of FDTD is not successful. Exiting program..." << endl;
                exit(EXIT_FAILURE);

            }

            cout << "Saving the simulated data into output file/s...." << endl;
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
                write_to_hdf5(file, string("/Amount of spacing in the left side (number of cells)"), sim_param.left_spacers);
                write_to_hdf5(file, string("/Amount of spacing in the right side (number of cells)"), sim_param.right_spacers);
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
