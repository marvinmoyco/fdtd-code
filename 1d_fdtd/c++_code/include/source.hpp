#ifndef SOURCE
#define SOURCE


//#include "common.hpp"



class Source{
    public:
        source_parameters source_param;
        source_output_d source_output;
        computational_domain source_comp_dom;
        

        Source(simulation_parameters sim_param,computational_domain comp_dom)
        {
            source_param.fmax = sim_param.fmax;
            source_param.dz = sim_param.dz;
            source_param.dt = sim_param.dt;
            source_param.Nz = sim_param.Nz;
            source_param.sim_time = sim_param.sim_time;
            source_param.source_type = sim_param.source_type;
            source_comp_dom.n = comp_dom.n;
            source_comp_dom.mu = comp_dom.mu;
            source_comp_dom.epsilon = comp_dom.epsilon;
            source_comp_dom.injection_point = comp_dom.injection_point;
        }

        int GaussianSource(double t0_coeff = 6.0,
                           double prop_coeff = 6.0,
                           double tau_coeff = 12.0,
                           double nmax = 1,
                           double nsrc = 1)
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
            source_output.Hsrc = -sqrt(source_comp_dom.n(source_comp_dom.injection_point))*exp(-pow(source_output.t_H,2));
            
            //cout << "Sizes of the electric and magnetic field component:" << endl;
            //cout << "Esrc: " << source_output.Esrc.size() << " | Hsrc: " << source_output.Hsrc.size() << endl;
            return 0;
        }
        
        int SinusoidalSource(double t0_coeff = 3,
                             double prop_coeff = 3.0,
                             double tau_coeff = 3,
                             double nmax = 1,
                             double nsrc = 1)
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
            view(source_output.Hsrc,range(0,t_condition_size+1)) = -exp(-pow(t_H_initial,2))*(sin(2*numeric_constants<double>::PI*source_param.fmax*t_condition));
            //Computing the electric field and magnetic field component of the source after t0
            view(source_output.Esrc,range(t_condition_size,source_param.Nt)) = (sin(2*numeric_constants<double>::PI*source_param.fmax*view(source_output.t,range(t_condition_size,source_param.Nt))));
            view(source_output.Hsrc,range(t_condition_size,source_param.Nt)) = -sqrt(source_comp_dom.n(source_comp_dom.injection_point))*(sin(2*numeric_constants<double>::PI*source_param.fmax*view(source_output.t,range(t_condition_size,source_param.Nt))));

            return 0;
        }


        int SquareWaveSource(double delay = 0,
                             double width = 0,
                             double nmax = 1,
                             double nsrc = 1)
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
            int end_index = ceil((source_param.tau + source_param.t0/4)/source_param.dt);
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
        
        int ModulatedSineSource(double t0_coeff = 3,
                                double prop_coeff = 3.0,
                                double tau_coeff = 3,
                                double nmax = 1,
                                double nsrc = 1)
        {
            //Calculate the necessary variables
            initialize(t0_coeff,prop_coeff,tau_coeff,nmax);
            cout << "========================================================================" << endl;
            cout << "t0: "<< source_param.t0 << " | tau: " << source_param.tau << " | T: " << source_param.sim_time << endl;

            //Calculate Nt 
            source_param.Nt = ceil(source_param.sim_time/source_param.dt);
            cout << "dt: " << source_param.dt << " seconds" << " | Nt: " << source_param.Nt << " iterations"<< endl;

            
            //source_output.t = arange(0.0,source_param.Nt*source_param.dt,source_param.dt);
            
            source_output.t = linspace<double>(0,source_param.Nt*source_param.dt,source_param.Nt) ;
            //Computing the time input for electric field component
            source_output.t_E = (source_output.t - source_param.t0)/source_param.tau;
            
            //Computing the electric field component of the source
            source_output.Esrc = (sin(2*numeric_constants<double>::PI*source_param.fmax*source_output.t))*(exp(-pow(source_output.t_E,2)));
            

            //Computing the time input for the magnetic field component
            double adj_H = (nsrc*source_param.dz/(2*c_0)) + (source_param.dt/2);
            source_output.t_H = ((source_output.t - source_param.t0)+adj_H)/source_param.tau;
            
            //Computing the magnetic field component of the source
            source_output.Hsrc = -sqrt(source_comp_dom.n(source_comp_dom.injection_point))*(sin(2*numeric_constants<double>::PI*source_param.fmax*source_output.t))*(exp(-pow(source_output.t_H,2)));
            
           
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
            else if(source_param.source_type == "modulatedsine")
            {
                source_param.tau = 0.5/source_param.fmax;
            }
            else{
                cout << "ERROR: Incorrect source type!" << endl;
                return -1;
            }

            source_param.t0 = t0_coeff*source_param.tau;
            source_param.t_prop = (nmax*source_param.Nz*source_param.dz)/c_0;
            //double initial_total_time = tau_coeff*source_param.tau +(prop_coeff*source_param.t_prop);
            double initial_total_time = (prop_coeff*source_param.t_prop);
            //cout << "initial_total: " << initial_total_time << " vs. simparam_sim_time: " << source_param.sim_time << endl;
            //Get the smaller total sim time to save memory
            source_param.sim_time = initial_total_time;



            return 0;
        }


};



#endif