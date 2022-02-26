#ifndef SUBDOMAIN
#define SUBDOMAIN

//#include "common.hpp"

class Subdomain
{
    public:
        subdomain_fields s_fields;
        computational_domain subdomain;
        subdomain_data subdomain_param;
        source_output_d subdomain_source;
        save_data_subdomain subdomain_output;
        double left_ghost_cell = 0.0;
        double right_ghost_cell = 0.0;
        
        
        //Initializing error vectors for convergence...
        xtensor<double,1> E_error{0,0};
        xtensor<double,1> H_error{0,0}; 


        //Constructor
        Subdomain(simulation_parameters sim_param, 
                  source_output_d sources,
                  computational_domain domain,
                  unsigned int size, 
                  unsigned int overlap_size, 
                  unsigned int non_overlap_size, 
                  unsigned int id,
                  auto mu,
                  auto epsilon)
        {
            //Transfer necessary data from Simulation class to Subdomain class...
            subdomain_param.subdomain_id = id;
            subdomain_param.subdomain_size = size;
            subdomain_param.overlap = overlap_size;
            subdomain_param.non_overlap_size = non_overlap_size;
            subdomain_param.dz = sim_param.dz;
            subdomain_param.dt = sim_param.dt;
            subdomain_param.Nz = sim_param.Nz;
            subdomain_param.Nt = sim_param.Nt;
            subdomain_param.num_subdomains = sim_param.num_subdomains;

            subdomain_param.left_spacers = sim_param.left_spacers;
            subdomain_param.right_spacers = sim_param.right_spacers;
            subdomain_param.injection_point = sim_param.injection_point;

            //subdomain = domain;
             //Store the computed mu and epsilon for each subdomain (for computation of update coeff)...
            subdomain.mu = mu;
            subdomain.epsilon = epsilon;
            //cout << "Created a new subdomain!" << endl;
            //cout << "Injection point: " << subdomain_param.injection_point << endl;
            //Initializing the fields and vectors
            
            s_fields.E = zeros<double>(subdomain.epsilon.shape());
            s_fields.H = zeros<double>(subdomain.epsilon.shape());
            s_fields.m_E = zeros<double>(subdomain.epsilon.shape());
            s_fields.m_H = zeros<double>(subdomain.epsilon.shape());
            //cout << "Field shapes: " << endl;
            //cout << "E shape: (" <<  s_fields.E.shape()[0] << ",) | H shape: (" << s_fields.H.shape()[0] << ",)" << endl; 
            //s_fields.m_E = (c_0*sim_param.dt)/(comp_domain.epsilon*sim_param.dz);
            //s_fields.m_H = (c_0*sim_param.dt)/(comp_domain.mu*sim_param.dz);

            //Convert the injection point into index of subdomain


            //Compute the update coefficients...
            if(subdomain_param.subdomain_id == 0)
            {
                //Do not include the extra padding at the start..
                view(s_fields.m_E,range(subdomain_param.overlap,_)) = (c_0*subdomain_param.dt)/(view(subdomain.epsilon,range(subdomain_param.overlap,_))*subdomain_param.dz);
                view(s_fields.m_H,range(subdomain_param.overlap,_)) = (c_0*subdomain_param.dt)/(view(subdomain.mu,range(subdomain_param.overlap,_))*subdomain_param.dz);
            }
            else if(subdomain_param.subdomain_id == subdomain_param.num_subdomains - 1)
            {
                //Do not include the extra padding at the end..
                view(s_fields.m_E,range(_,s_fields.m_E.shape()[0] - 1- subdomain_param.overlap)) = (c_0*subdomain_param.dt)/(view(subdomain.epsilon,range(_,s_fields.m_E.shape()[0] - 1
                                                                                                    - subdomain_param.overlap))*subdomain_param.dz);
                view(s_fields.m_H,range(_,s_fields.m_H.shape()[0] - 1- subdomain_param.overlap)) = (c_0*subdomain_param.dt)/(view(subdomain.mu,range(_,s_fields.m_H.shape()[0] - 1
                                                                                                    - subdomain_param.overlap))*subdomain_param.dz);
            
            }
            else
            {
                 
                s_fields.m_E = (c_0*subdomain_param.dt)/(subdomain.epsilon*subdomain_param.dz);
                s_fields.m_H = (c_0*subdomain_param.dt)/(subdomain.mu*subdomain_param.dz);
            
            }

            if(subdomain_param.subdomain_id == subdomain_param.num_subdomains - 1)
            {
                subdomain_source = sources;

                //Convert the injection point index (from the whole comp domain to the subdomain) by adding the amount of padding in the left side.
                subdomain_param.injection_point += subdomain_param.overlap;
            }


            
            

        }


        subdomain_fields simulate(int curr_iteration = 0,string boundary_condition = "", string excitation_method = "")
        {
            /*
            First step in the Schwarz Method - Solving the PDE on the subdomain, in this case the FDTD algorithm space loop.
            */

            //Check whether the pre-processing is successful...
            if(subdomain_param.preprocessed == false)
            {
                cout << "Pre-processing has failed. Exiting program..." << endl;
                exit(EXIT_FAILURE);
            }

            //Check the boundary condition if you are a subdomain at the end...
            subdomain_param.boundary_condition = boundary_condition;
            subdomain_param.excitation_method = excitation_method;
            //cout << "Boundary condition" << subdomain_param.boundary_condition << " | Excitation method: " << subdomain_param.excitation_method << endl;
            
            
            //cout << "Start of the FDTD Space Loop..." << endl;  
            
                
            //Start of the FDTD Space....
        

            //Initialize variable indices
            unsigned int start = 0;
            unsigned int stop = 0;
            if(subdomain_param.subdomain_id == 0) //If the subdomain is the 1st one...
            {
                start = subdomain_param.overlap;
                stop = s_fields.E.shape(0);
            }
            else if(subdomain_param.subdomain_id == subdomain_param.num_subdomains - 1) // If it is the last..
            {
                //cout << "Getting the indices for the last subdomain..." << endl;
                start = 0;
                stop = subdomain_param.non_overlap_size + subdomain_param.overlap -1;
            }
            else{ //If it is in between the 1st and last subdomain...

                start = 0;
                stop = s_fields.E.shape(0);

            }
            //cout << "Indices (start,stop): (" << start << "," << stop << ")" << endl;

            // Step 1: Store boundary data for the 1st subdomain (for the external boundary data) 
            if(subdomain_param.subdomain_id == 0)
            {
                if(boundary_condition == "pabc")
                {
                    s_fields.E(start) = s_fields.H_start.front();
                    s_fields.H_start.pop_front();
                    s_fields.H_start.push_back(s_fields.E(start+1));
                }
                else if(subdomain_param.boundary_condition == "dirichlet")
                {
                    s_fields.E(start) = 0;
                }
                
            }
            else{ 
               
                // Add PABC boundary conditions to LEFT INT boundaries 
                s_fields.E(start) = s_fields.H_start.front();
                s_fields.H_start.pop_front();
                s_fields.H_start.push_back(s_fields.E(start+1));

                 // Use the ghost cells here by updating the leftmost index (0) using the update equation
                s_fields.E(start) = s_fields.E(start) + (s_fields.m_E(start)*(s_fields.H(start) - left_ghost_cell ));
           
            }
           
            // Step 2: Update the H vector from E
            view(s_fields.H,range(start,stop-1)) = view(s_fields.H,range(start,stop-1)) + (view(s_fields.m_H,range(start,stop-1)))*(view(s_fields.E,range(start+1,stop)) - view(s_fields.E,range(start,stop-1)));

            // Step 3: Update source excitation (applicable only when the subdom is the 1st one)
            if(subdomain_param.subdomain_id == subdomain_param.num_subdomains - 1) //Insert the Hsrc when  you are at the 1st subdomain...
            {
                if(subdomain_param.excitation_method == "tfsf")
                {
                    //cout << "Inserting source in H TFSF" << endl;
                    s_fields.H(subdomain_param.injection_point-1) -= (s_fields.m_H(subdomain_param.injection_point-1)*subdomain_source.Esrc(curr_iteration));
                }
            }

            // Step 4: Update E from H
            view(s_fields.E,range(start+1,stop)) =  view(s_fields.E,range(start+1,stop)) + (view(s_fields.m_E,range(start+1,stop))*(view(s_fields.H,range(start+1,stop))-view(s_fields.H,range(start,stop-1))));

            // Step 5: Store H boundary terms
            if(subdomain_param.subdomain_id == subdomain_param.num_subdomains - 1)
            {
                if(subdomain_param.boundary_condition == "pabc")
                {
                    s_fields.H(stop-1) = s_fields.E_end.front();
                    s_fields.E_end.pop_front();
                    s_fields.E_end.push_back(s_fields.H(stop-3));
                    //s_fields.H(stop) = 0;
                }
                else if(subdomain_param.boundary_condition == "dirichlet")
                {
                    s_fields.H(stop-1) = 0;
                }
                
            }
            else
            {
                
                // Add PABC boundary conditions to RIGHT INT boundaries
                s_fields.H(stop-1) = s_fields.E_end.front();
                s_fields.E_end.pop_front();
                s_fields.E_end.push_back(s_fields.H(stop-2));

                // Use the ghost cells here by updating the rightmost index (n) using the update equation
                s_fields.H(stop-1) = s_fields.H(stop-1) + (s_fields.m_H(stop-1)*( right_ghost_cell - s_fields.E(stop-1)));
           
            }

            
            // Step 6: Update E source excitaiton
            if(subdomain_param.subdomain_id == subdomain_param.num_subdomains - 1)
            {
                if(subdomain_param.excitation_method == "hard")
                {
                    //cout << "Inserting source in Hard method" << endl;
                    s_fields.E(subdomain_param.injection_point) = subdomain_source.Esrc(curr_iteration);
                
                }
                else if(subdomain_param.excitation_method == "soft")
                {
                    //cout << "Inserting source in soft method" << endl;
                    s_fields.E(subdomain_param.injection_point) += subdomain_source.Esrc(curr_iteration);
                
                }
                else if(subdomain_param.excitation_method == "tfsf")
                {
                    //cout << "Inserting source in E TFSF" << endl;
                    s_fields.E(subdomain_param.injection_point) -= (s_fields.m_E(subdomain_param.injection_point)*subdomain_source.Hsrc(curr_iteration));
                
                }   
            }
            
           

            // Step 7: Saving the current snapshot in time in the output matrices
       
            if(curr_iteration > 0)
            {
                //If this is not the first iteration, use vstack to stack the data vertically
                subdomain_output.E = vstack(
                                            xtuple(subdomain_output.E,atleast_2d(s_fields.E))
                );
                subdomain_output.H = vstack(
                                            xtuple(subdomain_output.H,atleast_2d(s_fields.H))
                );

              
            }
            else
            {
                //If this is the first iteration, save the data to the xtensor
                subdomain_output.E = atleast_2d(s_fields.E);
                subdomain_output.H = atleast_2d(s_fields.H);
                //In the 1st and last subdomain, get the computed Reflectance and Transmittance respectively
             
            }


            return s_fields;
        }

        bool transfer_boundary_data(Subdomain& adj_subdomain,string side = "right",string method = "old")
        {
            //Intermediary variable for transferring data
            double boundary_data_E[2] = {0.0,0.0}; //1st element = This subdomain, 2nd element = Adjacent subdomain
            double boundary_data_H[2] = {0.0,0.0};
            
            if(side == "left")
            {
                 //Transfer the boundary data ONLY IF the subdomains are not the first one.
                if(subdomain_param.subdomain_id > 0 )
                {
                    if(method == "old")
                    {
                        //for E-fields
                        //Get the boundary data from both subdomains...
                        boundary_data_E[0] = s_fields.E(0);
                        boundary_data_E[1] = adj_subdomain.s_fields.E(adj_subdomain.s_fields.E.size()-1);

                        //Transfer the boundary E-Field data...
                        s_fields.E(0) = adj_subdomain.s_fields.E(adj_subdomain.s_fields.E.size()-adj_subdomain.subdomain_param.overlap);
                        adj_subdomain.s_fields.E(s_fields.E.size()-1) = adj_subdomain.s_fields.E(adj_subdomain.subdomain_param.overlap);
                        

                        s_fields.E(subdomain_param.overlap) = boundary_data_E[1];
                        adj_subdomain.s_fields.E(s_fields.E.size()-adj_subdomain.subdomain_param.overlap) = boundary_data_E[0];
                        

                        //for H-fields
                        //Get the boundary data from both subdomains...
                        boundary_data_H[0] = s_fields.H(0);
                        boundary_data_H[1] = adj_subdomain.s_fields.H(adj_subdomain.s_fields.H.size()-1);

                        //Transfer the boundary E-Field data...
                        s_fields.H(0) = adj_subdomain.s_fields.H(adj_subdomain.s_fields.H.size()-adj_subdomain.subdomain_param.overlap);
                        adj_subdomain.s_fields.H(s_fields.H.size()-1) = adj_subdomain.s_fields.H(adj_subdomain.subdomain_param.overlap);
                        

                        s_fields.H(subdomain_param.overlap) = boundary_data_H[1];
                        adj_subdomain.s_fields.H(s_fields.H.size()-adj_subdomain.subdomain_param.overlap) = boundary_data_H[0];
                    }
                    else if(method == "new")
                    {
                        //boundary_data_E[1] = subdomain_output.E
                    }
                    
                    
                }
                return true;
            }
            else if(side == "right")
            {
                //Transfer the boundary data ONLY IF the subdomains are not the last one.
                if(subdomain_param.subdomain_id < subdomain_param.num_subdomains )
                {
                    if(method == "old")
                    {
                        //for E-fields
                        //Get the boundary data from both subdomains...
                        boundary_data_E[0] = s_fields.E(s_fields.E.size()-1); //Correct
                        boundary_data_E[1] = adj_subdomain.s_fields.E(0);
                        //cout << "For Electric Field data..." << endl;
                        //cout << "Before transferring the data..." << endl;
                        //cout << "Subdomain " + subdomain_param.subdomain_id << " Rightmost value: " << s_fields.E(s_fields.E.size()-1) << endl;
                        //cout << "Subdomain " + adj_subdomain.subdomain_param.subdomain_id << " Leftmost value: " << adj_subdomain.s_fields.E(0) << endl;

                        //Transfer the boundary E-Field data...
                        s_fields.E(s_fields.E.size()-1) = adj_subdomain.s_fields.E(adj_subdomain.subdomain_param.overlap-1);
                        adj_subdomain.s_fields.E(0) = s_fields.E(s_fields.E.size()-subdomain_param.overlap);

                        s_fields.E(s_fields.E.size()-subdomain_param.overlap) = boundary_data_E[1];
                        adj_subdomain.s_fields.E(adj_subdomain.subdomain_param.overlap-1) = boundary_data_E[0];

                        //cout << "After transferring the data..." << endl;
                        //cout << "Subdomain " + subdomain_param.subdomain_id << " Rightmost value: " << s_fields.E(s_fields.E.size()-subdomain_param.overlap) << endl;
                        //cout << "Subdomain " + adj_subdomain.subdomain_param.subdomain_id << " Leftmost value: " <<  adj_subdomain.s_fields.E(adj_subdomain.subdomain_param.overlap) << endl;
                        //for H-fields
                        //Get the boundary data from both subdomains...
                        boundary_data_H[0] = s_fields.H(s_fields.H.size()-1); //Correct
                        boundary_data_H[1] = adj_subdomain.s_fields.H(0);
                        //cout << "For Electric Field data..." << endl;
                        //cout << "Before transferring the data..." << endl;
                        //cout << "Subdomain " + subdomain_param.subdomain_id << " Rightmost value: " << s_fields.E(s_fields.E.size()-1) << endl;
                        //cout << "Subdomain " + adj_subdomain.subdomain_param.subdomain_id << " Leftmost value: " << adj_subdomain.s_fields.E(0) << endl;

                        //Transfer the boundary E-Field data...
                        s_fields.E(s_fields.H.size()-1) = adj_subdomain.s_fields.H(adj_subdomain.subdomain_param.overlap-1);
                        adj_subdomain.s_fields.H(0) = s_fields.H(s_fields.H.size()-subdomain_param.overlap);

                        s_fields.H(s_fields.H.size()-subdomain_param.overlap) = boundary_data_H[1];
                        adj_subdomain.s_fields.H(adj_subdomain.subdomain_param.overlap-1) = boundary_data_H[0];
                    }
                    else if(method == "new")
                    {

                    }
                    

                }
                //cout << "This subdomain: " << s_fields.E << endl;
                //cout << "Adjacent subdomain: " << adj_subdomain.s_fields.E << endl;
                return true;

            }
            else{
                cout << "Incorrect input arguments." << endl;
                return false;
            }
        }


        bool compute_L2norm(string side="")
        {
            /*
                Computes the L2 norm using the linalg module in xtensor...
                using xtensor-BLAS...
            */

            if(side == "left")
            {
                E_error(0) = linalg::norm(view(s_fields.E,range(0,subdomain_param.overlap)),2);
                H_error(0) = linalg::norm(view(s_fields.H,range(0,subdomain_param.overlap)),2);
            }
            else if(side == "right")
            {
                E_error(1) = linalg::norm(view(s_fields.E,range(subdomain_param.overlap+subdomain_param.subdomain_size,s_fields.E.size())),2);
                H_error(1) = linalg::norm(view(s_fields.H,range(subdomain_param.overlap+subdomain_param.subdomain_size,s_fields.H.size())),2);
            }


            return true;
        }

        

        xtensor<double,1> getOverlapValues(char field, string side = "")
        {
            xtensor<double,1> overlapValues = zeros<double>({subdomain_param.overlap});
            if(field == 'E')
            {
                if(side == "left")
                {
                    overlapValues = view(s_fields.E,range(0,subdomain_param.overlap)); 
                }
                else if(side == "right")
                {
                    //cout << "Getting the overlap values now..." << endl;
                   // cout << "Overlap: " << subdomain_param.overlap << " | Subdomain size: " << subdomain_param.subdomain_size << endl;
                    //cout << "Non overlap size: " << subdomain_param.non_overlap_size << " | Efield size: " << s_fields.E.size() << endl;
                    overlapValues = view(s_fields.E,range(subdomain_param.overlap+subdomain_param.non_overlap_size,s_fields.E.size()));
                }
            } 
            else if(field == 'H')
            {
                if(side == "left")
                {
                    overlapValues = view(s_fields.H,range(0,subdomain_param.overlap)); 
                }
                else if(side == "right")
                {
                    overlapValues = view(s_fields.H,range(subdomain_param.overlap+subdomain_param.non_overlap_size,s_fields.H.size()));
                }
            }
            //cout << "Overlap shape: " << overlapValues.shape()[0] << endl;
            return overlapValues; 
        }
};

#endif