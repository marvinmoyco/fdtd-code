/*
    1D FDTD C++ Code
    
    A simple 1-dimensional FDTD simulator in Free space using Gaussian Pulse as source excitation.

    Based on the lecture videos of Dr. Rumpf (link: https://empossible.net/academics/emp5304/)
    and the book "ELECTROMAGNETIC SIMULATION USING THE FDTD METHOD WITH PYTHON"


    ffmpeg command: ffmpeg -f image2 -framerate 200 -i E_H_FFT_images_%07d.jpeg -s 1920x1080 output.mp4

*/
#include "simulation.hpp"
save_data output;

/*hre
* Input Arguments for compiled program:
* (1) file_path of input csv or string "manual" for loading of input parameters
* (2) Boundary condition - "dirichlet" for Dirichlet Boundary Condition or "pabc" for Perfectly Absorbing Boundary Condition
* (3) source excitation - "hard" for Hard Source Method, "soft" for Soft Source Method, and "tfsf" for Total Field/Scatter Field Method
* (4) custom name for output - Any arbitrary string for the filename of the output file
* (5) output file type - "csv" for comma separated values, "npy" for Numpy Array File, and "hdf5" for Hierarchal Data Format 5 file.
* (6) output directory - directory path to store the log file (.log) and the output file (file type depends on (5))
* (7) algorithm - "fdtd" for basic FDTD algorithm or "fdtd-schwarz" for  FDTD-Schwarz Algorithm
* (8) num_subdomains - value greater than 1 and less than 64 to determine how many threads/subdomains will be used
* (9) multithread - true for using openMP and false for serial version
* (10) overlap size - if 0-1, this is a percentage from Nz, if greater than 1, this is number of cells as the overlap
* (11) comprehensive - will the saved data be comprehensive or not? (boolean value) 
* (12) username/name - username or name
* (13) simulation description - short description about the simulation
*
*/
int main(int argc, char* argv[])
{
  auto overall_start = chrono::high_resolution_clock::now();
  //cout << " Argv[1] :  " << argv[1] << endl;
  //Create a simulation object...
  Simulation sim(argv[1]);
  
  //Perform pre-processing...
  computational_domain check = sim.create_comp_domain(0,0,100,stoul(argv[8]),stod(argv[10]),stoul(argv[9]),argv[7]);

  /*for(int i = 0; i<1000; i++)
  {
    if(i % 500 == 0)
    {







      int y = i;
      cout << "Reached " << y << endl;







      
    }
    sim.update_sim_param(i,i+1);
    //printing dz and dt...
    cout << "=================================================" << endl;
    cout << "Dz: [ ";
    for (double i: sim.sim_param.dz_list)
      cout << i << ", ";
    cout << " ]" << endl;
    cout << "Dt: [ ";
    for (double i: sim.sim_param.dt_list)
      cout << i << ", ";
    cout << " ]" << endl;
  }*/


  //Start the simulation...
  auto start_time = chrono::high_resolution_clock::now();
  output = sim.simulate(argv[2],argv[3]);
  auto end_time = chrono::high_resolution_clock::now();

  chrono::duration<double, std::milli> sim_duration = end_time - start_time;
  //Get the algo time in milliseconds
  sim.output.algo_time = sim_duration.count();

  cout << "Electric field shape: (" << output.E.shape()[0] << "," << output.E.shape()[1] << ")" << endl;
  cout << "Magnetic field shape: (" << output.H.shape()[0] << "," << output.H.shape()[1] << ")" << endl;

  auto overall_end = chrono::high_resolution_clock::now();
  chrono::duration<double, std::milli> overall_duration = overall_end - overall_start;

  sim.output.overall_time = overall_duration.count();
  //Save the simulation data....
  sim.save_to_file(argv[4],argv[5],argv[6],argv[11],argv[12],argv[13]);

  cout << "END OF SIMULATION" << endl;



  return 0;
  

}