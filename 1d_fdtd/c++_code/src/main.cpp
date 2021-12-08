/*
    1D FDTD C++ Code
    
    A simple 1-dimensional FDTD simulator in Free space using Gaussian Pulse as source excitation.

    Based on the lecture videos of Dr. Rumpf (link: https://empossible.net/academics/emp5304/)
    and the book "ELECTROMAGNETIC SIMULATION USING THE FDTD METHOD WITH PYTHON"


    ffmpeg command: ffmpeg -f image2 -framerate 200 -i E_H_FFT_images_%07d.jpeg -s 1920x1080 output.mp4

*/
#include "simulation.hpp"
save_data output;

/*
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
  
  //Create a simulation object...
  Simulation sim(argv[1]);

  //Perform pre-processing...
  computational_domain check = sim.create_comp_domain(0,0,100,stoul(argv[8]),stod(argv[10]),argv[9],argv[7]);

  //Start the simulation...
  output = sim.simulate(argv[2],argv[3]);

  //Save the simulation data....
  sim.save_to_file(argv[4],argv[5],argv[6],argv[11],argv[12],argv[13]);




  return 0;
  

}