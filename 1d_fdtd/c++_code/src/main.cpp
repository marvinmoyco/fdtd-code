/*
    1D FDTD C++ Code
    
    A simple 1-dimensional FDTD simulator in Free space using Gaussian Pulse as source excitation.

    Based on the lecture videos of Dr. Rumpf (link: https://empossible.net/academics/emp5304/)
    and the book "ELECTROMAGNETIC SIMULATION USING THE FDTD METHOD WITH PYTHON"


    ffmpeg command: ffmpeg -f image2 -framerate 200 -i E_H_FFT_images_%07d.jpeg -s 1920x1080 output.mp4

*/
#include "simulation.hpp"
save_data output;

//Arguments (1) file_path of input csv (2) Boundary condition (3) source excitation (4) custom name for output (5) output file type (6) output directory

int main(int argc, char* argv[])
{
  
  //Create a simulation object...
  Simulation sim(argv[1]);

  //Perform pre-processing...
  computational_domain check = sim.create_comp_domain(0,0,100,1,false);

  //Start the simulation...
  output = sim.simulate(argv[2],argv[3]);

  //Save the simulation data....
  sim.save_to_file(argv[4],argv[5],argv[6]);




  return 0;
  

}