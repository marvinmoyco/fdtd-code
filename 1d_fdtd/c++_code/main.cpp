/*
    1D FDTD C++ Code
    
    A simple 1-dimensional FDTD simulator in Free space using Gaussian Pulse as source excitation.

    Based on the lecture videos of Dr. Rumpf (link: https://empossible.net/academics/emp5304/)
    and the book "ELECTROMAGNETIC SIMULATION USING THE FDTD METHOD WITH PYTHON"


*/
#include "fdtd.hpp"
save_data output;

//Arguments (1) file_path of input csv (2) Boundary condition (3) source excitation (4) filename output

int main(int argc, char* argv[])
{
  //cout << " Argc: " << argc << endl;
  //Simulation sim("./csv/input/input.csv");
  Simulation sim(argv[1]);
  computational_domain check = sim.create_comp_domain(0,0,1000);
  if(check.check != -1)
  {
    //output = sim.simulate("pabc","tfsf");
    //sim.save_to_file("marvin");
    output = sim.simulate(argv[2],argv[3]);
    sim.save_to_file(argv[4],argv[5]);
    return 0;
  }


  return 0;
  

}