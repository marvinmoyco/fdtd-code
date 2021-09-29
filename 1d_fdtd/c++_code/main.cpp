/*
    1D FDTD C++ Code
    
    A simple 1-dimensional FDTD simulator in Free space using Gaussian Pulse as source excitation.

    Based on the lecture videos of Dr. Rumpf (link: https://empossible.net/academics/emp5304/)
    and the book "ELECTROMAGNETIC SIMULATION USING THE FDTD METHOD WITH PYTHON"


*/
#include "fdtd.hpp"
save_data output;
 
int main()
{
  Simulation sim("./csv/input/input.csv");
  int check = sim.create_comp_domain(0,0);
  if(check != -1)
  {
    output = sim.simulate("pabc","hard");
    sim.save_to_file("marvin");
    return 0;
  }
  return 0;
  

}