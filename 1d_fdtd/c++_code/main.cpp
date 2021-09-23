/*
    1D FDTD C++ Code
    
    A simple 1-dimensional FDTD simulator in Free space using Gaussian Pulse as source excitation.

    Based on the lecture videos of Dr. Rumpf (link: https://empossible.net/academics/emp5304/)
    and the book "ELECTROMAGNETIC SIMULATION USING THE FDTD METHOD WITH PYTHON"


*/
#include "fdtd.hpp"

 
int main()
{
  Simulation sim("input.csv");
  sim.create_comp_domain(0,0);
  return 0;
  

}