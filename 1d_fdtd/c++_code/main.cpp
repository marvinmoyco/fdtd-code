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
  Simulation sim("./csv/input.csv");
  sim.create_comp_domain(0,0);
  output = sim.simulate("dirichlet","no source");
  unsigned long int row_i = 2;
  unsigned long int col = sim.sim_param.Nt+1;
  //Create 2D xtensor for source
  xtensor<double,2> source;
  source.resize({row_i,col});
  row(source,0) = sim.sim_source_fields.Esrc;
  row(source,1) = sim.sim_source_fields.Hsrc;

   //Write to csv the source H component
  ofstream out_source;
  out_source.open("./csv/source.csv");
  dump_csv(out_source,source);

  //Write to csv the E fields
  ofstream out_stream;
  out_stream.open("./csv/electric_field.csv");
  dump_csv(out_stream,output.E);
  //Write to csv the H fields
  ofstream out_stream2;
  out_stream2.open("./csv/magnetic_field.csv");
  dump_csv(out_stream2,output.H);
  return 0;
  

}