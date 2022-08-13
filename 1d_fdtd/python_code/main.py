from FDTD.simulation import *
from FDTD.plot import *


sim1 = Simulation("input.csv")

sim1.init_comp_domain(spacer = 0,
                        inj_point = 0,
                        n_subdom = 1,
                        overlap = 5,
                        multithread = False,
                        algo="fdtd")

sim1.simulate(boundary_condition="pabc",excitation_method="tfsf")

plot1 = Plot(simulation=sim1,type='html',n_frame=2,save=False,output_path='../../../outputs/',read=False,output="sample_output")


plot1.plot_html()

#sim1.save_sim(name="sample",type="hdf5",output_dir="../../../outputs/",username="marvin",description="sample output data for python version")


