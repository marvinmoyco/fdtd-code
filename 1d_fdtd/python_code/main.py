from FDTD.simulation import *

sim1 = Simulation("input.csv")

sim1.init_comp_domain(spacer = 0,
                        inj_point = 0,
                        n_subdom = 4,
                        overlap = 5,
                        multithread = False,
                        algo="fdtd-schwarz")

sim1.simulate(boundary_condition="pabc",excitation_method="tfsf")

sim1.save_sim(name="sample",type="hdf5",output_dir="../../../outputs/",username="marvin",description="sample output data for python version")


sim1.plot_fields()