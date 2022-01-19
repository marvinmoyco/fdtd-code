from FDTD.simulation import *

sim1 = Simulation("input.csv")

sim1.init_comp_domain(spacer = 0,
                        inj_point = 0,
                        n_subdom = 1,
                        overlap = 0,
                        multithread = False,
                        algo="fdtd")

sim1.simulate(boundary_condition="dirichlet",excitation_method="soft")

sim1.save_sim(name="sample",type="hdf5",output_dir="../../../outputs/",username="marvin",description="sample output data for python version")