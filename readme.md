# FDTD Codes for ECE 198

### 1-D FDTD
* Plane wave propagation simulation using 1-D FDTD 
    -The current code can simulate plane wave propagation by injectin soft source and no absorbing boundary conditions.





##### Notes:
* Libraries used are header-only type so no need to compile them. Just include the directories of the necessary libraries in the makefile before compiling.
* Migrated the FDTD Code to C++ but still currently using matplotlib for the plotting.

##### Current Issues:

###### FDTD
1. Sinusoidal source propagate backwards (towards the left).
2. Grid dispersion (I think) is present since weird oscillations in the grid is observed to start at the end of the computational domain. (SOLVED)
3. PABC still reflects some energy and do propagate all of the energy outside the computational domain.
4. Integration of gnuplot to the current code is still underway. (Scrapped. Will continue to use matplotlib to plot the data.)

#### Libraries Used:

**C++ Libraries:**
* [xtensor](https://github.com/xtensor-stack/xtensor)
* [xtl](https://github.com/xtensor-stack/xtl)

**Python Libraries:**
* [numpy](https://numpy.org/doc/stable/)
* [matplotlib](https://matplotlib.org/stable/contents.html)