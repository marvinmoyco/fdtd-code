# FDTD Codes for ECE 198

### 1-D FDTD
* **FEATURES**
    1. Multiple Boundary Conditions - Dirichlet Boundary Condition or Perfectly Absorbing Boundary Condition.
    2. Multiple Sources - Gaussian Pulse, Sine Wave, Square Wave, Modulated Sine Wave.
    3. Multiple Source Excitation Method - Hard Source, Soft Source, and Total Field/Scatter Field
    4. Basic FDTD Algorithm - Conventional algorithm (leap-frog method)
    5. FDTD-Schwarz Method - Combined algorithm of FDTD-Schwarz in serial and parallel mode (using OpenMP).





##### Notes:
* Libraries used are header-only type so no need to compile them. Just include the directories of the necessary libraries in the makefile before compiling.
* Main code uses C++ for generating simulation data but data visualization uses Python 3 (matplotlib/plotly)

##### Current Issues:

###### FDTD
1. Sinusoidal source propagate backwards (towards the left) when using TF/SF.

###### FDTD-Schwarz Method
1. Sinusoidal source propagate backwards when using TF/SF.
2. Algorithm limited up to 8 subdomains per simulation only because it causes CPU bottleneck if increased to more than 8 subdomains.
3. Reflections generated in the internal boundary of each subdomain. This causes errors in the simulation where some energy is generated out of nowhere.


#### Libraries Used:

**C++ Libraries:**
* [xtensor](https://github.com/xtensor-stack/xtensor) - Main scientific computation library in C++ (very similar to Numpy).
* [xtl](https://github.com/xtensor-stack/xtl) - Dependencies of all xtensor libraries.
* [xtensor-io](https://github.com/xtensor-stack/xtensor-io) - For I/O operations like saving simulation data into HDF5 files.

    * [HgihFive](https://github.com/BlueBrain/HighFive) - C++ library for manipulating HDF5 files (header-only library)

        * **libhdf5** - can be installed using *sudo apt install libhdf5-dev*

* [xtensor-fftw](https://github.com/xtensor-stack/xtensor-fftw.git) - Used in a FFT operation in the simulation.

    * [FFTW 3 Library](http://www.fftw.org/) - can be installed using *sudo apt install libfftw3-dev*

* [xtensor-blas](https://github.com/xtensor-stack/xtensor-blas) - Used in computing L2 norm (using linalg module) 

    * [OpenBlas](https://www.openblas.net/) - can be installed using *sudo apt install libopenblas-dev*

    * [LAPACK](http://www.netlib.org/lapack/) - can be installed using *sudo apt install liblapack-dev*

**Python Libraries:**
* [numpy](https://numpy.org/doc/stable/) - Main numerical library for linear algebra and scientific computing in Python.
* [matplotlib](https://matplotlib.org/stable/contents.html) - Basic plotting library based on MATLAB.
    * [ffmpeg](https://www.ffmpeg.org/) - Required library when using matplotlib for creating animation (video) of simulation
* [plotly](https://plotly.com/python/) - Interactive graphics library for Python.
