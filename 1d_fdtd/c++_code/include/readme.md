# Structure of the different header files in FDTD-Schwartz Library

The different classes are divided into different files to prevent one long file and for easier debugging.

* **common.hpp** -  As the file says, it is a common file that is utilized by all header files. It contains all of the header dependencies and struct and global variable declarations.
* **simulation.hpp** - Main header file that is included in the main.cpp. This header file contains the Simulation class and calls upon the other header files. Note: This is the only header that is needed to be included in the C++ file (cpp).
* **source.hpp** - This header file contains the Source class. All of the computations regarding create source excitation in FDTD simulation are done here.
* **subdomain.hpp** - This header file contains the Subdomain class. All of the methods and processes that are done by each subdomain (except for the main simulate method) are done here.