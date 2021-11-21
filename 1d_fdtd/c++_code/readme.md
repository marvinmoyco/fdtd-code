#Instructions to simulate in a terminal

1. Make sure that the makefile is modified so that the libraries and header files are properly searched by the program.
2. To run the program, you have two options: (1) run the compiled program located at ./src/bin or (2) use the bash script "run_fdtd.sh".
3. The compiled program has the following input arguments (in this order):
    1. File path (Absolute/Relative path) of the input file (in csv format).
    2. Boundary conditions to be used in the simulation ('dirichlet' - Dirichlet Boundary Conditions or 'pabc' - Perfectly Absorbing Boundary Condition for 1D simulations only).
    3. Source Excitation method to be used in the simulation ('hard' - Hard Source Method, 'soft' - Soft Source Method, and 'tfsf' - Total Field/Scatter Field Method).
    4. Custom name of the output file (any string of characters will suffice).
    5. Output file type ('csv' - csv files of the simulation data (least comprehensive), 'npy' - Numpy array file (a 3D array composed of simulation data), and 'hdf5' - Hierarchial Data Format 5 (Similar to Python dictionaries with key-value pairs which makes this the most comprehensive output file)).
    6. Output directory for the output file of the simulation (Abosulte/Relative path).
4. To use the bash script, make sure that the following variables are set appropriately:
    1. fdtd_out - Name of the compiled C++ program.
    2. fdtd_bin - Path of the compiled C++ program (Absolute/Relative path).
    3. input_file - Name of the input csv file.
    4. input_dir - Path of the input csv file (Absolute/Relative).
    5. output_dir - Path of the output file (csv/npy/hdf5).
5. The output of the C++ program are composed of 2 files:
    1. Log file - txt file containing the cout commands of the C++ program and the terminal outputs of the Python file.
    2. Output file - Main output file of the simulation that contains the data from the program.