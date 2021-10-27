CXX= g++
xtensor_dir = ../../../libs/xtensor/include
xtl_dir = ../../../libs/xtl/include
#plotlib_dir = ../../../libs/matplotlib-cpp
#Do not comment this line if matplotlib-cpp is used.
python_dir = /usr/include/python3.8
numpy_dir = /usr/local/lib/python3.8/dist-packages/numpy/core/include
STD_VER = c++17
CPPFLAGS= -std=$(STD_VER) -Wall -O0 -I $(xtensor_dir) -I $(xtl_dir)
output_name = main
source = main.cpp
headers = fdtd.hpp

$(output_name).out: $(output_name).o
	$(CXX) $(output_name).o -o $(output_name).out


$(output_name).o: $(source) $(headers)
	$(CXX) $(source) -c $(CPPFLAGS) -o $(output_name).o

.PHONY: clean

clean:
	rm $(output_name).o $(output_name) core