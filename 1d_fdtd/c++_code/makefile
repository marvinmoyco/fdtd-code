CXX= g++
xtensor_dir = ../../../libs/xtensor/include
xtl_dir = ../../../libs/xtl/include
STD_VER = c++14
CPPFLAGS= -std=$(STD_VER) -Wall -O0 -I $(xtensor_dir) -I $(xtl_dir)
output_name = main
source = main.cpp
headers = fdtd.hpp

$(output_name).out: $(output_name).o
	$(CXX) $(output_name).o -o $(output_name).out


$(output_name).o: $(source) $(headers)
	$(CXX) -c $(CPPFLAGS) $(source) -o $(output_name).o

.PHONY: clean

clean:
	rm $(output_name).o $(output_name) core