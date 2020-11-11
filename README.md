# MC_Simulation

Developed by Christopher Konow and Ziyao Li

This package contains all of the code necessary to run stochastic Monte Carlo simulations of the Survival model for zebrafish skin pattern formation on a growing domain. It is based in the C++ programming language. A detailed description of the algorithm used for the simulations is contained in the STAR Methods section of our manuscript "Influence of survival, promotion, and growth on pattern formation in zebrafish skin". 

To begin a simulation, you must first update the export path in the first line of the code to an existing folder. You can also then change the parameters in the section immediately following it to see different behaviors discussed in the text. Although all of the simulations are performed on a rectangular domain growing in one dimension (axial growth), there are options in the source code files (src) for either trapezoidal growth (similar to that of the zebrafish tailfin) or two-dimensional square growth. However, these functions may require some tinkering of the main.cpp file. Each simulation will output a desired number of .csv files containing the two-dimensional ZebraArray, which is what we consider the zebrafish skin lattice. In addition, at the end of each simulation a "conditions.csv" file will be exported that contain a list of all parameters used in the simulation.

To visualize the exported .csv array files, use the zebrafish_plot package also located in the EpsteinLab Github to convert the arrays into .png images and animations. Note that the zebrafish_plot package runs using Python.

This simulation package runs using C++, and therefore requires an appropriate compiler (we used GNU 7.5.0). In addition, we strongly recommend using CMake to quickly and easily compile the code into the appropriate Makefile. We have included our CMakeLists.txt file in this package for ease of use. To run the simulation, update the folder and parameters to the desired values, make sure CMake is installed, and then run the following lines of code:

    cmake ./
    make
    ./MC_simulation

Note that this can vary based on your setup.

If you have any questions on this code, please contact Chris Konow at ckonow@brandeis.edu.

