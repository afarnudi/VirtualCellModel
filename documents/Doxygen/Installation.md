# Installing VCM #

VCM is a software package written in C++. It currently runs on MacOS and Linux. To build it you would need a c++ compiler like (g++, llvm, etc). VCM uses OpenMM as its main engin. If you don't have OpenMM please read the next section.

STEP 1: Clone the Membrane_OBJ project from our [github] page.

STEP 2: Change directory to Membrane_OBJ.

STEP 3: Run makefile by typing 'make'. You may need to modify the make file depending on the openmm installation directory on your machine (Default: /usr/local/openmm/).

## OpenMM ##
The Vertual Cell software package takes advantage of the OpenMM as its main molecular dynamics engine. Please use the instructions on the [OpenMM] webpage to install OpenMM. 

Note: Please install the C++ API. Also note that if you have already installed OpenMM with 'Conda' on your system you still need to install the C++ API from the source code (or a install.sh if one is provided by the developer for your OS). During the OpenMM installation you may need to install additional software packages (OpenCL, FFTW,SWIG, etc) depending on your choice of platforms (GPU, CPU) to run the package.

You can clone the Vertual Cell package from our [github] and use the makefile to build the package. Please note that you may need to modify the make file depending on the openmm installation directory on your machine. 


[OpenMM]: http://openmm.org "OpenMM"
[Gmsh]: http://gmsh.info "Gmsh"
[blender]: https://www.blender.org "blender"
[github]: https://github.com/afarnudi/Membrane_OBJ "github"

