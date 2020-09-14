## Installing VCM ##

VCM is a software package written in C++. It currently runs on MacOS and Linux. To build it you need a c++ compiler (g++, llvm, etc). VCM uses [OpenMM] as its main engin. If you don't have OpenMM please read the next section. If you already have OpenMM installed, you need to:

STEP 1: Install the [Boost] library

STEP 2: Clone the Membrane_OBJ project from our [github] page.

STEP 3: Change directory to Membrane_OBJ.

STEP 4: Run makefile by typing 'make'. You may need to modify the make file depending on the openmm installation directory on your machine if OpenMM is installed in a differnet directory than /usr/local/openmm/.

STEP 5: Run ./VCM -h for help.

# Update VCM #
To update VCM:

STEP 1: Delete the binary output file

STEP 2: Enter:  
``` console
make clean;make -j4
```

# OpenMM #
The Vertual Cell Model package takes advantage of OpenMM as its main molecular dynamics engine. Please use the instructions on the [OpenMM] webpage to install OpenMM. 

Note: Please install the C++ API. If you have already installed OpenMM with 'Conda' on your system you still need to install the C++ API from the source code (more information on [OpenMM] website). You may need to install additional software packages (OpenCL, FFTW,SWIG, etc) depending on your choice of platforms (GPU, CPU).


[OpenMM]: http://openmm.org "OpenMM"
[Gmsh]: http://gmsh.info "Gmsh"
[blender]: https://www.blender.org "blender"
[github]: https://github.com/afarnudi/Membrane_OBJ "github"
[Boost]: https://www.boost.org
