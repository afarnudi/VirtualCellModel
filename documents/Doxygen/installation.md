Installation          {#Installation}
============
# OpenMM
The Vertual Cell software package takes advantage of the OpenMM as its main molecular dynamics engine. Please use the instructions on the [OpenMM] webpage to install OpenMM. 

Note: Please install the C++ API. Also note that if you have installed OpenMM with 'Conda' you still need to install the C++ API from the source code (or a install.sh if one is provided by the developer for your OS).During the OpenMM installation you may need to install additional software packages (OpenCL, FFTW,SWIG, etc) depending on your choice of platforms (GPU, CPU) to run the package.

# Generating Mesh files
The Vertual Cell software package can import particle position and bond information from  meshfiles created by [Gmsh] \(export set to version 2\) and [blender]. You can also generate your own mesh information in one of the two supported formats.

# The vertual Cell Package
You can clone the Vertual Cell package from our [github] and use the makefile to build the package. Please note that you may need to modify the make file depending on the openmm installation directory on your machine. 


[OpenMM]: http://openmm.org "OpenMM"
[Gmsh]: http://gmsh.info "Gmsh"
[blender]: https://www.blender.org "blender"
[github]: https://github.com/afarnudi/Membrane_OBJ "github"

