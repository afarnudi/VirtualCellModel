## The Virtual Cell Model ##

The Virtual Cell Model (VCM) software package is a Molecular Dynamics (MD) simulation software. The core concept of VCM is to get particle coordinates and the potentials between them and solve Newton's equations to calculate the dynamics of the system. 

VCM's design focuses on simulating linear chains of particles (such as long polymers), particles on 2D triangular meshes (vesicles), and particles on 3D Meshes (solid objects). The combination of these inputs can be used to build any general system for simulation porpuses. With OpenMM as the engine, VCM can fully utilise the resources of the host machine (Multiple CPU cores and/or GPUs).
Please visit the [VCM] website for more information.

# Download and install
Required softwares:
[openMM]
[Boost]

Step 1: Clone.
Step 2: Change directory to Membrane_OBJ and make the project.
```console
cd Membrane_OBJ
make -j4
```

Please visit our [installation] page for more information.

[installation]: https://afarnudi.github.io/Membrane_OBJ/md__doxygen__installation.html
[VCM]: https://afarnudi.github.io/Membrane_OBJ/index.html
[openMM] http://openmm.org "OpenMM"
[Boost]https://www.boost.org
