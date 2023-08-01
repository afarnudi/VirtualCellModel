The Virtual Cell Model
======================

The Virtual Cell Model (VCM) software package is a Molecular Dynamics (MD) simulation software. The core concept of VCM is to get particle coordinates and the potentials between them and solve Newton's equations to calculate the dynamics of the system. 

VCM's design focuses on simulating linear chains of particles (such as long polymers), particles on 2D triangular meshes (vesicles), and particles on 3D Meshes (solid objects). The combination of these inputs can be used to build any general system for simulation purposes. With OpenMM as the engine, VCM can fully utilise the resources of the host machine (Multiple CPU cores and/or GPUs).
Please visit the VCM_ website for more information.

Requirements:
-------------
OpenMM_ (Python installation)

Installation
------------
```
pip install VCM
```
For quick and easy setup, we suggest installation of OpenMM through anaconda and use of anaconda's pip for VCM installation.

Please visit our installation_ page for more information.

.. _installation: https://afarnudi.github.io/VirtualCellModel/md__doxygen__installation.html

.. _VCM: https://afarnudi.github.io/VirtualCellModel/index.html

.. _openMM: http://docs.openmm.org/7.0.0/userguide/application.html