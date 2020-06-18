The Virtual Cell Model Software Package           {#mainpage}
============
Breif Introduction
--------------------
The Virtual Cell Model (VCM) software package is a Molecular Dynamics (MD) simulation software that  runs on MacOS and Linux. The core concept of VCM is to get particle coordinates and the potentials between them and solve Newton's equations to calculate the dynamics of the system. VCM's design focuses on simulating linear chains of particles (such as long polymers), particles on 2D triangular meshes (vesicles), and particles on 3D Meshes (solid objects). The combination of these inputs can be used to build any general system for simulation porpuses. With OpenMM as the engine, VCM can fully utilise the resources of the host machine (Multiple CPU cores and/or GPUs). For installation please  proceed to the installation page.

VCM was developed in the [Soft Condensed Matter Group] at [Sharif university of Technology] lead by [Prof Mohammad Reza Ejtehadi]. It is a unifying computational framework to create a multicomponent cell model that has the capability to predict changes in whole cell and cell nucleus characteristics (in terms of shape, direction, and chromatin conformation) on cell substrates.

The VCM utilises 4 basic classes: 1) The Membrane calss sets up particles that lie on 2D triangular meshes, b) The Chroamtin class, is used for 1D linear chains of cosecutive particles, c) The Actin, and d) The ECM both simulate particles on 3D meshe configurations. The Actin class designed to interact more closely with the Membrane class.





[Soft Condensed Matter Group]: http://softmatter.physics.sharif.edu "Soft Condensed Matter Group"
[Sharif university of Technology]: http://www.en.sharif.edu "Sharif English homepage"
[Prof. Mohammad Reza Ejtehadi]: http://sharif.edu/~ejtehadi/ "Prof Ejtehadi's homepage"
[GMSH]: http://gmsh.info "Gmsh homepage"

[^1]: The current version is compatible with the gmsh version II file style. The option is also available in gmsh versions 2 and above.
