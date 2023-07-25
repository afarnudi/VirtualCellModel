#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "write_functions.hpp"

void minimisation(MyOpenMMData*              omm,
                  MyAtomInfo                 all_atoms[],
                  int                     atom_count
                  ) {
    cout<<"Minimising the coordinates before running the simulation."<<endl;
    cout<<"Minimisation parameters:"<<endl;
    cout<<"\tMinimisation Tolerance: "<<generalParameters.MinimiseTolerance<<endl;
    cout<<"\tMinimisation Max Iterations: "<<generalParameters.MinimiseMaxIterations<<endl<<endl;
    OpenMM::LocalEnergyMinimizer::minimize(*(omm->context), generalParameters.MinimiseTolerance, generalParameters.MinimiseMaxIterations);
    double fake_time, fake_energyInKJ, fake_potential_energyInKJ;
    myGetOpenMMState(omm->context, fake_time, fake_energyInKJ, fake_potential_energyInKJ, all_atoms);
    writeOutputs(atom_count,0,all_atoms,0, fake_energyInKJ, fake_potential_energyInKJ, false);
    cout<<"\tMinimisation finished."<<endl<<endl;
//    exit(0);
}
