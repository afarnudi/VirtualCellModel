#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "write_functions.hpp"

void minimisation(MyOpenMMData*              omm,
                  MyAtomInfo                 all_atoms[],
                  Bonds*                     all_bonds
                  ) {
    cout<<"Minimising the coordinates before running the simulation."<<endl;
    cout<<"Minimisation parameters:"<<endl;
    cout<<"\tMinimisation Tolerance: "<<generalParameters.MinimiseTolerance<<endl;
    cout<<"\tMinimisation Max Iterations: "<<generalParameters.MinimiseMaxIterations<<endl<<endl;
    OpenMM::LocalEnergyMinimizer::minimize(*(omm->context));
    double fake_time, fake_energyInKJ, fake_potential_energyInKJ;
    myGetOpenMMState(omm->context, fake_time, fake_energyInKJ, fake_potential_energyInKJ, all_atoms);
    myWritePDBFrame(0, 0, 0, 0, all_atoms, all_bonds);
    string traj_name= generalParameters.trajectory_file_name+"_positionsAfterPushoff.xyz";
    ofstream writexyz(traj_name.c_str());
    for (int n=0; all_atoms[n].type != -1; n++) {
        writexyz<<all_atoms[n].posInNm[0]<<"\t"<<all_atoms[n].posInNm[1]<<"\t"<<all_atoms[n].posInNm[2]<<"\n";
    }
    writexyz.close();
//    exit(0);
    cout<<"\tMinimisation finished. New coordinates written to PDB."<<endl<<endl;
    
}
