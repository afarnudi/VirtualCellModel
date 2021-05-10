//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include "write_functions.hpp"


void writeOutputs(int                atom_count,
                  const MyAtomInfo   all_atoms[],
                  double             time,
                  double             energyInKJ,
                  double             potential_energyInKJ
                  ){
    if (generalParameters.usingBackupCheckpoint) {
        if (generalParameters.WantXYZ) {
            writeXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ,false);
        }
        if (generalParameters.WantVelocity) {
            writeVelocitiesFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ, false);
        }
    } else {
        if (generalParameters.WantXYZ) {
            writeXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ,true);
        }
        if (generalParameters.WantVelocity) {
            writeVelocitiesFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ, true);
        }
    }
}
