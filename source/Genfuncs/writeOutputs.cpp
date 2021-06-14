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
                  int                frameNum,
                  const MyAtomInfo   all_atoms[],
                  double             time,
                  double             energyInKJ,
                  double             potential_energyInKJ,
                  bool               usingBackupCheckpoint
                  ){
    
    if (generalParameters.WantPDB) {
        myWritePDBFrame(frameNum, time, energyInKJ, potential_energyInKJ, all_atoms, !usingBackupCheckpoint);
    }
    if (generalParameters.WantXYZ) {
        writeXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ, !usingBackupCheckpoint);
    }
    if (generalParameters.WantVelocity) {
        writeVelocitiesFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ, !usingBackupCheckpoint);
    }
    if (generalParameters.WantXYZbin) {
        writeXYZbinFrame(all_atoms, generalParameters.precision, !usingBackupCheckpoint);
    }
    if (generalParameters.WantVelocityBin) {
        writeVELbinFrame(all_atoms, generalParameters.precision, !usingBackupCheckpoint);
    }
    if (generalParameters.WantTPKBin) {
        writeTPKbinFrame(time, energyInKJ, potential_energyInKJ, generalParameters.precision, !usingBackupCheckpoint);
    }
    
}
