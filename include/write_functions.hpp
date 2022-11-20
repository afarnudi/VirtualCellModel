//
//  write_functions.hpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#ifndef write_functions_hpp
#define write_functions_hpp

#include <stdio.h>
#include "ECM.h"
#include "Membrane.h"
#include <iomanip>
#include <limits>
#include "OpenMM_structs.h"

struct xyzStashInfo{
    std::string xyzframes;
    std::string velframes;
    std::string xyzbinframes;
    std::string velbinframes;
    std::string tpkbinframes;
};



void Results (ECM ecm, string label, char* buffer);
void Results (Membrane membrane, string label, char* buffer);
void ParticleGeneratingReport (char* buffer, Membrane particle );
void MembraneGeneratingReport (char* buffer, Membrane membrane );
void EcmGeneratingReport (char* buffer, ECM ecm );
void checkingForce (Membrane membrane, int MD_Step, char* buffer);

void writeOutputs(int                atom_count,
                  int                frameNum,
                  const MyAtomInfo   atoms[],
                  double             time,
                  double             energyInKJ,
                  double             potential_energyInKJ,
                  bool               usingBackupCheckpoint
                 );
void appendOutputs(int                atom_count,
                  int                frameNum,
                  const MyAtomInfo   atoms[],
                  double             time,
                  double             energyInKJ,
                  double             potential_energyInKJ,
                   xyzStashInfo     &xyzStash
                 );
void flushOutputs(xyzStashInfo      xyzStash);

/**                               PDB FILE WRITER
 * Given state data, output a single frame (pdb "model") of the trajectory.
 */
void myWritePDBFrame(int                frameNum,
                     double             timeInPs,
                     double             energyInKcal,
                     double             potential_energy,
                     const MyAtomInfo   atoms[],
                     bool               copyFromBuffer);
/**                               PSF FILE WRITER
 * Given state data at the beginning of the simulation, output a protein structure file.
 */
void myWritePSF(int   num_of_atoms,
                int   num_of_bonds,
                const MyAtomInfo   atoms[],
                const Bonds        bonds[]);



/**                               PDB FILE WRITER
 * Given state data, output a single frame (pdb "model") of the trajectory.
 */
void writeXYZFrame  (int                atom_count,
                     const MyAtomInfo   atoms[],
                     double             time,
                     double             energyInKJ,
                     double             potential_energyInKJ,
                     bool buff
                     );
std::string stashXYZFrame  (int                atom_count,
                     const MyAtomInfo   atoms[],
                     double             time,
                     double             energyInKJ,
                     double             potential_energyInKJ
                     );

void writeVelocitiesFrame  (int                atom_count,
                            const MyAtomInfo   atoms[],
                            double             time,
                            double             energyInKJ,
                            double             potential_energyInKJ,
                            bool buff
                            );
std::string stashVelocitiesFrame  (int                atom_count,
                            const MyAtomInfo   atoms[],
                            double             time,
                            double             energyInKJ,
                            double             potential_energyInKJ
                            );

void writeXYZbinFrame(const MyAtomInfo atoms[],
                      string precision,
                      bool copyFromBuffer);
std::string  stashXYZbinFrame(const MyAtomInfo atoms[],
                              string precision);
void writeVELbinFrame(const MyAtomInfo atoms[],
                      string precision,
                      bool copyFromBuffer);
std::string stashVELbinFrame(const MyAtomInfo atoms[],
                      string precision);
void writeTPKbinFrame(double             time,
                      double             energyInKJ,
                      double             potential_energyInKJ,
                      string precision,
                      bool copyFromBuffer);
std::string stashTPKbinFrame(double             time,
                      double             energyInKJ,
                      double             potential_energyInKJ,
                      string precision);

#endif /* write_functions_hpp */
