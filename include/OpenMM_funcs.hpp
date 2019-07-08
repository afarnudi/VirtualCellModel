//
//  General_constants.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 26/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef OpenMM_funcs_hpp
#define OpenMM_funcs_hpp

#include "OpenMM_structs.h"
#include "Membrane.h"


/** This function and an opaque structure are used to interface our main
 * programme with OpenMM without the main programme having any direct interaction
 * with the OpenMM API.
 * This function initilises the openmm system + contex + forces.
 */
MyOpenMMData* myInitializeOpenMM(const MyAtomInfo               atoms[],
                                 double                         stepSizeInFs,
                                 std::string&                   platformName,
                                 Bonds*                         bonds,
                                 Dihedrals*                     dihedrals,
                                 std::vector<std::set<int> >    &membrane_set,
                                 std::vector<std::vector<int> > interaction_map);




/** -----------------------------------------------------------------------------
 *                     TAKE MULTIPLE STEPS USING OpenMM
 * -----------------------------------------------------------------------------
 */
void          myStepWithOpenMM(MyOpenMMData*, int numSteps);

/** -----------------------------------------------------------------------------
 *                     COPY STATE BACK TO CPU FROM OPENMM
 * -----------------------------------------------------------------------------
 */
void          myGetOpenMMState(MyOpenMMData*, bool wantEnergy,
                               double& time, double& energy,
                               MyAtomInfo atoms[]);
/** -----------------------------------------------------------------------------
 *                     DEALLOCATE OpenMM OBJECTS
 * -----------------------------------------------------------------------------
 */
void          myTerminateOpenMM(MyOpenMMData*);

/**                               PDB FILE WRITER
 * Given state data, output a single frame (pdb "model") of the trajectory.
 */
void myWritePDBFrame(int frameNum, double timeInPs, double energyInKcal,
                const MyAtomInfo atoms[], std::string traj_name);

/**Relay the position information of the membrane nodes to other data structures ready to pass to OpenMM handles.*/
MyAtomInfo* convert_membrane_position_to_openmm(Membrane mem);
/**Relay the bond information of the membrane nodes to other data structures ready to pass to OpenMM handles.*/
Bonds* convert_membrane_bond_info_to_openmm(Membrane mem);
/**Relay the dihedral angle (triangle-triangle angle) information of the membrane triangle to other data structures ready to pass to OpenMM handles.*/
Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane &mem);

#endif
