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
/**Relay the position information of the membrane nodes to other data structures ready to pass to OpenMM handles.*/
MyAtomInfo* convert_membrane_position_to_openmm(Membrane mem);
/**Relay the bond information of the membrane nodes to other data structures ready to pass to OpenMM handles.*/
Bonds* convert_membrane_bond_info_to_openmm(Membrane mem);
/**Relay the dihedral angle (triangle-triangle angle) information of the membrane triangle to other data structures ready to pass to OpenMM handles.*/
Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane mem);

#endif
