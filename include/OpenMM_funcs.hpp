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
#include "ECM.h"
#include "Actin.h"
#include "Chromatin.h"


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
                                 std::vector<std::set<int> >    &actin_set,
                                 std::vector<std::set<int> >    &ecm_set,
                                 std::vector<std::set<int> >    &chromatin_set,
                                 std::vector<std::vector<int> > interaction_map);




/** -----------------------------------------------------------------------------
 *                     TAKE MULTIPLE STEPS USING OpenMM
 * -----------------------------------------------------------------------------
 */
void          myStepWithOpenMM(MyOpenMMData*,
                               int numSteps);

/** -----------------------------------------------------------------------------
 *                     COPY STATE BACK TO CPU FROM OPENMM
 * -----------------------------------------------------------------------------
 */
void          myGetOpenMMState(MyOpenMMData*,
                               bool         wantEnergy,
                               bool         wantForce,
                               double&      time,
                               double&      energy,
                               MyAtomInfo   atoms[]);
/** -----------------------------------------------------------------------------
 *                     DEALLOCATE OpenMM OBJECTS
 * -----------------------------------------------------------------------------
 */
void          myTerminateOpenMM(MyOpenMMData*);

/**
 * Calculate the energy for the membrane bacteria problem using the surface equation.
 */
void calc_energy(vector<Membrane>     mem,
                 MyAtomInfo           atoms[]);

/**
 * Calculate the energy for the membrane bacteria problem usinf perpendicular vectors.
 */
void calc_energy_2(vector<Membrane>     mem,
                   MyAtomInfo           atoms[]);



/**                               PDB FILE WRITER
 * Given state data, output a single frame (pdb "model") of the trajectory.
 */
void myWritePDBFrame(int                frameNum,
                     bool               wantforce,
                     double             timeInPs,
                     double             energyInKcal,
                     const MyAtomInfo   atoms[],
                     std::string        traj_name);

/**Relay Membrane class's atom information to other data structures ready to pass to OpenMM handles.*/
void OpenMM_membrane_info_relay (vector<Membrane>       membranes,
                                 vector<std::set<int> > &membrane_set,
                                 MyAtomInfo*            all_atoms,
                                 Bonds*                 all_bonds,
                                 Dihedrals*             all_dihedrals,
                                 int                    &atom_count,
                                 int                    &bond_count,
                                 int                    &dihe_count);
/**Relay the position information of the membrane nodes to other data structures ready to pass to OpenMM handles.*/
MyAtomInfo* convert_membrane_position_to_openmm(Membrane mem);
/**Relay the bond information of the membrane nodes to other data structures ready to pass to OpenMM handles.*/
Bonds* convert_membrane_bond_info_to_openmm(Membrane mem);
/**Relay the dihedral angle (triangle-triangle angle) information of the membrane triangle to other data structures ready to pass to OpenMM handles.*/
Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane &mem);

/**Relay Actin class's atom information to other data structures ready to pass to OpenMM handles.*/
void OpenMM_Actin_info_relay (vector<Actin>          acts,
                              vector<std::set<int> > &act_set,
                              MyAtomInfo*            all_atoms,
                              Bonds*                 all_bonds,
                              Dihedrals*             all_dihedrals,
                              int                    &atom_count,
                              int                    &bond_count,
                              int                    &dihe_coun);
/**Relay the position information of the Actin nodes to other data structures ready to pass to OpenMM handles.*/
MyAtomInfo* convert_Actin_position_to_openmm(Actin act);
/**Relay the bond information of the Actin nodes to other data structures ready to pass to OpenMM handles.*/
Bonds* convert_Actin_bond_info_to_openmm(Actin act);


/**Relay ECM class's atom information to other data structures ready to pass to OpenMM handles.*/
void OpenMM_ECM_info_relay (vector<ECM>             ecms,
                            vector<std::set<int> >  &ecm_set,
                            MyAtomInfo*             all_atoms,
                            Bonds*                  all_bonds,
                            Dihedrals*              all_dihedrals,
                            int                     &atom_count,
                            int                     &bond_count,
                            int                     &dihe_count);
/**Relay the position information of the ECM nodes to other data structures ready to pass to OpenMM handles.*/
MyAtomInfo* convert_ECM_position_to_openmm(ECM ecm);
/**Relay the bond information of the ECM nodes to other data structures ready to pass to OpenMM handles.*/
Bonds* convert_ECM_bond_info_to_openmm(ECM ecm);

void OpenMM_Chromatin_info_relay (vector<Chromatin>         chromos,
                                  vector<std::set<int> >    &chromo_set,
                                  MyAtomInfo*               all_atoms,
                                  Bonds*                    all_bonds,
                                  Dihedrals*                all_dihedrals,
                                  int                       &atom_count,
                                  int                       &bond_count,
                                  int                       &dihe_coun);
/**Relay the position information of the Actin nodes to other data structures ready to pass to OpenMM handles.*/
MyAtomInfo* convert_Chromatin_position_to_openmm(Chromatin chromo);
/**Relay the bond information of the Actin nodes to other data structures ready to pass to OpenMM handles.*/
Bonds* convert_Chromatin_bond_info_to_openmm(Chromatin chromo);



/**Creates a vector of node pairs (excluded bonds) from the system bond list using class labels.*/
std::vector< std::pair< int, int > > exclusion_list_generator(Bonds*      bonds,
                                                              std::string label_1,
                                                              std::string label_2);

/**Initiate the Lenard Jones 12 6 interaction for sets of class atoms.*/
void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<std::set<int> >                set_1,
                              vector<std::set<int> >                set_2,
                              int                                   set_1_index,
                              int                                   set_2_index);

void init_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                      const MyAtomInfo                      atoms[],
                                      vector<std::set<int> >                set_1,
                                      vector<std::set<int> >                set_2,
                                      int                                   set_1_index,
                                      int                                   set_2_index);

#endif
