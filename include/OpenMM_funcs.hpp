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
                                 TimeDependantData*             time_dependant_data,
                                 Bonds*                         bonds,
                                 Dihedrals*                     dihedrals,
                                 std::vector<std::set<int> >    &membrane_set,
                                 std::vector<std::set<int> >    &actin_set,
                                 std::vector<std::set<int> >    &ecm_set,
//                                 std::vector<std::set<int> >    &chromatin_set,
                                 std::vector<std::vector<std::set<int> >  >    &chromatin_set,
                                 std::vector<std::vector<int> > interaction_map);




/** -----------------------------------------------------------------------------
 *                     TAKE MULTIPLE STEPS USING OpenMM
 * -----------------------------------------------------------------------------
 */
void          myStepWithOpenMM(MyOpenMMData*,
                               TimeDependantData*,
                               MyAtomInfo atoms[],
                               int        numSteps,
                               int&        total_step);

/** -----------------------------------------------------------------------------
 *                     COPY STATE BACK TO CPU FROM OPENMM
 * -----------------------------------------------------------------------------
 */
void          myGetOpenMMState(MyOpenMMData*,
                               double&      time,
                               double&      energy,
                               double&      potential_energy,
                               MyAtomInfo   atoms[]);

/** -----------------------------------------------------------------------------
 *                     COPY STATE(ONLY POSITIONS) BACK TO CPU FROM OPENMM
 * -----------------------------------------------------------------------------
 */
void Cheap_GetOpenMMState(MyOpenMMData*,
                          MyAtomInfo atoms[]);

/** -----------------------------------------------------------------------------
 *                     Update System parameters
 * -----------------------------------------------------------------------------
 */
void Kelvin_Voigt_update(MyOpenMMData*,
                      TimeDependantData*);

void force_update(MyOpenMMData*,
                         TimeDependantData*);

void hill_update(MyOpenMMData*,
TimeDependantData*,  MyAtomInfo atoms[]);

void kf_update(MyOpenMMData*,
TimeDependantData*,  MyAtomInfo atoms[]);


/** -----------------------------------------------------------------------------
 *                     BOND LENGTH
 * -----------------------------------------------------------------------------
 */
std::vector<double> dist_calc(TimeDependantData*,
                              MyAtomInfo atoms[],
                              int        bondtype);


std::vector<double> Nominal_length_calc(TimeDependantData*,
                                        int bondtype);

/** -----------------------------------------------------------------------------
 *                   Monte_Carlo
 * -----------------------------------------------------------------------------
 */
void          Monte_Carlo_Reinitialize(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, Membrane &mem, MyAtomInfo atoms[],int &MC_total_tries,int &Accepted_Try_Counter, double &MC_Acceptance_Rate);

                               
/** -----------------------------------------------------------------------------
 *                     DEALLOCATE OpenMM OBJECTS
 * -----------------------------------------------------------------------------
 */
void          myTerminateOpenMM(MyOpenMMData*,
                                TimeDependantData*);



/**                               PDB FILE WRITER
 * Given state data, output a single frame (pdb "model") of the trajectory.
 */
void myWritePDBFrame(int                frameNum,
                     double             timeInPs,
                     double             energyInKcal,
                     double             potential_energy,
                     const MyAtomInfo   atoms[],
                     const Bonds        bonds[]);

/**                               PDB FILE WRITER
 * Given state data, output a single frame (pdb "model") of the trajectory.
 */
void writeXYZFrame  (int                atom_count,
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
Bonds* convert_Actin_bond_info_to_openmm(Actin act, MyAtomInfo* atoms);


void OpenMM_ActMem_info_relay (vector<Actin>          acts,
                               vector<Membrane>       membranes,
                               Bonds*                 all_bonds,
                               int                    mem_atom_count,
                               int                    &bond_count);


Bonds* convert_ActMem_bond_info_to_openmm(Actin act, int k);



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
Bonds* convert_ECM_bond_info_to_openmm(ECM ecm , MyAtomInfo* atoms);

void OpenMM_Chromatin_info_relay (vector<Chromatin>                 chromos,
//                                  vector<std::set<int> >  &chromo_set,
                                  vector<vector<std::set<int> > >  &chromo_set,
                                  MyAtomInfo*                       all_atoms,
                                  Bonds*                            all_bonds,
                                  Dihedrals*                        all_dihedrals,
                                  int                              &atom_count,
                                  int                              &bond_count,
                                  int                              &dihe_coun);
/**Relay the position information of the Actin nodes to other data structures ready to pass to OpenMM handles.*/
MyAtomInfo* convert_Chromatin_position_to_openmm(Chromatin chromo);
/**Relay the bond information of the Actin nodes to other data structures ready to pass to OpenMM handles.*/
Bonds* convert_Chromatin_bond_info_to_openmm(Chromatin chromo);



/**Creates a vector of node pairs (excluded bonds) from the system bond list using class labels.*/
std::vector< std::pair< int, int > > exclusion_list_generator(Bonds*      bonds);
/**Add excluded bonds to Custom non bonded force from exclude_list.*/
void add_exclusion(OpenMM::CustomNonbondedForce* custom_bond,
                   std::vector< std::pair< int, int > > exclude_list);

/**Initiate the Lenard Jones 12 6 interaction for sets of class atoms.*/
void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<std::set<int> >                set_1,
                              vector<std::set<int> >                set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              string                                set_1_name,
                              string                                set_2_name);
//overload for chromatin class interactions contatining different node types with other classes.
void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<vector<std::set<int> > >       set_1,
                              vector<std::set<int> >                set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              string                                set_1_name,
                              string                                set_2_name);
//overload for chromatin class interactions contatining different node types with other classes.
void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              std::set<int>                         set_1,
                              vector<std::set<int> >                set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              int                                   chromo_ind,
                              string                                set_1_name,
                              string                                set_2_name);
//overload for inter chromatin class interactions contatining different node types.
void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<vector<std::set<int> > >       set_1,
                              vector<vector<std::set<int> > >       set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              int                                   sub_set_1,
                              int                                   sub_set_2,
                              string                                set_1_name,
                              string                                set_2_name);

/**Initiate External force for a set of class atoms.*/
bool init_ext_force(vector<OpenMM::CustomExternalForce*> &ext_force,
                              const MyAtomInfo                      atoms[],
                              vector<std::set<int> >                set_1,
                              int                                   set_1_index,
                              string                                set_name);


/**Initiate the stronger Lenard Jones 12 6 interaction for sets of class atoms.*/
void initdouble_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<std::set<int> >                set_1,
                              vector<std::set<int> >                set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              string                                set_1_name,
                              string                                set_2_name);


/**initiate LJ 4-2*/
void init_LJ_4_2_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
const MyAtomInfo                      atoms[],
vector<std::set<int> >                set_1,
vector<std::set<int> >                set_2,
int                                   set_1_index,
int                                   set_2_index,
string                                set_1_name,
string                                set_2_name);



void init_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                      const MyAtomInfo                      atoms[],
                                      vector<std::set<int> >                set_1,
                                      vector<std::set<int> >                set_2,
                                      int                                   set_1_index,
                                      int                                   set_2_index,
                                      string                                set_1_name,
                                      string                                set_2_name);

void init_Modified_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                      const MyAtomInfo                      atoms[],
                                      vector<std::set<int> >                     set_1,
                                      vector<std::set<int> >                     set_2,
                                      int                                   set_1_index,
                                      int                                   set_2_index,
                                      string                                set_1_name,
                                      string                                set_2_name);
                                      
//overload for chromatin class interactions contatining different node types with other classes.
void init_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                      const MyAtomInfo                      atoms[],
                                      vector<vector<std::set<int> > >       set_1,
                                      vector<std::set<int> >                set_2,
                                      int                                   set_1_index,
                                      int                                   set_2_index,
                                      string                                set_1_name,
                                      string                                set_2_name);
//overload for inter chromatin class interactions contatining different node types.
void init_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                      const MyAtomInfo                       atoms[],
                                      vector<vector<std::set<int> > >             set_1,
                                      vector<vector<std::set<int> > >             set_2,
                                      int                                    set_1_index,
                                      int                                    set_2_index,
                                      string                                 set_1_name,
                                      string                                 set_2_name);

//overload for inter chromatin class interactions contatining different node types.
void init_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                      const MyAtomInfo                       atoms[],
                                      vector<vector<std::set<int> > >        set_1,
                                      vector<vector<std::set<int> > >        set_2,
                                      int                                    set_1_index,
                                      int                                    set_2_index,
                                      int                                    sub_set_1,
                                      int                                    sub_set_2,
                                      string                                 set_1_name,
                                      string                                 set_2_name);

OpenMM::State getCurrentState(MyOpenMMData*  omm,
                              bool        wantEnergy,
                              double&     timeInPs,
                              bool        wantforce);
                        
void setNewState(MyOpenMMData*      omm,
                 bool            wantEnergy,
                 double&         energyInKcal,
                 MyAtomInfo      atoms[]);

/**Set the interaction map inputs to forces in the OpenMM system.*/
void set_interactions(const MyAtomInfo                       atoms[],
                      Bonds*                                 bonds,
                      vector<std::set<int> >                &membrane_set,
                      vector<std::set<int> >                &actin_set,
                      vector<std::set<int> >                &ecm_set,
                      vector<vector<std::set<int> > >       &chromatin_set,
                      vector<vector<int> >                   interaction_map,
                      vector<OpenMM::CustomExternalForce*>  &ext_force,
                      vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                      vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                      OpenMM::System                        &system
                      );
/**Add particles to the system and forces.
 *Specify the atoms and their properties:
 *(1) System needs to know the masses.
 *(2) NonbondedForce needs charges,van der Waals properties (in MD units!).
 *(3) Collect default positions for initializing the simulation later.
 */
void add_particles_to_system_and_forces(const MyAtomInfo                       atoms[],
                                        vector<OpenMM::Vec3>                  &initialPosInNm,
                                        vector<OpenMM::Vec3>                  &initialVelInNmperPs,
                                        vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                                        vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                        OpenMM::System                        &system);
/**Define bonded forces and specify the involved atoms.
 */
void set_bonded_forces(Bonds*                                 bonds,
                       OpenMM::HarmonicBondForce*            &HarmonicBond,
                       OpenMM::HarmonicBondForce*            &Kelvin_VoigtBond,
                       vector<OpenMM::CustomBondForce*>      &X4harmonics,
                       vector<OpenMM::CustomBondForce*>      &FENEs,
                       vector<OpenMM::CustomBondForce*>      &Contractiles,
                       vector<OpenMM::CustomBondForce*>      &KFs,
                       vector<OpenMM::CustomBondForce*>      &hill_bonds,
                       vector<OpenMM::CustomBondForce*>      &Harmonic_minmax,
                       TimeDependantData*                    &time_dependant_data,
                       OpenMM::System                        &system
                       );
/**Set dihedral forces for triangle pair interactions.
 */
void set_dihedral_forces(Dihedrals*                                 dihedrals,
                         vector<OpenMM::CustomCompoundBondForce*>  &DihedralForces,
                         OpenMM::System                            &system
                         );
#endif
