#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;

void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<set<int> >                     set_1,
                              vector<set<int> >                     set_2,
                              int                                   set_1_index,
                              int                                   set_2_index){
    
    set<int> :: iterator it_1 = set_1[set_1_index].begin();
    set<int> :: iterator it_2 = set_2[set_2_index].begin();
    
    LJ_12_6_interactions.push_back(new OpenMM::CustomNonbondedForce("epsilon*((sigma/r)^12-2*(sigma/r)^6)"));
    int index = LJ_12_6_interactions.size()-1;
    
    
    
    LJ_12_6_interactions[index]->addGlobalParameter("sigma",   0.5*( atoms[*it_1].sigma_LJ_12_6
                                                                    + atoms[*it_2].sigma_LJ_12_6 )
                                                    * OpenMM::NmPerAngstrom);
    
    
    LJ_12_6_interactions[index]->addGlobalParameter("epsilon",  sqrt(atoms[*it_1].epsilon_LJ_12_6
                                                                     * atoms[*it_2].epsilon_LJ_12_6)
                                                    * OpenMM::KJPerKcal
                                                    * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    LJ_12_6_interactions[index]->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    LJ_12_6_interactions[index]->setCutoffDistance(1.5 * ( atoms[*it_1].sigma_LJ_12_6
                                                          + atoms[*it_2].sigma_LJ_12_6 )
                                                   * OpenMM::NmPerAngstrom);
    
    
    
    LJ_12_6_interactions[index]->addInteractionGroup(set_1[set_1_index], set_2[set_2_index]);
    
}
