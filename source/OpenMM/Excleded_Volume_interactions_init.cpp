#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

//const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;

void init_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                              const MyAtomInfo                      atoms[],
                              vector<set<int> >                     set_1,
                              vector<set<int> >                     set_2,
                              int                                   set_1_index,
                              int                                   set_2_index){
    
    set<int> :: iterator it_1 = set_1[set_1_index].begin();
    set<int> :: iterator it_2 = set_2[set_2_index].begin();
    
    ExcludedVolumes.push_back(new OpenMM::CustomNonbondedForce("10*(sigma/r)^6"));
    int index = ExcludedVolumes.size()-1;
    ExcludedVolumes[index]->addGlobalParameter("sigma",   0.5*( atoms[*it_1].radius
                                                               + atoms[*it_2].radius )
                                               * OpenMM::NmPerAngstrom);
    ExcludedVolumes[index]-> setNonbondedMethod(    OpenMM::CustomNonbondedForce::CutoffNonPeriodic );
    ExcludedVolumes[index]->setCutoffDistance(1.5* ( atoms[*it_1].radius
                                                    + atoms[*it_2].radius )
                                              * OpenMM::NmPerAngstrom);
    ExcludedVolumes[index]->addInteractionGroup(set_1[set_1_index], set_2[set_2_index]);
}
