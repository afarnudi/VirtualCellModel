#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

//const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;

void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                      atoms[],
                          vector<set<int> >                     set_1,
                          vector<set<int> >                     set_2,
                          int                                   set_1_index,
                          int                                   set_2_index,
                          string                                set_1_name,
                          string                                set_2_name){

    set<int> :: iterator it_1 = set_1[set_1_index].begin();
    set<int> :: iterator it_2 = set_2[set_2_index].begin();



    string sigma   = "sigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string epsilon = to_string(4*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature) ;
    string potential = epsilon+"*(("+sigma+"/r)^12-("+sigma+"/r)^6+0.25)";
    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));
    int index = WCAs.size()-1;

//    if(atoms[*it_1].radius >= atoms[*it_2].radius )
//    {
//        ExcludedVolumes[index]->addGlobalParameter(sigma,    atoms[*it_1].radius);
//    }
//    else
//    {
//        ExcludedVolumes[index]->addGlobalParameter(sigma,    atoms[*it_2].radius);
//    }
    WCAs[index]->addGlobalParameter(sigma,      atoms[*it_1].radius+atoms[*it_2].radius );

    if (generalParameters.Periodic_condtion_status) {
        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
    } else {
        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }

    WCAs[index]-> setCutoffDistance( 1.1224620483 * (atoms[*it_1].radius+atoms[*it_2].radius) );
    
     


    WCAs[index]-> addInteractionGroup(set_1[set_1_index], set_2[set_2_index]);
}



void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                       atoms[],
                          vector<vector<set<int> > >             set_1,
                          vector<set<int> >                      set_2,
                          int                                    set_1_index,
                          int                                    set_2_index,
                          string                                 set_1_name,
                          string                                 set_2_name,
                          bool                                   use_max_radius){

    set<int> compined_set;
    for (int i=0; i<set_1[set_1_index].size(); i++) {
        compined_set.insert(set_1[set_1_index][i].begin(),set_1[set_1_index][i].end());
    }
    
    set<int> :: iterator it_1 = compined_set.begin();
    set<int> :: iterator it_2 = set_2[set_2_index].begin();
    
    string sigma   = "sigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string epsilon = to_string(4*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature) ;
    string potential = epsilon+"*(("+sigma+"/r)^12-("+sigma+"/r)^6+0.25)";
    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));
    int index = WCAs.size()-1;
    
    if (use_max_radius) {
        WCAs[index]->addGlobalParameter(sigma,   2*max( atoms[*it_1].radius
                                                                 , atoms[*it_2].radius ) );
    } else {
        WCAs[index]->addGlobalParameter(sigma,   ( atoms[*it_1].radius
                                                                 + atoms[*it_2].radius ) );
    }
    
    
    if (generalParameters.Periodic_condtion_status) {
        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
    } else {
        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }
    
    if (use_max_radius) {
        WCAs[index]-> setCutoffDistance( pow(2, 1./6) * 2 * max( atoms[*it_1].radius
                                                          , atoms[*it_2].radius ) );
    } else {
        WCAs[index]-> setCutoffDistance( pow(2, 1./6) * 2 *( atoms[*it_1].radius
                                                        + atoms[*it_2].radius ) );
    }
    
    WCAs[index]-> addInteractionGroup(compined_set, set_2[set_2_index]);
}
//
void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                       atoms[],
                          vector<vector<set<int> > >             set_1,
                          vector<vector<set<int> > >             set_2,
                          int                                    set_1_index,
                          int                                    set_2_index,
                          string                                 set_1_name,
                          string                                 set_2_name){


    set<int> compined_set_1;
    set<int> compined_set_2;
    
    for (int i=0; i<set_1[set_1_index].size(); i++) {
        compined_set_1.insert(set_1[set_1_index][i].begin(),set_1[set_1_index][i].end());
    }
    for (int i=0; i<set_2[set_2_index].size(); i++) {
        compined_set_2.insert(set_1[set_2_index][i].begin(),set_1[set_2_index][i].end());
    }
    
    set<int> :: iterator it_1 = compined_set_1.begin();
    set<int> :: iterator it_2 = compined_set_2.begin();
    
    
    string sigma   = "sigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string epsilon = to_string(4*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature) ;
    string potential = epsilon+"*(("+sigma+"/r)^12-("+sigma+"/r)^6+0.25)";
    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));
    
    
    int index = WCAs.size()-1;
    WCAs[index]->addGlobalParameter(sigma,   ( atoms[*it_1].radius
                                                                 + atoms[*it_2].radius ));
    
    if (generalParameters.Periodic_condtion_status) {
        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
    } else {
        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }
    WCAs[index]-> setCutoffDistance( pow(2, 1./6) * ( atoms[*it_1].radius
                                                         + atoms[*it_2].radius ) );
    
    WCAs[index]-> addInteractionGroup(compined_set_1, compined_set_2);
}


















