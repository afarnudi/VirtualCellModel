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
                              int                                   set_2_index,
                              string                                set_1_name,
                              string                                set_2_name){
    
    set<int> :: iterator it_1 = set_1[set_1_index].begin();
    set<int> :: iterator it_2 = set_2[set_2_index].begin();
    
    
    
    string epsilon = "epsilon" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string sigma = "sigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string potential = epsilon + "*((" + sigma + "/r)^12-2*(" + sigma + "/r)^6)" ;
    
    //cout<<potential << '\n' ;
    
    
    LJ_12_6_interactions.push_back(new OpenMM::CustomNonbondedForce(potential));
    int index = LJ_12_6_interactions.size()-1;
    
    
    
    LJ_12_6_interactions[index]->addGlobalParameter(sigma,   0.5*( atoms[*it_1].sigma_LJ_12_6
                                                                    + atoms[*it_2].sigma_LJ_12_6 ) );
    
    
    LJ_12_6_interactions[index]->addGlobalParameter(epsilon,  sqrt(atoms[*it_1].epsilon_LJ_12_6
                                                                     * atoms[*it_2].epsilon_LJ_12_6) );
    LJ_12_6_interactions[index]->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    LJ_12_6_interactions[index]->setCutoffDistance(2.5 * ( atoms[*it_1].sigma_LJ_12_6
                                                          + atoms[*it_2].sigma_LJ_12_6 ) );
    
    
    
    LJ_12_6_interactions[index]->addInteractionGroup(set_1[set_1_index], set_2[set_2_index]);
    
}

void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<vector<set<int> > >            set_1,
                              vector<set<int> >                     set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              string                                set_1_name,
                              string                                set_2_name){
    
    set<int> :: iterator it_1 = set_1[set_1_index][0].begin();
    set<int> :: iterator it_2 = set_2[set_2_index].begin();
    
    set<int> compined_set;
    for (int i=0; i<set_1[set_1_index].size(); i++) {
        compined_set.insert(set_1[set_1_index][i].begin(),set_1[set_1_index][i].end());
    }
    
    string epsilon = "epsilon" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string sigma = "sigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string potential = epsilon + "*((" + sigma + "/r)^12-2*(" + sigma + "/r)^6)" ;
    
    //cout<<potential << '\n' ;
    
    
    LJ_12_6_interactions.push_back(new OpenMM::CustomNonbondedForce(potential));
    int index = LJ_12_6_interactions.size()-1;
    
    
    LJ_12_6_interactions[index]->addGlobalParameter(sigma,   0.5*( atoms[*it_1].sigma_LJ_12_6
                                                                  + atoms[*it_2].sigma_LJ_12_6 ) );
    LJ_12_6_interactions[index]->addGlobalParameter(epsilon,  sqrt(atoms[*it_1].epsilon_LJ_12_6
                                                                 * atoms[*it_2].epsilon_LJ_12_6) );
    LJ_12_6_interactions[index]->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    
    LJ_12_6_interactions[index]->setCutoffDistance(2.5 * ( atoms[*it_1].sigma_LJ_12_6
                                                          + atoms[*it_2].sigma_LJ_12_6 ) );
    
    
    
    LJ_12_6_interactions[index]->addInteractionGroup(compined_set, set_2[set_2_index]);
    
}

void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              set<int>                              set_1,
                              vector<set<int> >                     set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              int                                   chromo_ind,
                              string                                set_1_name,
                              string                                set_2_name){
    
    set<int> :: iterator it_1 = set_1.begin();
    set<int> :: iterator it_2 = set_2[set_2_index].begin();
    
    
    string epsilon = "epsilon" + set_1_name + std::to_string(set_1_index) + std::to_string(chromo_ind) + set_2_name + std::to_string(set_2_index) ;
    string sigma   = "sigma"   + set_1_name + std::to_string(set_1_index) + std::to_string(chromo_ind) + set_2_name + std::to_string(set_2_index) ;
    string potential = epsilon + "*((" + sigma + "/r)^12-2*(" + sigma + "/r)^6)" ;
    
    //cout<<potential << '\n' ;
    
    
    LJ_12_6_interactions.push_back(new OpenMM::CustomNonbondedForce(potential));
    int index = LJ_12_6_interactions.size()-1;
    
    cout<<"sigma   = "<< 0.5*( atoms[*it_1].sigma_LJ_12_6 + atoms[*it_2].sigma_LJ_12_6 ) << endl;
    cout<<"epsilon = "<< sqrt(atoms[*it_1].epsilon_LJ_12_6* atoms[*it_2].epsilon_LJ_12_6) << endl;
    
    LJ_12_6_interactions[index]->addGlobalParameter(sigma,   0.5*( atoms[*it_1].sigma_LJ_12_6 +
                                                                   atoms[*it_2].sigma_LJ_12_6 ) );
    LJ_12_6_interactions[index]->addGlobalParameter(epsilon,  sqrt(atoms[*it_1].epsilon_LJ_12_6
                                                                 * atoms[*it_2].epsilon_LJ_12_6) );
    LJ_12_6_interactions[index]->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    
    LJ_12_6_interactions[index]->setCutoffDistance(2.5 * ( atoms[*it_1].sigma_LJ_12_6
                                                          + atoms[*it_2].sigma_LJ_12_6 ) );
    
    
    
    LJ_12_6_interactions[index]->addInteractionGroup(set_1, set_2[set_2_index]);
    
}

void init_LJ_12_6_interaction(vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                              const MyAtomInfo                      atoms[],
                              vector<vector<set<int> > >            set_1,
                              vector<vector<set<int> > >            set_2,
                              int                                   set_1_index,
                              int                                   set_2_index,
                              int                                   sub_set_1,
                              int                                   sub_set_2,
                              string                                set_1_name,
                              string                                set_2_name){
    
    
    set<int> :: iterator it_1 = set_1[set_1_index][sub_set_1].begin();
    set<int> :: iterator it_2 = set_2[set_2_index][sub_set_2].begin();
    
    string epsilon = "epsilon" + set_1_name + std::to_string(set_1_index) + std::to_string(sub_set_1) + set_2_name + std::to_string(set_2_index) + std::to_string(sub_set_2) ;
    string sigma = "sigma" + set_1_name + std::to_string(set_1_index) + std::to_string(sub_set_1) + set_2_name + std::to_string(set_2_index) + std::to_string(sub_set_2) ;
    string potential = epsilon + "*((" + sigma + "/r)^12-2*(" + sigma + "/r)^6)" ;
    
    LJ_12_6_interactions.push_back(new OpenMM::CustomNonbondedForce(potential));
    
    int index = LJ_12_6_interactions.size()-1;
    
    LJ_12_6_interactions[index]->addGlobalParameter(sigma,   0.5*( atoms[*it_1].sigma_LJ_12_6
                                                                    + atoms[*it_2].sigma_LJ_12_6 ) );
    
    
    LJ_12_6_interactions[index]->addGlobalParameter(epsilon,  sqrt(atoms[*it_1].epsilon_LJ_12_6
                                                                     * atoms[*it_2].epsilon_LJ_12_6) );
    LJ_12_6_interactions[index]->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    LJ_12_6_interactions[index]->setCutoffDistance(2.5 * ( atoms[*it_1].sigma_LJ_12_6
                                                          + atoms[*it_2].sigma_LJ_12_6 ) );
    
    
    
    LJ_12_6_interactions[index]->addInteractionGroup(set_1[set_1_index][sub_set_1], set_2[set_2_index][sub_set_2]);
    
}
















bool init_ext_force(vector<OpenMM::CustomExternalForce*> &ext_force,
                              const MyAtomInfo                      atoms[],
                              vector<set<int> >                     set_1,
                              int                                   set_1_index,
                              string                                set_name){
    bool is_force=false;
    set<int> :: iterator it_1 = set_1[set_1_index].begin();
    int force_model=atoms[*it_1].ext_force_model;
    string kx = "kx" + set_name + std::to_string(set_1_index) ;
    string ky = "ky" + set_name + std::to_string(set_1_index) ;
    string kz = "kz" + set_name + std::to_string(set_1_index) ;
    string potential = "(-1)* " + kx + " *x + " + "(-1)* " + ky + " *y + " + "(-1)* " + kz + " *z ";
    
    switch (force_model) {
        //constant force
        case 1:
        {
            is_force=true;
            ext_force.push_back(new OpenMM::CustomExternalForce(potential));
            
            int index = ext_force.size()-1;
            
            ext_force[index]->addGlobalParameter(kx,  atoms[*it_1].ext_force_constants[0] );
            ext_force[index]->addGlobalParameter(ky,  atoms[*it_1].ext_force_constants[1] );
            ext_force[index]->addGlobalParameter(kz,  atoms[*it_1].ext_force_constants[2] );
            
            for(set<int>::iterator it=set_1[set_1_index].begin(); it != set_1[set_1_index].end(); ++it){
                ext_force[index]->addParticle(*it);
            }
            
        }
            break;
            
        default:
            break;
    }
    return is_force;
}

