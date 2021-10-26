#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "OpenMM_WCA.hpp"

//const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;

double get_WCA_CutoffDistance(double    radius_1,
                              double    radius_2,
                              bool      use_max_radius){
//    1.1224620483 * (radius_1+radius_2)
    double cutoff;
    if (use_max_radius) {
        cutoff = pow(2, 1./6) * 2 * max( radius_1, radius_2);
    } else {
        cutoff = pow(2, 1./6) *( radius_1 + radius_2 );
    }
    return cutoff;
}

double get_WCA_sigma(double    radius_1,
                     double    radius_2,
                     bool      use_max_radius){
    double sigma;
    if (use_max_radius) {
        sigma =  2 * max( radius_1, radius_2);
    } else {
        sigma = ( radius_1 + radius_2 );
    }
    return sigma;
}


string generate_WCA_potential(int       set_1_index,
                              int       set_2_index,
                              string    set_1_name,
                              string    set_2_name){
    
    string sigma   = generate_parameter_name("sigma", set_1_index, set_2_index, set_1_name,  set_2_name);
    string epsilon = to_string(4*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature) ;
    string potential = epsilon+"*(("+sigma+"/r)^12-("+sigma+"/r)^6+0.25)";
    return potential;
}

string generate_WCA_potential(string    set_1_name,
                              string    set_2_name){
    
    string sigma   = generate_parameter_name("sigma", set_1_name,  set_2_name);
    string epsilon = to_string(4*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature) ;
    string potential = epsilon+"*(("+sigma+"/r)^12-("+sigma+"/r)^6+0.25)";
    return potential;
}


void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                      atoms[],
                          set<int>                              set_1,
                          set<int>                              set_2,
                          string                                set_1_name,
                          string                                set_2_name,
                          bool                                  use_max_radius){

    double atom1_radius = atoms[*set_1.begin()].radius;
    double atom2_radius = atoms[*set_2.begin()].radius;
    
    string sigma   = generate_parameter_name("sigma", set_1_name,  set_2_name);
    string potential = generate_WCA_potential( set_1_name, set_2_name);
    
    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));
    int index = int(WCAs.size())-1;

    WCAs[index]->addGlobalParameter(sigma, get_WCA_sigma(atom1_radius, atom2_radius, use_max_radius) );

    if (generalParameters.Periodic_condtion_status) {
        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
    } else {
        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }

    WCAs[index]-> setCutoffDistance( get_WCA_CutoffDistance(atom1_radius, atom2_radius, use_max_radius) );
                                    
    WCAs[index]-> addInteractionGroup(set_1, set_2);
}




//void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
//                          const MyAtomInfo                      atoms[],
//                          vector<set<int> >                     vec_set_1,
//                          vector<set<int> >                     vec_set_2,
//                          int                                   set_1_index,
//                          int                                   set_2_index,
//                          string                                set_1_name,
//                          string                                set_2_name){
//
//    set<int> set_1 = vec_set_1[set_1_index];
//    set<int> set_2 = vec_set_2[set_2_index];
//    
//    double atom1_radius = atoms[*set_1.begin()].radius;
//    double atom2_radius = atoms[*set_2.begin()].radius;
//    
//    string sigma   = generate_parameter_name("sigma", set_1_index, set_2_index, set_1_name,  set_2_name);
//    string potential = generate_WCA_potential(set_1_index, set_2_index, set_1_name, set_2_name);
//    
//    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));
//    int index = int(WCAs.size())-1;
//
//    WCAs[index]->addGlobalParameter(sigma, get_WCA_sigma(atom1_radius, atom2_radius, false) );
//
//    if (generalParameters.Periodic_condtion_status) {
//        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
//    } else {
//        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
//    }
//
//    WCAs[index]-> setCutoffDistance( get_WCA_CutoffDistance(atom1_radius, atom2_radius, false) );
//                                    
//    WCAs[index]-> addInteractionGroup(set_1, set_2);
//}
//
//
//
//void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
//                          const MyAtomInfo                       atoms[],
//                          vector<vector<set<int> > >             vec_vec_set_1,
//                          vector<set<int> >                      vec_set_2,
//                          int                                    set_1_index,
//                          int                                    set_2_index,
//                          string                                 set_1_name,
//                          string                                 set_2_name,
//                          bool                                   use_max_radius){
//
//    set<int> set_1 = get_flat_set(vec_vec_set_1, set_1_index);
//    set<int> set_2 = vec_set_2[set_2_index];
//    
//    double atom1_radius = atoms[*set_1.begin()].radius;
//    double atom2_radius = atoms[*set_2.begin()].radius;
//    
//    string sigma   = generate_parameter_name("sigma", set_1_index, set_2_index, set_1_name,  set_2_name);
//    string potential = generate_WCA_potential(set_1_index, set_2_index, set_1_name, set_2_name);
//    
//    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));
//    int index = int(WCAs.size())-1;
//    
//    WCAs[index]->addGlobalParameter(sigma,   get_WCA_sigma(atom1_radius, atom2_radius, use_max_radius) );
//    
//    
//    if (generalParameters.Periodic_condtion_status) {
//        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
//    } else {
//        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
//    }
//    
//    WCAs[index]-> setCutoffDistance( get_WCA_CutoffDistance(atom1_radius, atom2_radius, use_max_radius) );
//    
//    WCAs[index]-> addInteractionGroup(set_1, set_2);
//}
////
//void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
//                          const MyAtomInfo                       atoms[],
//                          vector<vector<set<int> > >             vec_vec_set_1,
//                          vector<vector<set<int> > >             vec_vec_set_2,
//                          int                                    set_1_index,
//                          int                                    set_2_index,
//                          string                                 set_1_name,
//                          string                                 set_2_name){
//
//
//    set<int> set_1 = get_flat_set(vec_vec_set_1, set_1_index);
//    set<int> set_2 = get_flat_set(vec_vec_set_2, set_2_index);
//    
//    double atom1_radius = atoms[*set_1.begin()].radius;
//    double atom2_radius = atoms[*set_2.begin()].radius;
//    
//    string sigma   = generate_parameter_name("sigma", set_1_index, set_2_index, set_1_name,  set_2_name);
//    string potential = generate_WCA_potential(set_1_index, set_2_index, set_1_name, set_2_name);
//    
//    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));
//    int index = int(WCAs.size())-1;
//    
//    WCAs[index]->addGlobalParameter(sigma,   get_WCA_sigma(atom1_radius, atom2_radius, false) );
//    
//    if (generalParameters.Periodic_condtion_status) {
//        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
//    } else {
//        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
//    }
//    
//    WCAs[index]-> setCutoffDistance( get_WCA_CutoffDistance(atom1_radius, atom2_radius, false) );
//    
//    WCAs[index]-> addInteractionGroup(set_1, set_2);
//}


















