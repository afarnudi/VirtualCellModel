#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

//const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;

//void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
//                          const MyAtomInfo                      atoms[],
//                          vector<set<int> >                     set_1,
//                          vector<set<int> >                     set_2,
//                          int                                   set_1_index,
//                          int                                   set_2_index,
//                          string                                set_1_name,
//                          string                                set_2_name){
//
//    set<int> :: iterator it_1 = set_1[set_1_index].begin();
//    set<int> :: iterator it_2 = set_2[set_2_index].begin();
//
//
//
//    string sigma = "sigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
//    //string potential = "(" + sigma + "/(r-" + sigma + "))^12";
//    string potential = "100*((" + sigma + "/r)^12 - (" + sigma + "/r)^6 + 0.25)";
//    ExcludedVolumes.push_back(new OpenMM::CustomNonbondedForce(potential));
//    int index = ExcludedVolumes.size()-1;
//
////    if(atoms[*it_1].radius >= atoms[*it_2].radius )
////    {
////        ExcludedVolumes[index]->addGlobalParameter(sigma,    atoms[*it_1].radius);
////    }
////    else
////    {
////        ExcludedVolumes[index]->addGlobalParameter(sigma,    atoms[*it_2].radius);
////    }
//    ExcludedVolumes[index]->addGlobalParameter(sigma,    0.5*(atoms[*it_1].radius+atoms[*it_2].radius));
//
//    if (generalParameters.Periodic_condtion_status) {
//        ExcludedVolumes[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
//    } else {
//        ExcludedVolumes[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
//    }
//
//    if(atoms[*it_1].radius >= atoms[*it_2].radius )
//    {
//        ExcludedVolumes[index]-> setCutoffDistance( 1.1224 * atoms[*it_1].radius);
//    }
//    else
//    {
//        ExcludedVolumes[index]-> setCutoffDistance( 1.1224 * atoms[*it_2].radius);
//    }
//
//
//    ExcludedVolumes[index]-> addInteractionGroup(set_1[set_1_index], set_2[set_2_index]);
//}



//void init_Modified_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
//                                      const MyAtomInfo                      atoms[],
//                                      vector<set<int> >                     set_1,
//                                      vector<set<int> >                     set_2,
//                                      int                                   set_1_index,
//                                      int                                   set_2_index,
//                                      string                                set_1_name,
//                                      string                                set_2_name){
//
//    set<int> :: iterator it_1 = set_1[set_1_index].begin();
//    set<int> :: iterator it_2 = set_2[set_2_index].begin();
//
//
//
//    string sigma = "sigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
//    string potential = "10*(" + sigma + "/r)^12";
//
//    ExcludedVolumes.push_back(new OpenMM::CustomNonbondedForce(potential));
//    int index = ExcludedVolumes.size()-1;
//    ExcludedVolumes[index]->addGlobalParameter(sigma,   0.5*( atoms[*it_1].radius
//                                                        + atoms[*it_2].radius ) );
////    ExcludedVolumes[index]->addPerParticleParameter("sigma");
////    cout<<"sigma = 0.5*("<<atoms[*it_1].radius<<" + "<<atoms[*it_2].radius<<") = "<<0.5*( atoms[*it_1].radius
////                                                                                       + atoms[*it_2].radius )<<endl;
//    if (generalParameters.Periodic_condtion_status) {
//        ExcludedVolumes[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
//    } else {
//        ExcludedVolumes[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
//    }
//
//    ExcludedVolumes[index]-> setCutoffDistance( 2.5 * ( atoms[*it_1].radius
//                                                        + atoms[*it_2].radius ) );
//
//    ExcludedVolumes[index]-> addInteractionGroup(set_1[set_1_index], set_2[set_2_index]);
//}
//
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

    string epsilon = "WCAepsilon" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string sigma = "WCAsigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string cutoff = "WCAcutoff" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string potential = "4*"+epsilon+"*(("+sigma+"/r)^12-("+sigma+"/r)^6+0.25)*step("+cutoff+"-r)";

    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));


    int index = WCAs.size()-1;
    if (use_max_radius) {
        WCAs[index]->addGlobalParameter(sigma,   max( atoms[*it_1].radius
                                                    , atoms[*it_2].radius ) );
    } else {
        WCAs[index]->addGlobalParameter(sigma,       ( atoms[*it_1].radius
                                                     + atoms[*it_2].radius ) );
    }
    WCAs[index]->addGlobalParameter(epsilon,   generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature);

    if (generalParameters.Periodic_condtion_status) {
        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
    } else {
        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }

    if (use_max_radius) {
        WCAs[index]-> setCutoffDistance( pow(2, 1./6) * max( atoms[*it_1].radius
                                                          , atoms[*it_2].radius ) );
        WCAs[index]->addGlobalParameter(cutoff,   pow(2, 1./6) * max( atoms[*it_1].radius
                                                                   , atoms[*it_2].radius ) );
    } else {
        WCAs[index]-> setCutoffDistance( pow(2, 1./6)  *( atoms[*it_1].radius
                                                        + atoms[*it_2].radius ) );
        WCAs[index]->addGlobalParameter(cutoff,   pow(2, 1./6) *( atoms[*it_1].radius
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


    string epsilon = "WCAepsilon" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string sigma = "WCAsigma" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string cutoff = "WCAcutoff" + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    string potential = "4*"+epsilon+"*(("+sigma+"/r)^12-("+sigma+"/r)^6+0.25)*step("+cutoff+"-r)";

    WCAs.push_back(new OpenMM::CustomNonbondedForce(potential));


    int index = WCAs.size()-1;
    WCAs[index]->addGlobalParameter(sigma,   ( atoms[*it_1].radius
                                             + atoms[*it_2].radius ));
    WCAs[index]->addGlobalParameter(epsilon,   generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature);
    if (generalParameters.Periodic_condtion_status) {
        WCAs[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
    } else {
        WCAs[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }
    WCAs[index]-> setCutoffDistance( pow(2,1./6) *( atoms[*it_1].radius
                                                  + atoms[*it_2].radius ) );
    WCAs[index]->addGlobalParameter(cutoff,   pow(2, 1./6) *( atoms[*it_1].radius
                                                            + atoms[*it_2].radius ) );
    WCAs[index]-> addInteractionGroup(compined_set_1, compined_set_2);
}
//
//
//void init_Excluded_volume_interaction(vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
//                                      const MyAtomInfo                       atoms[],
//                                      vector<vector<set<int> > >             set_1,
//                                      vector<vector<set<int> > >             set_2,
//                                      int                                    set_1_index,
//                                      int                                    set_2_index,
//                                      int                                    sub_set_1,
//                                      int                                    sub_set_2,
//                                      string                                 set_1_name,
//                                      string                                 set_2_name){
//
//
//    set<int> :: iterator it_1 = set_1[set_1_index][sub_set_1].begin();
//    set<int> :: iterator it_2 = set_2[set_2_index][sub_set_2].begin();
//
//
//    string sigma = "sigma" + set_1_name + std::to_string(set_1_index) + std::to_string(sub_set_1) + set_2_name + std::to_string(set_2_index) + std::to_string(sub_set_2);
//    string potential = "(" + sigma + "/(r-" + sigma + "))^12";
//
//    ExcludedVolumes.push_back(new OpenMM::CustomNonbondedForce(potential));
//
//
//    int index = ExcludedVolumes.size()-1;
//    ExcludedVolumes[index]->addGlobalParameter(sigma,   0.5*( atoms[*it_1].radius
//                                                                 + atoms[*it_2].radius ) );
//    if (generalParameters.Periodic_condtion_status) {
//        ExcludedVolumes[index]-> setNonbondedMethod( OpenMM::CustomNonbondedForce::CutoffPeriodic);
//    } else {
//        ExcludedVolumes[index]-> setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
//    }
//    ExcludedVolumes[index]-> setCutoffDistance( 2.5 * ( atoms[*it_1].radius
//                                                         + atoms[*it_2].radius ) );
//
//    ExcludedVolumes[index]-> addInteractionGroup(set_1[set_1_index][sub_set_1], set_2[set_2_index][sub_set_2]);
//}
