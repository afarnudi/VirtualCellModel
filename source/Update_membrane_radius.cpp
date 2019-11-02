#include "General_class_functions.h"
#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;

bool check_for_membrane_update(vector<Membrane>    &membranes,
                               double               time)
{
    bool update_mem = false;
    for (int mem_index=0; mem_index<membranes.size(); mem_index++) {
        if (time < membranes[mem_index].get_End_update_time_in_Ps() &&
            time > membranes[mem_index].get_Begin_update_time_in_Ps()   ) {
            if (membranes[mem_index].get_new_node_radius() != -1) {
                update_mem = true;
            } else if (membranes[mem_index].get_new_nominal_length() != -1){
                update_mem = true;
            }
        }
    }
    return update_mem;
}


double updated_sigam_value(Membrane &mem, double time){
    double t1, t2, r, rnew, dt, a, b;
    // All calculations in Fs
    t1   = mem.get_Begin_update_time_in_Ps()*1000;
    t2   = mem.get_End_update_time_in_Ps()*1000;
    r    = mem.get_node_radius();
    rnew = mem.get_new_node_radius();
    dt   = t2-t1;
    // sigma_ev = a * time + b
    a    = (rnew - r)/dt;
    b    = (r*t2 - rnew* t1)/dt;
    
    return  a * time * 1000 + b;
}
double updated_length_value(Membrane &mem, double time){
    double t1, t2, r, rnew, dt, a, b;
    // All calculations in Fs
    t1   = mem.get_Begin_update_time_in_Ps()*1000;
    t2   = mem.get_End_update_time_in_Ps()*1000;
    r    = mem.get_node_radius();
    rnew = mem.get_new_node_radius();
    dt   = t2-t1;
    // sigma_ev = a * time + b
    a    = (rnew - r)/dt;
    b    = (r*t2 - rnew* t1)/dt;
    
    return  a * time * 1000 + b;
}
/** -----------------------------------------------------------------------------
 */
void updateOpenMMforces(vector<Membrane>      &membranes,
                        MyOpenMMData*          omm,
                        double                 time,
                        MyAtomInfo             atoms[],
                        Bonds*                 bonds,
                        vector<set<int> >      &membrane_set,
                        vector<set<int> >      &actin_set,
                        vector<set<int> >      &ecm_set,
                        vector<set<int> >      &chromatin_set,
                        vector<vector<int> >   interaction_map)
{
    
    int mem_count=0;
    
//    double t1, t2, r, rnew, dt, a, b;
    
    for (int i=0; i<membranes.size(); i++) {
        //Update node radius
        if (membranes[i].get_new_node_radius() != -1) {
            
            vector<double> sigma_ev(1, updated_sigam_value(membranes[i], time));
            
            for (int node_index = mem_count; node_index < membranes[i].get_num_of_nodes() + mem_count; node_index++) {
                
                atoms[node_index].radius = sigma_ev[0];
                
                sigma_ev[0] *= OpenMM::NmPerAngstrom;
                
                for (int j=0; j<omm->EV.size(); j++) {
                    
                    omm->EV[j]->setParticleParameters(node_index, sigma_ev);
                }
            }
        }
        vector<double> params(1,0);
        omm->EV[1]->getParticleParameters(0, params);
        cout<<"param = "<<params[0]<<endl;
        
        //Update spring nominal length
//        if (membranes[i].get_new_nominal_length() != -1) {
//
//            vector<double> sigma_ev(1, updated_sigam_value(membranes[i], time));
//
//            for (int node_index = mem_count; node_index < membranes[i].get_num_of_nodes() + mem_count; node_index++) {
//
//                atoms[node_index].radius = sigma_ev[0];
//                sigma_ev[0] *= OpenMM::NmPerAngstrom;
//
//                for (int j=0; j<omm->EV.size(); j++) {
//                    omm->EV[j]->setParticleParameters(node_index, sigma_ev);
//                }
//            }
//        }
        
        mem_count += membranes[i].get_num_of_nodes();
    }
    
    cout<<"cutoff before = "<<omm->EV[0]->getCutoffDistance() <<endl;
    for (int i=0; i< omm->EV.size(); i++) {
        omm->EV[i]->updateParametersInContext(*omm->context);
    }
    cout<<"cutoff after = "<<omm->EV[0]->getCutoffDistance() <<endl;
//    exit(EXIT_SUCCESS);
    int mem_bond_count=0;
    for (int i=0; i<membranes.size(); i++) {
        if (membranes[i].get_new_node_radius() != -1) {
            for (int k=mem_bond_count; k < mem_bond_count+ membranes[i].get_num_of_node_pairs(); k++) {
                int atom1, atom2 ;
                double length, stiffness;
                omm->harmonic->getBondParameters(k, atom1, atom2, length, stiffness);
                omm->harmonic->setBondParameters(k, atom1, atom2, length*0.999, stiffness);
            }
        }
        mem_bond_count += membranes[i].get_num_of_node_pairs();
    }
    
    omm->harmonic->updateParametersInContext(*omm->context);
    
}
