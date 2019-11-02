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
        if (membranes[mem_index].get_num_of_nodes() != -1 &&
            time > membranes[mem_index].get_Begin_update_time_in_Ps()   ) {
            if (time < membranes[mem_index].get_End_update_time_in_Ps()) {
                update_mem = true;
            }
        }
    }
    return update_mem;
}


/** -----------------------------------------------------------------------------
 */
void updateOpenMMforces(vector<Membrane>                &membranes,
                        MyOpenMMData*                    omm,
                        double                           time,
                        MyAtomInfo                       atoms[],
                        Bonds*                           bonds,
                        vector<set<int> >               &membrane_set)
{
    
    int mem_count=0;
    
    double t1, t2, r, rnew, dt, a, b;
    
    for (int i=0; i<membranes.size(); i++) {
        if (membranes[i].get_new_node_radius() != -1) {
            if (time > membranes[i].get_Begin_update_time_in_Ps() &&
                time < membranes[i].get_End_update_time_in_Ps()         )
            {
                // All calculations in Fs
                t1   = membranes[i].get_Begin_update_time_in_Ps()*1000;
                t2   = membranes[i].get_End_update_time_in_Ps()*1000;
                r    = membranes[i].get_node_radius();
                rnew = membranes[i].get_new_node_radius();
                dt   = t2-t1;
                // sigma_ev = a * time + b
                a    = (rnew - r)/dt;
                b    = (r*t2 - rnew* t1)/dt;
                
                vector<double> sigma_ev;
                
                sigma_ev.push_back( a * time * 1000 + b);
//                for (int node_index = mem_count; node_index < membranes[i].get_num_of_nodes() + mem_count; node_index++) {
//
//                    atoms[node_index].radius = sigma_ev[0];
//                    sigma_ev[0] *= OpenMM::NmPerAngstrom;
//
//                    for (int j=0; j<omm->EV.size(); j++) {
//                        omm->EV[j]->setParticleParameters(node_index, sigma_ev);
//                    }
//                }
                
            }
        }
        
        mem_count += membranes[i].get_num_of_nodes();
    }
    
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
