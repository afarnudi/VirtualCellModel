#include "General_class_functions.h"
#include "Membrane.h"
#include "Chromatin.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "Interaction_table.hpp"

using OpenMM::Vec3;
using std::vector;
using std::set;

bool check_for_membrane_update(vector<Membrane>    &membranes,
                               double               time
                               )
{
    
    bool update_mem = false;
    for (int mem_index=0; mem_index<membranes.size(); mem_index++) {
        if (membranes[mem_index].get_update_status() &&
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
                        vector<Chromatin>                chromos,
                        MyOpenMMData*                    omm,
                        double                           time,
                        MyAtomInfo                       atoms[],
                        Bonds*                           bonds,
                        vector<set<int> >               &membrane_set,
                        NonBondInteractionMap            interaction_map)
{
    
    int mem_count=0;
    
    double t1, t2, r, rnew, dt, a, b;
    double new_node_radius=0;
    int mem_bond_count=0;
    
    for (int i=0; i<membranes.size(); i++) {
        if (membranes[i].get_update_status()) {
            if (time > membranes[i].get_Begin_update_time_in_Ps() &&
                time < membranes[i].get_End_update_time_in_Ps()         )
            {
                // All calculations in Fs
                t1   = membranes[i].get_Begin_update_time_in_Ps()*1000;
                t2   = membranes[i].get_End_update_time_in_Ps()*1000;
                r    = membranes[i].get_average_Membrane_radius();
                rnew = membranes[i].get_new_Membrane_radius();
                dt   = t2-t1;
                // sigma_ev = a * time + b
                a    = (rnew - r)/dt;
                b    = (r*t2 - rnew* t1)/dt;
                
                
                string sigma = "sigma" + generalParameters.Membrane_label + std::to_string(i) + generalParameters.Membrane_label + std::to_string(i) ;
                double mem_radius_at_current_time = (a * time * 1000 + b);
                double new_sig = 0.5*membranes[i].get_avg_node_dist()*mem_radius_at_current_time/r;
                new_node_radius = new_sig;
//                cout<<"new_node_radius "<<new_node_radius<<endl;
                if (interaction_map.get_interacton_type(i,i)!="None") {
                    omm->context->setParameter(sigma, new_sig);
                }
                
                
                for (int ch=0; ch<chromos.size(); ch++) {
                    if (interaction_map.get_interacton_type(ch+1,0)=="Excluded-Volume") {
                        sigma = "sigma" + generalParameters.Chromatin_label + std::to_string(ch) + generalParameters.Membrane_label + std::to_string(i) ;
//                        new_sig = (new_sig + chromos[ch].get_node_radius())*0.5;
                        omm->context->setParameter(sigma, new_sig);
                    } else if (interaction_map.get_interacton_type(ch+1,0)=="Lennard-Jones"){
                        for (int chind=0; chind < chromos[ch].get_num_of_node_types(); chind++) {
                            sigma = "sigma" + generalParameters.Chromatin_label + std::to_string(ch) + std::to_string(chind) + generalParameters.Membrane_label + std::to_string(i) ;
                            omm->context->setParameter(sigma, (new_sig + chromos[ch].get_node_radius())*0.5);
                        }
                        
                    }
                    
                }
                
                for (int k=mem_bond_count; k < mem_bond_count+ membranes[i].get_num_of_node_pairs(); k++) {
                    int atom1, atom2 ;
                    double length, stiffness;
                    omm->harmonic->getBondParameters(k, atom1, atom2, length, stiffness);
                    double old_length = membranes[i].get_node_pair_Nominal_Length_in_Nm(k);
                    double scaled_length = old_length*mem_radius_at_current_time/r;
                    
                    omm->harmonic->setBondParameters(k, atom1, atom2, scaled_length, stiffness);
                }
                mem_bond_count += membranes[i].get_num_of_node_pairs();
            }
        }
        
        mem_count += membranes[i].get_num_of_nodes();
    }

    omm->harmonic->updateParametersInContext(*omm->context);
    
}

void expand(vector<Chromatin>                chromos,
            MyOpenMMData*                    omm)
{
    for (int ch=0; ch<chromos.size(); ch++) {
        string epsilon = "epsilon" + generalParameters.Chromatin_label + std::to_string(ch) + generalParameters.Membrane_label + std::to_string(0) ;
        string sigma   = "sigma"   + generalParameters.Chromatin_label + std::to_string(ch) + generalParameters.Membrane_label + std::to_string(0) ;
        omm->context->setParameter(epsilon, 0);
        omm->context->setParameter(sigma, 1000);
        
    }
}
