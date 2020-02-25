#include "Chromatin.h"

using std::cout;
using std::endl;

int Chromatin::generate_virtual_sites(void){
    vector<vector<double> > temp_Node_Position;
    temp_Node_Position = Node_Position;
    
    num_virtual_sites_per_bond = (bond_length-2*Node_radius)/(2*bond_radius);
    
    
    if (num_virtual_sites_per_bond == 0 ) {
        cout<<"bond radius too large. Checking if optimisation flag is on...\n";
        if (optimise_bond_radius) {
            cout<<"optimisation on.\n";
            bond_radius = (bond_length-2*Node_radius)/2.;
            cout<<"One virtual site with radius "<<bond_radius<<" will be used to prevent bonds from going through each other.\n";
            num_virtual_sites_per_bond =1;
        } else {
            cout<<"optimisation flag off. Virtual sites will not be used. Bonds  might cross through each other during the simulation.\n";
            return 0;
        }
    } else {
        double gaps = ( bond_length - 2*Node_radius - num_virtual_sites_per_bond*bond_radius)/(num_virtual_sites_per_bond+1);
        cout<<"bond will have "<<gaps<<" Nm gaps between them. Checking if optimisation flag is on...\n";
        if (optimise_bond_radius) {
            cout<<"optimisation on.\n";
            num_virtual_sites_per_bond++;
            bond_radius = (bond_length-2*Node_radius)/(2.*(num_virtual_sites_per_bond));
            cout<<"Virtual site with radius "<<bond_radius<<" to close the gap between virtual sites between  nodes.\n";
        } else {
            cout<<"optimisation flag off.\n"<<num_virtual_sites_per_bond<<" Virtual sites with radius of "<<bond_radius<<" Nm will be used. Smaller particles (including other virtual sites) might cross through the bonds.\n";
        }
    }
    GenConst::ChromatinVirtualSites = true;
    vector<double> virtual_sites_fill(3, 0);
    
    vector<vector<double > > temp_pos = Node_Position;
    vector<vector<double > > temp_vel = Node_Velocity;
    vector<vector<double > > temp_for = Node_Force;
    Node_Position.clear();
    Node_Velocity.clear();
    Node_Force.clear();
    Num_of_Nodes=0;
    
    Node_Position.push_back(temp_pos[0]);
    Node_Velocity.push_back(temp_vel[0]);
    Node_Force.push_back(temp_for[0]);
    Num_of_Nodes++;
    for (int i = 1; i<temp_pos.size(); i++) {
        for (int j=0; j<num_virtual_sites_per_bond; j++) {
            Node_Position.push_back(virtual_sites_fill);
            Node_Velocity.push_back(virtual_sites_fill);
            Node_Force.push_back(virtual_sites_fill);
            Num_of_Nodes++;
        }
        Node_Position.push_back(temp_pos[i]);
        Node_Velocity.push_back(temp_vel[i]);
        Node_Force.push_back(temp_for[i]);
        Num_of_Nodes++;
    }
    
    cout<<"Num_of_Nodes : "<<Num_of_Nodes<<endl;
    return 1;
}

void Chromatin::calculate_extra_virtual_bonds(void){
    vector<int> bond_pairs(2);
    int num_of_real_nodes = (Num_of_Nodes+num_virtual_sites_per_bond)/(num_virtual_sites_per_bond+1);
    cout<<"num_of_real_nodes = "<<num_of_real_nodes<<endl;
    vector<int> temp_vsite_and_bindings(3);
    vector<double> temp_vsite_weights(2);
    double vsite_distance_steps = (bond_length-2*Node_radius)/(num_virtual_sites_per_bond+1);
    
    
    if (num_of_real_nodes == 2) {
        //num_of_extra_VS_bonds = num_of_real_nodes-1 + ( (num_virtual_sites_per_bond+3)*num_virtual_sites_per_bond)/2;
        
        for (int i=1; i<num_virtual_sites_per_bond+1; i++) {
            bond_pairs[0]=0;
            bond_pairs[1]=i;
            virtual_bond_pairs.push_back(bond_pairs);
            
            temp_vsite_and_bindings[0] = i;
            temp_vsite_and_bindings[1] = 0;
            temp_vsite_and_bindings[2] = num_virtual_sites_per_bond+1;
            Vsite_and_bindings.push_back(temp_vsite_and_bindings);
            
            
            temp_vsite_weights[0] = (Node_radius + i*vsite_distance_steps )/bond_length;
            temp_vsite_weights[1] = 1 - temp_vsite_weights[0];
            Vsite_binding_weights.push_back(temp_vsite_weights);
        }
        for (int i=1; i<num_virtual_sites_per_bond+1; i++) {
            for (int j=i+1; j<num_virtual_sites_per_bond+2; j++) {
                bond_pairs[0]=i;
                bond_pairs[1]=j;
                virtual_bond_pairs.push_back(bond_pairs);
            }
        }
        
        num_of_extra_VS_bonds = int(virtual_bond_pairs.size());
        num_of_total_bonds = num_of_real_nodes - 1 + num_of_extra_VS_bonds;
        
    } else if (num_of_real_nodes >= 3){
        //num_of_extra_VS_bonds = num_of_real_nodes-1 + (num_of_real_nodes-2)*(num_virtual_sites_per_bond*(3*num_virtual_sites_per_bond+7)/2)+num_virtual_sites_per_bond*(num_virtual_sites_per_bond+3)/2;
        
        for (int realnode=0; realnode<num_of_real_nodes-2; realnode++) {
            //real node with the virtual sites on the same bond
            for (int vsite=realnode*(num_virtual_sites_per_bond+1)+1; vsite<(realnode+1)*(num_virtual_sites_per_bond+1); vsite++) {
                
                bond_pairs[0]=realnode*(num_virtual_sites_per_bond+1);
                bond_pairs[1]=vsite;
                virtual_bond_pairs.push_back(bond_pairs);
                
                temp_vsite_and_bindings[0] = vsite;
                temp_vsite_and_bindings[1] = realnode*(num_virtual_sites_per_bond+1);
                temp_vsite_and_bindings[2] = (realnode+1)*(num_virtual_sites_per_bond+1);
                Vsite_and_bindings.push_back(temp_vsite_and_bindings);
                
                
                temp_vsite_weights[0] = (Node_radius + (vsite-temp_vsite_and_bindings[1])*vsite_distance_steps )/bond_length;
                temp_vsite_weights[1] = 1 - temp_vsite_weights[0];
                Vsite_binding_weights.push_back(temp_vsite_weights);
            }
            for (int vsite1=realnode*(num_virtual_sites_per_bond+1)+1; vsite1<(realnode+1)*(num_virtual_sites_per_bond+1); vsite1++) {
                for (int vsite2=vsite1+1; vsite2<(realnode+1)*(num_virtual_sites_per_bond+1)+1; vsite2++) {
                    
                    bond_pairs[0]=vsite1;
                    bond_pairs[1]=vsite2;
                    virtual_bond_pairs.push_back(bond_pairs);
                }
            }
            for (int vsite=(realnode+1)*(num_virtual_sites_per_bond+1)+1; vsite<(realnode+2)*(num_virtual_sites_per_bond+1); vsite++) {
                
                bond_pairs[0]=realnode*(num_virtual_sites_per_bond+1);
                bond_pairs[1]=vsite;
                virtual_bond_pairs.push_back(bond_pairs);
            }
            for (int vsite1=realnode*(num_virtual_sites_per_bond+1)+1; vsite1<(realnode+1)*(num_virtual_sites_per_bond+1); vsite1++) {
                for (int vsite2=(realnode+1)*(num_virtual_sites_per_bond+1)+1; vsite2<(realnode+2)*(num_virtual_sites_per_bond+1)+1; vsite2++) {
                    
                    bond_pairs[0]=vsite1;
                    bond_pairs[1]=vsite2;
                    virtual_bond_pairs.push_back(bond_pairs);
                }
            }
        }
        for (int vsite=(num_of_real_nodes-2)*(num_virtual_sites_per_bond+1)+1; vsite<(num_of_real_nodes-1)*(num_virtual_sites_per_bond+1); vsite++) {
            
            bond_pairs[0]=(num_of_real_nodes-2)*(num_virtual_sites_per_bond+1);
            bond_pairs[1]=vsite;
            virtual_bond_pairs.push_back(bond_pairs);
            
            temp_vsite_and_bindings[0] = vsite;
            temp_vsite_and_bindings[1] = (num_of_real_nodes-2)*(num_virtual_sites_per_bond+1);
            temp_vsite_and_bindings[2] = Num_of_Nodes-1;
            Vsite_and_bindings.push_back(temp_vsite_and_bindings);
            
            
            temp_vsite_weights[0] = (Node_radius + (vsite-temp_vsite_and_bindings[1])*vsite_distance_steps )/bond_length;
            temp_vsite_weights[1] = 1 - temp_vsite_weights[0];
            Vsite_binding_weights.push_back(temp_vsite_weights);
        }
        
        for (int vsite1=(num_of_real_nodes-2)*(num_virtual_sites_per_bond+1)+1; vsite1<(num_of_real_nodes-1)*(num_virtual_sites_per_bond+1); vsite1++) {
            for (int vsite2=vsite1+1; vsite2<(num_of_real_nodes-1)*(num_virtual_sites_per_bond+1)+1; vsite2++) {
                
                bond_pairs[0]=vsite1;
                bond_pairs[1]=vsite2;
                virtual_bond_pairs.push_back(bond_pairs);
            }
        }
        
        num_of_extra_VS_bonds = int(virtual_bond_pairs.size());
        num_of_total_bonds = num_of_real_nodes - 1 + num_of_extra_VS_bonds;
    }
}
