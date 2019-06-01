//
//  interaction.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void Chromatin_Membrane_hard_sphere(Chromatin &chromo, Membrane &Mem){
    double le1,lmin;
    double deltax,deltay,deltaz,Node_distance,force;
    double interaction_strength=0.05*GenConst::MD_T*GenConst::K, temp_potential_energy=0;
    int Node_A, Node_B;
    
    lmin=Mem.Max_node_pair_length/2;
    le1=lmin*2;

    for (int i=0; i<chromo.return_num_of_nodes(); i++) {

        Node_B=i;
        
        for (int j=0; j<chromo.Membrane_neighbour_node[i].size(); j++) {
            
            Node_A=chromo.Membrane_neighbour_node[i][j];
            
            deltax=chromo.return_node_position(Node_B, 0)-Mem.return_node_position(Node_A, 0);
            deltay=chromo.return_node_position(Node_B, 1)-Mem.return_node_position(Node_A, 1);
            deltaz=chromo.return_node_position(Node_B, 2)-Mem.return_node_position(Node_A, 2);
            
            Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
            
            force=0.0;
            
            
            if(Node_distance < le1   &&  Node_distance > lmin  )  // repulsive force
            {
                double exp_le1=exp(1.0/(Node_distance-le1));
                force = ( (interaction_strength*exp_le1)/(Node_distance-lmin) )*( 1.0/(Node_distance-lmin)+1.0/( (Node_distance-le1)*(Node_distance-le1) ) );
                temp_potential_energy = interaction_strength*exp_le1/(Node_distance-lmin);
            }
//            else if(force<-1000   ||  Node_distance<lmin )
//            {
//                double c=1500;
//                force = -2.0*c*Node_distance/(lmin) + 3*c;
//                temp_potential_energy=   c*Node_distance*Node_distance/lmin -3*c*Node_distance;
//            }
            
            if (force!=0) {
                force=force/Node_distance;
                Mem.add_to_force( -force*deltax, Node_A, 0);
                Mem.add_to_force( -force*deltay, Node_A, 1);
                Mem.add_to_force( -force*deltaz, Node_A, 2);

                chromo.add_to_force( force*deltax, Node_B, 0);
                chromo.add_to_force( force*deltay, Node_B, 1);
                chromo.add_to_force( force*deltaz, Node_B, 2);
            }
            // implimentation of forces:
           
        }
    }
}


void Chromatin_Membrane_neighbour_finder(Chromatin& chromo, Membrane Mem){
    
    chromo.Membrane_neighbour_node.clear();
    chromo.Membrane_neighbour_node.resize(chromo.return_num_of_nodes());
    double threshold_dist=(chromo.return_node_radius()+Mem.Average_node_pair_length);
    for (int i=0; i<chromo.return_num_of_nodes(); i++) {
        for (int j=0; j<Mem.return_num_of_nodes(); j++) {
            
            double delta_x=chromo.return_node_position(i, 0)-Mem.return_node_position(j, 0);
            double delta_y=chromo.return_node_position(i, 1)-Mem.return_node_position(j, 1);
            double delta_z=chromo.return_node_position(i, 2)-Mem.return_node_position(j, 2);
            double dist=sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
            if (dist<threshold_dist) {
                chromo.Membrane_neighbour_node[i].push_back(j);
//                cout<<i<<"\t"<<j<<endl;
            }
        }
        
    }

}
