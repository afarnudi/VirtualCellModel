//
//  interaction.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void Chromatin_Membrane_hard_sphere(Chromatin chromo, Membrane Mem){
    double le1,lmin;
    double deltax,deltay,deltaz,Node_distance,force;
    double interaction_strength=1000;
    int Node_A, Node_B;
    lmin=chromo.return_node_radius()+Mem.Average_node_pair_length/2.0;
    le1=lmin*2.15;
//    cout<<"in hard\n";
//    for (int i=0; i<chromo.return_num_of_nodes(); i++) {
//        for (int j=0; j<chromo.Membrane_neighbbour_node[i].size(); j++) {
//            cout<<i<<"\t"<<chromo.Membrane_neighbbour_node[i][j]<<"\n";
//        }
//
//    }
    for (int i=0; i<chromo.return_num_of_nodes(); i++) {
//        for (int j=0; j<chromo.Membrane_neighbbour_node[i].size(); j++) {
//            cout<<i<<"\t"<<chromo.Membrane_neighbbour_node[i][j]<<"\n";
//        }
        Node_B=i;
        for (int j=0; j<chromo.Membrane_neighbbour_node[i].size(); j++) {
            
            Node_A=chromo.Membrane_neighbbour_node[i][j];
            
            deltax=chromo.return_node_position(Node_B, 0)-Mem.return_node_position(Node_A, 0);
            deltay=chromo.return_node_position(Node_B, 1)-Mem.return_node_position(Node_A, 1);
            deltaz=chromo.return_node_position(Node_B, 2)-Mem.return_node_position(Node_A, 2);
            
            Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
            cout<<"lmin= "<<lmin<<"\tle1= "<<le1<<endl;
//            cout<<"Node_distance "<<Node_distance<<endl;
            force=0.0;
            double exp_le1=exp(1.0/(Node_distance-le1));
            
            if(Node_distance < le1   &&  Node_distance > lmin  )  // repulsive force
            {
                force= -(interaction_strength*exp_le1/(Node_distance-lmin))*( 1.0/(Node_distance-lmin) + 1.0/((Node_distance-le1)*(Node_distance-le1)));                 // force on i th from j
            }
            if(force<-1000.05   ||  Node_distance<lmin )
            {
                force =-1000.05-1000* ( 0.671965*Mem.return_node_radius() - Node_distance );
            }
            
            // implimentation of forces:
            Mem.add_to_force( -force*deltax/Node_distance, Node_A, 0);
            Mem.add_to_force( -force*deltay/Node_distance, Node_A, 1);
            Mem.add_to_force( -force*deltaz/Node_distance, Node_A, 2);
            
            chromo.add_to_force( force*deltax/Node_distance, Node_B, 0);
            chromo.add_to_force( force*deltay/Node_distance, Node_B, 1);
            chromo.add_to_force( force*deltaz/Node_distance, Node_B, 2);
        }
    }
}


void Chromatin_Membrane_neighbour_finder(Chromatin& chromo, Membrane Mem){
    
    chromo.Membrane_neighbbour_node.clear();
    chromo.Membrane_neighbbour_node.resize(chromo.return_num_of_nodes());
    double threshold_dist=(chromo.return_node_radius()+Mem.Average_node_pair_length/2.0)*2.15;
    for (int i=0; i<chromo.return_num_of_nodes(); i++) {
        for (int j=0; j<Mem.return_num_of_nodes(); j++) {
            
            double delta_x=chromo.return_node_position(i, 0)-Mem.return_node_position(j, 0);
            double delta_y=chromo.return_node_position(i, 1)-Mem.return_node_position(j, 1);
            double delta_z=chromo.return_node_position(i, 2)-Mem.return_node_position(j, 2);
            double dist=sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
            if (dist<threshold_dist) {
                chromo.Membrane_neighbbour_node[i].push_back(j);
//                cout<<i<<"\t"<<j<<endl;
            }
        }
        
    }
//    cout<<"in finder\n";
//    for (int i=0; i<chromo.return_num_of_nodes(); i++) {
//        for (int j=0; j<chromo.Membrane_neighbbour_node[i].size(); j++) {
//                cout<<i<<"\t"<<chromo.Membrane_neighbbour_node[i][j]<<"\n";
//        }
//
//    }
//    cout<<"\n";
}
