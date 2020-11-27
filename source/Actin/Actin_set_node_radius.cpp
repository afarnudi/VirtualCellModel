//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 24/09/2020.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"

using std::cout;
using std::endl;

void Actin::set_node_radius(void){
    Node_radius.clear();
    
    
    if (Node_radius_stat=="Av") {
        Node_radius.resize(Num_of_Nodes, 0.5*Average_node_pair_length);
        cout<<"Node radii set to half of bond average distances."<<endl;
    } else if (Node_radius_stat=="Au"){
        Node_radius.resize(Num_of_Nodes, 0);
        for (int i=0; i<Num_of_Nodes; i++) {
            double avgdist=0;
            for (int j=0; j<Node_neighbour_list[i].size(); j++) {
                //            cout<<"Node_neighbour_list["<<i<<"]["<<j<<"]="<<node_pair_vec[i][j]<<endl;
                int node_1 = i;
                int node_2 = Node_neighbour_list[node_1][j];
                int bond12 = Node_neighbour_list_respective_bond_index[node_1][j];
                
                double bond_vec[3]={Node_Position[node_1][0]-Node_Position[ node_2 ][0],
                                    Node_Position[node_1][1]-Node_Position[ node_2 ][1],
                                    Node_Position[node_1][2]-Node_Position[ node_2 ][2]};
                avgdist += vector_length(bond_vec);
            }
            Node_radius[i] = avgdist/Node_neighbour_list[i].size();
            
            
        }
    } else {
        Node_radius.resize(Num_of_Nodes, stod(Node_radius_stat));
    }
    
}
