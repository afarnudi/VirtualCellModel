//
//  Membrane_node_neighbour_list.cpp
//  Mem
//
//  Created by Ali Farnudi on 21/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Node_neighbour_list_constructor(){
    Node_neighbour_list.clear();
    Node_neighbour_list.resize(Num_of_Nodes);
    Node_neighbour_list_respective_bond_index.resize(Num_of_Nodes);
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        int node_1 = Node_Bond_list[i][0];
        int node_2 = Node_Bond_list[i][1];
        Node_neighbour_list[ node_1 ].push_back( node_2 );
        Node_neighbour_list[ node_2 ].push_back( node_1 );
        Node_neighbour_list_respective_bond_index[ node_1 ].push_back(i);
        Node_neighbour_list_respective_bond_index[ node_2 ].push_back(i);
    }
}

void Membrane::Bond_triangle_neighbour_list_constructor(){
    Bond_triangle_neighbour_indices.clear();
    Bond_triangle_neighbour_indices.resize(Num_of_Node_Pairs, -1);
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
    
        int node_1= Node_Bond_list[i][0];
        int node_2= Node_Bond_list[i][1];
        
        for (int j=0; j<Num_of_Triangle_Pairs; j++) {
            
            
            int triangle1_node_B = Triangle_Pair_Nodes[ j ][1];
            int triangle1_node_A = Triangle_Pair_Nodes[ j ][2];
            
            if ( (triangle1_node_A == node_1 && triangle1_node_B == node_2) ||
                 (triangle1_node_A == node_2 && triangle1_node_B == node_1) ) {
                Bond_triangle_neighbour_indices[i]=j;
                break;
            }
        }
        if (Bond_triangle_neighbour_indices[i]==-1) {
            string errorMessage = TWARN;
            errorMessage+="There is an error in the Membrane's 'Bond_triangle_neighbour_list_constructor'. it would seem that there is a node_pair that cannot be associated with two triangles. If the mesh is an open surface, this could be due to the bonds lying on the perimeter of the surface. This means that the voronoi area calculations will breakdown.\n";
            errorMessage+= TRESET;
            cout<<errorMessage<<endl;
            open_surface = true;
            WantGeometricProps = false;
//            throw std::runtime_error(errorMessage);
        }
    }
}
