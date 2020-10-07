#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>



double Membrane::calc_theta_angle_ABC(int node_A, int node_B, int node_C){
    /**A is the is the middle point of the angle*/
    double AB[3], AC[3], ACxAB[3];
    AB[0] = Node_Position[node_B][0] - Node_Position[node_A][0];
    AB[1] = Node_Position[node_B][1] - Node_Position[node_A][1];
    AB[2] = Node_Position[node_B][2] - Node_Position[node_A][2];
    
    AC[0] = Node_Position[node_C][0] - Node_Position[node_A][0];
    AC[1] = Node_Position[node_C][1] - Node_Position[node_A][1];
    AC[2] = Node_Position[node_C][2] - Node_Position[node_A][2];
    
    crossvector(ACxAB, AB, AC);
    return innerproduct(AB, AC)/vector_length(ACxAB);
}



void Membrane::calculate_surface_area_with_voronoi(){
    
    int triangle1_node_A,
    triangle1_node_B,
    triangle1_node_C,
    triangle2_node_D;
    
    update_COM_position();
    vector<vector<vector<double> > > adjacent_trangle_normal_vec_list;
    adjacent_trangle_normal_vec_list.resize(Num_of_Node_Pairs);
    
    vector<vector<double> > cot_theta_list;
    cot_theta_list.resize(Num_of_Node_Pairs);
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        cot_theta_list[i].resize(2,0);
        
        triangle1_node_C = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][0];
        triangle1_node_B = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][1];
        triangle1_node_A = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][2];
        triangle2_node_D = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][3];
        
        cot_theta_list[i][0] = calc_theta_angle_ABC(triangle1_node_C, triangle1_node_A, triangle1_node_B );
        cot_theta_list[i][1] = calc_theta_angle_ABC(triangle2_node_D, triangle1_node_A, triangle1_node_B );
        
        
        adjacent_trangle_normal_vec_list[i].resize(2);
        double AB[3], AC[3], AD[3];
        
        for (int j=0; j<3; j++) {
            AB[j] = Node_Position[triangle1_node_B][j] - Node_Position[triangle1_node_A][j];
            AC[j] = Node_Position[triangle1_node_C][j] - Node_Position[triangle1_node_A][j];
            AD[j] = Node_Position[triangle2_node_D][j] - Node_Position[triangle1_node_A][j];
        }
        
        
        
        
    }
    
    surface_area_voronoi=0;
    node_voronoi_area.clear();
    node_voronoi_area.resize(Num_of_Nodes,0);
    
    for (int i=0; i<Num_of_Nodes; i++) {
//        node_voronoi_area[i]=0;
        for (int j=0; j<Node_neighbour_list[i].size(); j++) {
            //            cout<<"Node_neighbour_list["<<i<<"]["<<j<<"]="<<node_pair_vec[i][j]<<endl;
            int node_1 = i;
            int node_2 = Node_neighbour_list[node_1][j];
            int bond12 = Node_neighbour_list_respective_bond_index[node_1][j];
            
            double bond_vec[3]={Node_Position[node_1][0]-Node_Position[ node_2 ][0],
                                Node_Position[node_1][1]-Node_Position[ node_2 ][1],
                                Node_Position[node_1][2]-Node_Position[ node_2 ][2]};
            node_voronoi_area[i] += 0.125 * vector_length(bond_vec)* vector_length(bond_vec)* (cot_theta_list[bond12 ][0] + cot_theta_list[bond12 ][1]);
        }
        surface_area_voronoi += node_voronoi_area[i];
        
        
    }
}
