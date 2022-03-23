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
    vector<vector<double> > bond_normal_vec_list;
    bond_normal_vec_list.resize(Num_of_Node_Pairs);
    
    cot_theta_list.clear();
    cot_theta_list.resize(Num_of_Node_Pairs);
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        cot_theta_list[i].resize(2,0);
        
        triangle1_node_C = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][0];
        triangle1_node_B = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][1];
        triangle1_node_A = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][2];
        triangle2_node_D = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][3];
        
        cot_theta_list[i][0] = calc_theta_angle_ABC(triangle1_node_C, triangle1_node_A, triangle1_node_B );
        cot_theta_list[i][1] = calc_theta_angle_ABC(triangle2_node_D, triangle1_node_A, triangle1_node_B );
        
        
        bond_normal_vec_list[i].resize(3,0);
        double AB[3], AC[3], AD[3], outwardvector[3];
        double ABxAC[3], ABxAD[3];
        
        for (int j=0; j<3; j++) {
            AB[j] = Node_Position[triangle1_node_B][j] - Node_Position[triangle1_node_A][j];
            AC[j] = Node_Position[triangle1_node_C][j] - Node_Position[triangle1_node_A][j];
            AD[j] = Node_Position[triangle2_node_D][j] - Node_Position[triangle1_node_A][j];
            outwardvector[j] = Node_Position[triangle1_node_A][j] - COM_position[j];
        }
        
        crossvector(ABxAC, AB, AC);
        crossvector(ABxAD, AB, AD);
        
        if (innerproduct(outwardvector, ABxAC)<0) {
            ABxAC[0]=-ABxAC[0];
            ABxAC[1]=-ABxAC[1];
            ABxAC[2]=-ABxAC[2];
        }
        if (innerproduct(outwardvector, ABxAD)<0) {
            ABxAD[0]=-ABxAD[0];
            ABxAD[1]=-ABxAD[1];
            ABxAD[2]=-ABxAD[2];
        }
        double veclength=0;
        for (int j=0; j<3; j++) {
            bond_normal_vec_list[i][j] = ABxAC[j]+ABxAD[j];
            veclength += bond_normal_vec_list[i][j]*bond_normal_vec_list[i][j];
        }
        veclength  =sqrt(veclength);
        for (int j=0; j<3; j++) {
            bond_normal_vec_list[i][j] /= veclength;
        }
        
        
        
    }
    
    surface_area_voronoi=0;
    node_voronoi_area.clear();
    node_voronoi_area.resize(Num_of_Nodes,0);
    
    node_voronoi_normal_vec.clear();
    node_voronoi_normal_vec.resize(Num_of_Nodes);
    
    for (int i=0; i<Num_of_Nodes; i++) {
        node_voronoi_normal_vec[i].resize(3,0);
        for (int j=0; j<Node_neighbour_list[i].size(); j++) {
            
            int node_1 = i;
            int node_2 = Node_neighbour_list[node_1][j];
            int bond12 = Node_neighbour_list_respective_bond_index[node_1][j];
            
            double bond_vec[3]={Node_Position[node_1][0]-Node_Position[ node_2 ][0],
                                Node_Position[node_1][1]-Node_Position[ node_2 ][1],
                                Node_Position[node_1][2]-Node_Position[ node_2 ][2]};
            node_voronoi_area[i] += 0.125 * vector_length(bond_vec)* vector_length(bond_vec)* (cot_theta_list[bond12 ][0] + cot_theta_list[bond12 ][1]);
            
            
            for (int k=0; k<3; k++) {
                node_voronoi_normal_vec[i][k] += bond_normal_vec_list[bond12][k];
            }
            
        }
        
        surface_area_voronoi += node_voronoi_area[i];
        double nodevec[3]={node_voronoi_normal_vec[i][0],
                           node_voronoi_normal_vec[i][1],
                           node_voronoi_normal_vec[i][2]};
        
        double veclength = vector_length(nodevec);
        for (auto &coord: node_voronoi_normal_vec[i]){
            coord/=veclength;
        }
        
    }
}
