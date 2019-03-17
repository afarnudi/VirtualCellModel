#include "Actin.h"

void Actin::Elastic_Force_Calculator(void)
{
    double delta_x, delta_y, delta_z, distance, initial_distance;
    int node_A, node_B;
    double  force=0;
    
    
    Total_Potential_Energy=0.0;
    
    for(int i=0 ;i<Num_of_Node_Pairs ; i++ ){
        
        node_A=Node_Bond_list[i][0];
        node_B=Node_Bond_list[i][1];
        
        initial_distance = Node_Bond_relaxed_length[i];
        
        delta_x = Node_Position[node_B][0]-Node_Position[node_A][0];
        delta_y = Node_Position[node_B][1]-Node_Position[node_A][1];
        delta_z = Node_Position[node_B][2]-Node_Position[node_A][2];
        
        double a[3]={delta_x, delta_y, delta_z};
        distance = vector_length(a);
        
        if (spring_model == 0) {
            force = Hookian(distance, initial_distance);
        } else if (spring_model == 1){
            force = Kelvin(distance, i);
        } else if (spring_model == 2){
            force = Maxwell(distance, i);
        }
        force/= distance;
        
        Node_Force[node_B][0] += force*delta_x;
        Node_Force[node_B][1] += force*delta_y;
        Node_Force[node_B][2] += force*delta_z;
        
        Node_Force[node_A][0] -= force*delta_x;
        Node_Force[node_A][1] -= force*delta_y;
        Node_Force[node_A][2] -= force*delta_z;
    }
}
