//
//  Membrane_spring_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"


void Membrane::excluded_volume(void){
    
    double distance, delta_x, delta_y, delta_z, force;
    
    
    for (int i=0; i<Num_of_Nodes-1; i++) {
        
        for (int j=i+1; j<Num_of_Nodes; j++) {
            force = 0;
            bool interact = true;
            for (int k=0; k<Node_neighbour_list[i].size(); k++) {
                if (Node_neighbour_list[i][k] == j) {
                    interact = false;
                }
            }
            
            if (interact) {
                
                delta_x = Node_Position[j][0] - Node_Position[i][0];
                delta_y = Node_Position[j][1] - Node_Position[i][1];
                delta_z = Node_Position[j][2] - Node_Position[i][2];
                
                double a[3]={delta_x, delta_y, delta_z};
                distance = vector_length(a)/Node_radius;
                
                if (distance < 4) {
                    double r = (distance - 1) * (distance - 1);
                    r *= r;
                    double c = 100;
                    force = c/r - c/(81);
                    
//                    cout<<"Node i :\n"
//                    <<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n"
//                    <<"Node j :\n"
//                    <<Node_Position[j][0]<<"\t"<<Node_Position[j][1]<<"\t"<<Node_Position[j][2]<<"\n";
                    
                    
                    force /= -distance;
                    Node_Force[i][0] += force*delta_x;
                    Node_Force[i][1] += force*delta_y;
                    Node_Force[i][2] += force*delta_z;
                    
                    Node_Force[j][0] += -force*delta_x;
                    Node_Force[j][1] += -force*delta_y;
                    Node_Force[j][2] += -force*delta_z;
                    
//                    cout<<"force : "<<force<<"\n"
//                    <<"Node force i :\n"
//                    <<force*delta_x<<"\t"<<force*delta_y<<"\t"<<force*delta_z<<"\n";
                }
                
            }
            
        }
    }
    
    
}
