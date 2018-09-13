//
//  interaction.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void interaction_1(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag){
    //Build neighbour list
    
    if (MD_Step% 1000 ==0) {
        if (!costume_interaction_flag) {
            update_neighbour_list(membrane, ecm, ECM_membrane_neighbour_list, costume_interaction_flag);
        } else {
            update_neighbour_list_2(membrane, ecm, ECM_membrane_neighbour_list);
        }
        
    }
    
    //Force Implimintation
    
//    double average_force_mag=0;
//    double num=0;
    for (int i=0; i<ecm.return_num_of_nodes(); i++) {
//        cout<<ECM_membrane_neighbour_list[i]<<endl;
        if (ECM_membrane_neighbour_list[i]!=-1) {
//            cout<<"i'm in\n";
            int mem_index=ECM_membrane_neighbour_list[i];
            double ecm_mem_node_distance=return_ecm_membrane_node_distance(membrane, mem_index, ecm, i);
            double delta_x=membrane.Node_Position[mem_index][0]-ecm.Node_Position[i][0];
            double delta_y=membrane.Node_Position[mem_index][1]-ecm.Node_Position[i][1];
            double delta_z=membrane.Node_Position[mem_index][2]-ecm.Node_Position[i][2];
            double reduced_radius=1.5*ecm.return_epsilon()/ecm_mem_node_distance;
            double reduced_radius_squared=reduced_radius*reduced_radius;
            double reduced_radius_cubed=reduced_radius_squared*reduced_radius;
            double reduced_radius_quint=reduced_radius_cubed*reduced_radius_squared;
            double force_magnitude=4.0*ecm.return_sigma()*(reduced_radius_quint- reduced_radius_cubed)/(1.5*ecm.return_epsilon()*ecm_mem_node_distance);
//            cout<<force_magnitude<<endl;
//            average_force_mag+=force_magnitude*delta_y;
//            num++;
//            ecm.ECM_Node_Force[i][0]+=force_magnitude*delta_x;
//            ecm.ECM_Node_Force[i][1]+=force_magnitude*delta_y;
//            ecm.ECM_Node_Force[i][2]+=force_magnitude*delta_z;
//            
            membrane.Membrane_Node_Force[mem_index][0] += -force_magnitude*delta_x;
            membrane.Membrane_Node_Force[mem_index][1] += -force_magnitude*delta_y;
            membrane.Membrane_Node_Force[mem_index][2] += -force_magnitude*delta_z;
//            cout<<membrane.Membrane_Node_Force[mem_index][0]<<"\t"<<            membrane.Membrane_Node_Force[mem_index][1]<<"\t"<<            membrane.Membrane_Node_Force[mem_index][2]<<endl;
            
        }
    }
//    cout<<"average_force_mag y direction="<<average_force_mag/num<<endl;
}

void interaction_2(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag){
    //Build neighbour list
    
    if (MD_Step% 1000 ==0) {
        if (!costume_interaction_flag) {
            triangle_update_neighbour_list(membrane, ecm, ECM_membrane_neighbour_list, costume_interaction_flag);
        } else {
            triangle_update_neighbour_list_2(membrane, ecm, ECM_membrane_neighbour_list);
        }
        
    }
    
    //Force Implimintation
    
    //    double average_force_mag=0;
    //    double num=0;
    for (int i=0; i<ecm.return_num_of_triangles(); i++) {
        //        cout<<ECM_membrane_neighbour_list[i]<<endl;
        if (ECM_membrane_neighbour_list[i]!=-1) {
            //            cout<<"i'm in\n";
            int mem_index=ECM_membrane_neighbour_list[i];
            
            double tri_com[3];
            double tri_mem_node_distance=return_triangle_membrane_distance(membrane, mem_index, ecm, i, tri_com);
            bool barrier=barrier_2(membrane, mem_index);
            
            if (barrier) {
                membrane.Membrane_Node_Velocity[mem_index][0]*=-1;
                membrane.Membrane_Node_Velocity[mem_index][1]*=-1;
                membrane.Membrane_Node_Velocity[mem_index][2]*=-1;
            }
            
            double delta_x=membrane.Node_Position[mem_index][0]-tri_com[0];
            double delta_y=membrane.Node_Position[mem_index][1]-tri_com[1];
            double delta_z=membrane.Node_Position[mem_index][2]-tri_com[2];
            
            double reduced_radius=1.5*ecm.return_epsilon()/tri_mem_node_distance;
            double reduced_radius_squared=reduced_radius*reduced_radius;
            double reduced_radius_cubed=reduced_radius_squared*reduced_radius;
            double reduced_radius_quint=reduced_radius_cubed*reduced_radius_squared;
            
            double force_magnitude=4.0*ecm.return_sigma()*(reduced_radius_quint- reduced_radius_cubed)/(1.5*ecm.return_epsilon()*tri_mem_node_distance);
            //            cout<<force_magnitude<<endl;
            //            average_force_mag+=force_magnitude*delta_y;
            //            num++;
            //            ecm.ECM_Node_Force[i][0]+=force_magnitude*delta_x;
            //            ecm.ECM_Node_Force[i][1]+=force_magnitude*delta_y;
            //            ecm.ECM_Node_Force[i][2]+=force_magnitude*delta_z;
            //
            membrane.Membrane_Node_Force[mem_index][0] += -force_magnitude*delta_x;
            membrane.Membrane_Node_Force[mem_index][1] += -force_magnitude*delta_y;
            membrane.Membrane_Node_Force[mem_index][2] += -force_magnitude*delta_z;
            //            cout<<membrane.Membrane_Node_Force[mem_index][0]<<"\t"<<            membrane.Membrane_Node_Force[mem_index][1]<<"\t"<<            membrane.Membrane_Node_Force[mem_index][2]<<endl;
            
        }
    }
    //    cout<<"average_force_mag y direction="<<average_force_mag/num<<endl;
}

void interaction_3(Membrane &membrane){
    for (int i=0; i<membrane.return_num_of_nodes(); i++) {
        double y=membrane.Node_Position[i][1];
        if (y<2){
            double epsilon=1;
            double sigma=10;
            double y_reduced=1.5*epsilon/y;
            double y_reduced_3=y_reduced*y_reduced*y_reduced, y_reduced_5=y_reduced_3*y_reduced*y_reduced;
            double force=4*sigma*(y_reduced_5-y_reduced_3)/(1.5*epsilon);
            membrane.Membrane_Node_Force[i][1] += -force;
            
        }
    }
}

void interaction_4(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag){
    //Build neighbour list
    
    if (MD_Step% 1000 ==0) {
        if (!costume_interaction_flag) {
            triangle_update_neighbour_list_new(membrane, ecm, ECM_membrane_neighbour_list, costume_interaction_flag);
        } else {
            triangle_update_neighbour_list_new_2(membrane, ecm, ECM_membrane_neighbour_list);
        }
        
    }
    
    //Force Implimintation
    
    //    double average_force_mag=0;
    //    double num=0;
    for (int i=0; i<ecm.return_num_of_triangles(); i++) {
        //        cout<<ECM_membrane_neighbour_list[i]<<endl;
        if (ECM_membrane_neighbour_list[i]!=-1) {
            //            cout<<"i'm in\n";
            int mem_index=ECM_membrane_neighbour_list[i];
            
            double tri_com[3];
            double tri_mem_node_distance=return_triangle_membrane_distance(membrane, mem_index, ecm, i, tri_com);
            bool barrier=barrier_2(membrane, mem_index);
            
            if (barrier) {
                membrane.Membrane_Node_Velocity[mem_index][0]*=-1;
                membrane.Membrane_Node_Velocity[mem_index][1]*=-1;
                membrane.Membrane_Node_Velocity[mem_index][2]*=-1;
            }
            
            double delta_x=membrane.Node_Position[mem_index][0]-tri_com[0];
            double delta_y=membrane.Node_Position[mem_index][1]-tri_com[1];
            double delta_z=membrane.Node_Position[mem_index][2]-tri_com[2];
            
            double reduced_radius=1.5*ecm.return_epsilon()/tri_mem_node_distance;
            double reduced_radius_squared=reduced_radius*reduced_radius;
            double reduced_radius_cubed=reduced_radius_squared*reduced_radius;
            double reduced_radius_quint=reduced_radius_cubed*reduced_radius_squared;
            
            double force_magnitude=4.0*ecm.return_sigma()*(reduced_radius_quint- reduced_radius_cubed)/(1.5*ecm.return_epsilon()*tri_mem_node_distance);
            //            cout<<force_magnitude<<endl;
            //            average_force_mag+=force_magnitude*delta_y;
            //            num++;
            //            ecm.ECM_Node_Force[i][0]+=force_magnitude*delta_x;
            //            ecm.ECM_Node_Force[i][1]+=force_magnitude*delta_y;
            //            ecm.ECM_Node_Force[i][2]+=force_magnitude*delta_z;
            //
            membrane.Membrane_Node_Force[mem_index][0] += -force_magnitude*delta_x;
            membrane.Membrane_Node_Force[mem_index][1] += -force_magnitude*delta_y;
            membrane.Membrane_Node_Force[mem_index][2] += -force_magnitude*delta_z;
            //            cout<<membrane.Membrane_Node_Force[mem_index][0]<<"\t"<<            membrane.Membrane_Node_Force[mem_index][1]<<"\t"<<            membrane.Membrane_Node_Force[mem_index][2]<<endl;
            
        }
    }
    //    cout<<"average_force_mag y direction="<<average_force_mag/num<<endl;
}


void Node_ecm_Barrier(Membrane &membrane, ECM &ecm, vector<int> ECM_membrane_neighbour_list)
{
    double triangle_COM_position[3]; // coordinates to the centre of mass of the triangles
    double node_ecm_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3], AB[3], AC[3]; // normal vector of membrane
    double node_triangle_distance_vector[3];
    double relevant_velocity[3];
    
    vector<int> membrane_list;
    for (int i=0; i<ecm.return_num_of_nodes(); i++) {
        if (ECM_membrane_neighbour_list[i]!=-1) {
            membrane_list.push_back(ECM_membrane_neighbour_list[i]);
        }
    }
    for (int i=0; i<membrane_list.size()-1; i++) {
        int temp_node=membrane_list[i];
        for (int j=i+1; j<membrane_list.size(); j++) {
            if (membrane_list[j]==temp_node) {
                membrane_list.erase(membrane_list.begin()+j);
                j--;
            }
        }
    }
    
    
    for (int i =0; i < ecm.return_num_of_triangles(); i++)
    {
        int node_A=ecm.ECM_triangle_list[i][0];
        int node_B=ecm.ECM_triangle_list[i][1];
        int node_C=ecm.ECM_triangle_list[i][2];
        triangle_COM_position[0]=(ecm.Node_Position[node_A][0] + ecm.Node_Position[node_B][0] + ecm.Node_Position[node_C][0])/3.0;
        triangle_COM_position[1]=(ecm.Node_Position[node_A][1] + ecm.Node_Position[node_B][1] + ecm.Node_Position[node_C][1])/3.0;
        triangle_COM_position[2]=(ecm.Node_Position[node_A][2] + ecm.Node_Position[node_B][2] + ecm.Node_Position[node_C][2])/3.0;
        
        for (int index=0; index<membrane_list.size(); index++)
        {
            int temp_node_index=membrane_list[index];
            node_ecm_distance_amplitude = sqrt( (triangle_COM_position[0]-membrane.Node_Position[temp_node_index][0])*(triangle_COM_position[0]- membrane.Node_Position[temp_node_index][0]) + (triangle_COM_position[1]- membrane.Node_Position[temp_node_index][1])*(triangle_COM_position[1]- membrane.Node_Position[temp_node_index][1]) + (triangle_COM_position[2]- membrane.Node_Position[temp_node_index][2])*(triangle_COM_position[2]- membrane.Node_Position[temp_node_index][2]) );
            if (  node_ecm_distance_amplitude < 0.9*ecm.return_epsilon()  )
            {
                AB[0]=ecm.Node_Position[ node_B][0]-ecm.Node_Position[node_A][0];
                AB[1]=ecm.Node_Position[ node_B][1]-ecm.Node_Position[ node_A][1];
                AB[2]=ecm.Node_Position[ node_B][2]-ecm.Node_Position[ node_A][2];
                AC[0]=ecm.Node_Position[ node_C][0]-ecm.Node_Position[ node_A][0];
                AC[1]=ecm.Node_Position[ node_C][1]-ecm.Node_Position[ node_A][1];
                AC[2]=ecm.Node_Position[ node_C][2]-ecm.Node_Position[ node_A][2];
                crossvector(ABxAC,AB,AC);
                
                node_triangle_distance_vector[0]=membrane.Node_Position[temp_node_index][0]-triangle_COM_position[0];
                node_triangle_distance_vector[1]=membrane.Node_Position[temp_node_index][1]-triangle_COM_position[1];
                node_triangle_distance_vector[2]=membrane.Node_Position[temp_node_index][2]-triangle_COM_position[2];
                
                relevant_velocity[0] = membrane.Membrane_Node_Velocity[temp_node_index][0];
                relevant_velocity[1] = membrane.Membrane_Node_Velocity[temp_node_index][1];
                relevant_velocity[2] = membrane.Membrane_Node_Velocity[temp_node_index][2];
                
                if    ( innerproduct(relevant_velocity,ABxAC)<0)
                {
                    membrane.Membrane_Node_Velocity[temp_node_index][0]*=-1;
                    membrane.Membrane_Node_Velocity[temp_node_index][1]*=-1;
                    membrane.Membrane_Node_Velocity[temp_node_index][2]*=-1;
                }//END OF:  if    ( (abs( perpendicular_distance )<Actin_Membrane_Radius_of_Hard_Sphere_Interaction) &&
            }//END OF: if (  actin_membrane_distance_amplitude < sqrt(0.43*a * 0.43*a +
            
        }//END OF: for (int actin_counter=Actin_Membrane_shared_num_of_Nodes;actin_counter<Actin_num_of_Nodes;
        
    }//END OF: for (int i =0; i < Membrane_num_of_Triangles  ; i++)
    
}
