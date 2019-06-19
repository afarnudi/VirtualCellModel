//
//  Interaction_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 22/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void Membrane_ECM_shared_node_force (ECM &ecm, Membrane &mem){
    
    int mem_nodes=mem.return_num_of_nodes();
    double force=0, temp_potential_energy=0;
    double delta_x=0, delta_y=0, delta_z=0, Node_distance=0;
    double sigma = 0.8;
    double epsilon = 4 * GenConst::K * GenConst::MD_T * mem.return_ECM_interaction_strength();
    
    for (int i=0; i<mem_nodes; i++) {
        if (mem.ECM_Node_neighbour_list[i].size() != 0) {
            int ecm_Node = mem.ECM_Node_neighbour_list[i][0];
            
            delta_x = mem.return_node_position(i, 0) - ecm.return_node_position(ecm_Node, 0);
            delta_y = mem.return_node_position(i, 1) - ecm.return_node_position(ecm_Node, 1);
            delta_z = mem.return_node_position(i, 2) - ecm.return_node_position(ecm_Node, 2);
            
            double a[3]={delta_x, delta_y, delta_z};
            Node_distance = vector_length(a);
            
            double r_1 = (Node_distance)/sigma;
            double r_3 = r_1*r_1*r_1;
            double r_5 = r_3*r_1*r_1;
            
            force = 2*epsilon*( 2.0/r_5 - 1.0/r_3 )/sigma;
            temp_potential_energy = epsilon * ( r_1*1.0/r_5 - r_1*1.0/r_3  );
            
            force/=Node_distance;
            
            mem.add_to_force(force*delta_x, i, 0);
            mem.add_to_force(force*delta_y, i, 1);
            mem.add_to_force(force*delta_z, i, 2);
//            cout<<"got one\n";
//            ecm.add_to_force(-force*delta_x, ecm_Node, 0);
//            ecm.add_to_force(-force*delta_y, ecm_Node, 1);
//            ecm.add_to_force(-force*delta_z, ecm_Node, 2);
        }
    }
}

//double return_ecm_membrane_node_distance(Membrane mem, int mem_node, ECM ecm, int ecm_node){
//    return sqrt( (mem.Node_Position[mem_node][0]-ecm.Node_Position[ecm_node][0])*(mem.Node_Position[mem_node][0]-ecm.Node_Position[ecm_node][0])+(mem.Node_Position[mem_node][1]-ecm.Node_Position[ecm_node][1])*(mem.Node_Position[mem_node][1]-ecm.Node_Position[ecm_node][1])+(mem.Node_Position[mem_node][2]-ecm.Node_Position[ecm_node][2])*(mem.Node_Position[mem_node][2]-ecm.Node_Position[ecm_node][2]) );
//}
//
//double return_triangle_membrane_distance(Membrane mem, int mem_node, ECM ecm, int tri_index, double tri_com[3]){
//    int node_A=ecm.Triangle_List[tri_index][0];
//    int node_B=ecm.Triangle_List[tri_index][1];
//    int node_C=ecm.Triangle_List[tri_index][2];
//    
//    tri_com[0]=(ecm.Node_Position[node_A][0]+ecm.Node_Position[node_B][0]+ecm.Node_Position[node_C][0])/3.0;
//    tri_com[1]=(ecm.Node_Position[node_A][1]+ecm.Node_Position[node_B][1]+ecm.Node_Position[node_C][1])/3.0;
//    tri_com[2]=(ecm.Node_Position[node_A][2]+ecm.Node_Position[node_B][2]+ecm.Node_Position[node_C][2])/3.0;
//    
//    return sqrt( (mem.Node_Position[mem_node][0]-tri_com[0])*(mem.Node_Position[mem_node][0]-tri_com[0])+(mem.Node_Position[mem_node][1]-tri_com[1])*(mem.Node_Position[mem_node][1]-tri_com[1])+(mem.Node_Position[mem_node][2]-tri_com[2])*(mem.Node_Position[mem_node][2]-tri_com[2]) );
//}

void particle_vesicle_shared_node_force (Membrane &particle, Membrane &vesicle){
    
    int particle_nodes=particle.return_num_of_nodes();
    double force=0, temp_potential_energy=0;
    double delta_x=0, delta_y=0, delta_z=0, Node_distance=0;
    double sigma = 0.8;
    double epsilon = 4 * GenConst::K * GenConst::MD_T * particle.return_vesicle_interaction_strength();
    
    for (int i=0; i<particle_nodes; i++) {
        if (particle.Vesicle_Node_neighbour_list[i].size() != 0) {
            int vesicle_Node = particle.Vesicle_Node_neighbour_list[i][0];
            delta_x = particle.return_node_position(i, 0) - vesicle.return_node_position(vesicle_Node, 0);
            delta_y = particle.return_node_position(i, 1) - vesicle.return_node_position(vesicle_Node, 1);
            delta_z = particle.return_node_position(i, 2) - vesicle.return_node_position(vesicle_Node, 2);
            
            double a[3]={delta_x, delta_y, delta_z};
            Node_distance = vector_length(a);
    
            double SigmaOverR= sigma/Node_distance;
            double repulsion= pow(SigmaOverR ,12);
            double attraction = pow(SigmaOverR , 6);
 
            
            force = -1*epsilon*( -12* repulsion + 6*attraction)/(Node_distance* Node_distance);
            temp_potential_energy = epsilon * (repulsion- attraction );
            
            
            particle.add_to_force(force*delta_x, i, 0);
            particle.add_to_force(force*delta_y, i, 1);
            particle.add_to_force(force*delta_z, i, 2);

            vesicle.add_to_force(-force*delta_x, vesicle_Node, 0);
            vesicle.add_to_force(-force*delta_y, vesicle_Node, 1);
            vesicle.add_to_force(-force*delta_z, vesicle_Node, 2);
        }
    }
}


void Vesicle_pointparticle_neighbour_finder (point_particle &particle, Membrane &vesicle){
    double Min_dist=1000;
    double delta_x, delta_y, delta_z, dist;
    int vesicle_nodes = vesicle.return_num_of_nodes();
    for (int i=0; i< vesicle_nodes ; i++){
        delta_x = particle.return_position(0) - vesicle.return_node_position(i, 0);
        delta_y = particle.return_position(1) - vesicle.return_node_position(i, 1);
        delta_z = particle.return_position(2) - vesicle.return_node_position(i, 2);
        double a[3]={delta_x, delta_y, delta_z};
        dist = vector_length(a);

        if (Min_dist > dist){
            Min_dist=dist;
            particle.Neighbournumber=i;
            particle.Neighbour[0]=delta_x;
            particle.Neighbour[1]=delta_y;
            particle.Neighbour[2]=delta_z;
            
            
        }
    }
    particle.Neighbour_dist=Min_dist;

        
}
void pointparticle_vesicle_shared_node_force (point_particle &particle, Membrane &vesicle){
    
    int vesicle_Node= particle.Neighbournumber;
    cout<< "vesicle Node"<<particle.Neighbournumber<<endl;
    double force=0, temp_potential_energy=0;
    double delta_x=0, delta_y=0, delta_z=0, distance=0;
    double sigma = 0.8;
    double epsilon = 4 *100* GenConst::K * GenConst::MD_T;
    double cut_off= 5;
    double SigmaOverR,repulsion,attraction;
    distance=particle.Neighbour_dist;
    if (distance < cut_off){
            
            SigmaOverR= sigma/distance;
            repulsion= pow(SigmaOverR ,12);
            attraction = pow(SigmaOverR , 6);
            delta_x= particle.Neighbour[0];
            delta_y= particle.Neighbour[1];
            delta_z= particle.Neighbour[2];
    }


            
            force = -1*epsilon*( -12* repulsion + 6*attraction)/(distance*distance);
            temp_potential_energy = epsilon * (repulsion- attraction );
            
            
            particle.add_to_force(force*delta_x ,0);
            particle.add_to_force(force*delta_y ,1);
            particle.add_to_force(force*delta_z ,2);

            vesicle.add_to_force(-force*delta_x, vesicle_Node, 0);
            vesicle.add_to_force(-force*delta_y, vesicle_Node, 1);
            vesicle.add_to_force(-force*delta_z, vesicle_Node, 2);
        
    
}

void pointparticle_pointparticle_interaction(point_particle &first,point_particle &second){
        
    double force=0, temp_potential_energy=0;
    double delta_x=0, delta_y=0, delta_z=0, distance=0;
    double sigma = first.return_P_P_sigma();
    double epsilon = first.return_P_P_epsilon()* GenConst::K * GenConst::MD_T;
    double cut_off= first.return_P_P_cut_off();
    

    double SigmaOverR,repulsion,attraction;
    delta_x= second.return_position(0)-first.return_position(0);
    delta_y= second.return_position(1)-first.return_position(1);
    delta_z= second.return_position(2)-first.return_position(2);
    double a[3]={delta_x, delta_y, delta_z};
    distance = vector_length(a);
    if (distance < cut_off){
            
            SigmaOverR= sigma/distance;
            repulsion= pow(SigmaOverR ,12);
            attraction = pow(SigmaOverR , 6);
    


            
    force = 12*epsilon*( repulsion-attraction )/(distance*distance);
    temp_potential_energy = epsilon * (repulsion- attraction );
            
            
            first.add_to_force(-force*delta_x ,0);
            first.add_to_force(-force*delta_y ,1);
            first.add_to_force(-force*delta_z ,2);
            second.add_to_force(force*delta_x ,0);
            second.add_to_force(force*delta_y ,1);
            second.add_to_force(force*delta_z ,2);
    }
        
}