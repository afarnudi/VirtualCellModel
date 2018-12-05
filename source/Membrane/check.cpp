//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::check(void){
    Min_node_pair_length=1000;
    Max_node_pair_length=0;
    Average_node_pair_length=0;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        double dist=0;
        dist=sqrt((Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])*(Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])+(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])*(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])+(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2])*(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2]));
        Average_node_pair_length+=dist;
        
        if (Min_node_pair_length>dist) {
            Min_node_pair_length=dist;
        }
        if (dist>Max_node_pair_length) {
            Max_node_pair_length=dist;
        }
    }
    Average_node_pair_length/=Num_of_Node_Pairs; cout<<"Max node distance="<<Max_node_pair_length<<"\tmin node distance="<<Min_node_pair_length<<"\tAverage node distance="<<Average_node_pair_length<<endl;
//    cout<<"=============================\n";
//    cout<<"spring_model= "<<spring_model<<endl;
//    cout<<"=============================\n";
    if ((Min_node_pair_length*2<Max_node_pair_length) && Bending_coefficient!=0) {
        cout<<"Initial node distances are not ready/optimised for triangle bending calculations. A few MD steps will be added to the beginning of the simulation to avoid programme break down.\n";
        node_distance_correction();
        calculate_mesh_properties();
//        exit(EXIT_FAILURE);
    }
}

void Membrane::node_distance_correction(void){
    double MD_relax_Steps=1000;
    double slope=(1.88*Min_node_pair_length-Max_node_pair_length)/MD_relax_Steps, max=Max_node_pair_length;
//    cout<<"spring coefficient= "<<Spring_coefficient<<endl;
    double temp_Damping_coefficient=Damping_coefficient;
    Damping_coefficient=0;
    for(int MD_Step=0 ;MD_Step<=MD_relax_Steps ; MD_Step++){
        //Setting the min angle of triangles to 20 dgrees or pi/9
        
        Max_node_pair_length=slope*MD_Step+max;
//        cout<<"Max_node_pair_length= "<<Max_node_pair_length<<endl;
        for (int i=0; i<100; i++) {
            MD_Evolution_beginning(GenConst::MD_Time_Step);
            
//            FENE();
            potential_1();
//            Bending_potetial();
            
            MD_Evolution_end(GenConst::MD_Time_Step);
            
            
        }
        Results("mem");
//        Thermostat_2(GenConst::MD_KT);
//        calculate_mesh_properties();
        
        
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
//    cout<<"Max_node_pair_length= "<<Max_node_pair_length<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Velocity[i][j]=0;
            Node_Force[i][j]=0;
        }
    }
    Damping_coefficient=temp_Damping_coefficient;
}

void Membrane::calculate_mesh_properties(void){
    double temp_min=1000;
    double temp_max=0;
    double temp_avg=0;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        double dist=0;
        dist=sqrt((Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])*(Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])+(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])*(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])+(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2])*(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2]));
        temp_avg+=dist;
        
        if (temp_min>dist) {
            temp_min=dist;
        }
        if (dist>temp_max) {
            temp_max=dist;
        }
    }
    temp_avg/=Num_of_Node_Pairs;
    cout<<"Max node distance="<<temp_max<<"\tmin node distance="<<temp_min<<"\tAverage node distance="<<temp_avg<<endl;
}


void Membrane::Results (string label)
{
    //output file names:
    
    string energy_file_name;
    string traj_file_name;
    
    traj_file_name="Results/PreResults";
    traj_file_name+=".xyz";
    //trajectory:
    
    ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    Trajectory <<Num_of_Nodes<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< Num_of_Nodes;j++) // saving trajectory
    {
//        if (j==5 || j==220) {
//            Trajectory << "b" <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
//        } else {
            Trajectory << label <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
//        }
    }
    
    
    
    
    //    //Energy:
    //    ofstream Potential_Energy;
    //    Potential_Energy.open(energy_file_name.c_str(), ios::app);
    //    Potential_Energy<<membrane.Total_Potential_Energy<<endl;
}
