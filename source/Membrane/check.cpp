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
    Average_node_pair_length/=Num_of_Node_Pairs;
    
    cout<<"Max node distance="<<Max_node_pair_length<<"\tmin node distance="<<Min_node_pair_length<<"\tAverage node distance="<<Average_node_pair_length<<endl;

    if ((Min_node_pair_length*2<Max_node_pair_length) && Bending_coefficient!=0) {
        if (!Relaxation) {
            cout<<"\nWarning:Initial node distances are not ready/optimised for triangle bending calculations. If the Membrane is not attached to an Actin or another network we strongly recommend turning on the 'Relaxation' flag or switching off the bending.\n\n";
        } else {
            cout<<"\nInitial node distances are not ready/optimised for triangle bending calculations. A few MD steps will be added to the beginning of the simulation to avoid programme break down.\n\n";
            Relaxation_Process_Model=2;
        }
    }
}
/*
void Membrane::node_distance_correction(void){
    cout<<"\nBeginnig the Relaxation\nProgress:\n";
     int progress=0;
    double MD_relax_Steps_1=2000;//2000;
    double MD_relax_Steps_2=2500;//20000;
    double MD_relax_Steps_3=500;
    string pov_relaxation_file_name ="Results/Relaxation/Relaxation_POV_"+GenConst::trajectory_file_name+"Membrane_"+to_string(index)+"_"+file_time;
    
    double slope=(Max_node_pair_length/1.8-Min_node_pair_length)/MD_relax_Steps_1, min=Min_node_pair_length;
//    cout<<"spring coefficient= "<<Spring_coefficient<<endl;
    double temp_Damping_coefficient=Damping_coefficient;
    double temp_Bending_coefficient=Bending_coefficient;
    Damping_coefficient=0.5;
    Bending_coefficient=70*GenConst::MD_T*GenConst::K;
    for(int MD_Step=0 ;MD_Step<=MD_relax_Steps_1 ; MD_Step++){
        //Setting the min angle of triangles to 20 dgrees or pi/9
        Min_node_pair_length=slope*MD_Step+min;

        for (int i=0; i<100; i++) {
            MD_Evolution_beginning(GenConst::MD_Time_Step);
            Relaxation_potential();
            Bending_potetial_2(0);
            MD_Evolution_end(GenConst::MD_Time_Step);
        }
        if (MD_Step!=0 && MD_Step%10==0) {
            relaxation_traj();
            Damper_check(MD_Step);}
        if (int(100*MD_Step/(MD_relax_Steps_1+ MD_relax_Steps_2))>progress){
            cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step*100<<"\r" << std::flush;
        progress+=5;}
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
//    check();
//    Damping_coefficient=2;
//    Bending_coefficient=35*GenConst::MD_T*GenConst::K;
//=======
//    }
    Damping_coefficient=1;
//>>>>>>> Thermo
    for(int MD_Step=0 ;MD_Step<=MD_relax_Steps_2 ; MD_Step++){
        //Setting the min angle of triangles to 20 dgrees or pi/9
//        Max_node_pair_length=slope*MD_Step+max;
        //        cout<<"Max_node_pair_length= "<<Max_node_pair_length<<endl;
        for (int i=0; i<100; i++) {
            MD_Evolution_beginning(GenConst::MD_Time_Step);
            Relaxation_potential();
            Bending_potetial_2(0);
            MD_Evolution_end(GenConst::MD_Time_Step);
        }
        if(MD_Step >= MD_relax_Steps_3)
        {Damping_coefficient=0.9*Damping_coefficient;}
        if (MD_Step!=0 && MD_Step%10==0) {
            relaxation_traj();
           // Damper_check(MD_Step+MD_relax_Steps_1);
        }
        if ( int( 100*(MD_Step+MD_relax_Steps_1) /(MD_relax_Steps_1+ MD_relax_Steps_2))>progress){
            cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step*100+ MD_relax_Steps_1*100<<"\r" << std::flush;
        progress+=5;}
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)

    double min_node_distance=10000000000000;
    double temp_dist;
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Velocity[i][j]=0;
            Node_Force[i][j]=0;
        }
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        temp_dist=vector_length(a);
        if (temp_dist < min_node_distance) {
            min_node_distance=temp_dist;
        }
    }
    min_radius_after_relaxation=min_node_distance;
    
    Damping_coefficient=temp_Damping_coefficient;
    Bending_coefficient=temp_Bending_coefficient;
    export_relaxed(0);
    write_pov_traj(pov_relaxation_file_name, to_string(index),  MD_relax_Steps_1+ MD_relax_Steps_2-1);
}*/

void Membrane::calculate_mesh_properties(void){
    double temp_min=1000;
    double temp_max=0;
    double temp_avg=0;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        double dist=0;
        dist=sqrt((Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])*(Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0]) + (Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])*(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1]) + (Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2])*(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2]));
        temp_avg+=dist;
        
        if (temp_min>dist) {
            temp_min=dist;
        }
        if (dist>temp_max) {
            temp_max=dist;
        }
    }
    temp_avg/=Num_of_Node_Pairs;
//    Max_node_pair_length=temp_max;
//    Min_node_pair_length=temp_min;
    Average_node_pair_length=temp_avg;
    cout<<"Distance properties after relaxation:\n";
    cout<<"Max node distance="<<temp_max<<"\tMin node distance="<<temp_min<<endl;
    cout<<"parameters set on:\n";
    cout<<"Max_node_pair_length="<<Max_node_pair_length<<"\tMin_node_pair_length="<<Min_node_pair_length<<endl;
}

double Membrane::Average_velocity_squared(){
    double average_velocity_squard=0;
    for (int i=0; i< Num_of_Nodes; i++){
        average_velocity_squard+=sqrt(Node_Velocity[i][0]*Node_Velocity[i][0] +Node_Velocity[i][1]*Node_Velocity[i][1] + Node_Velocity[i][2]*Node_Velocity[i][2] );
    }
    average_velocity_squard= average_velocity_squard/Num_of_Nodes;
    return(average_velocity_squard);
}

