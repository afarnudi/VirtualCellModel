//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"



void Chromatin::write_traj (string traj_name, string label){
    ofstream Trajectory;
    Trajectory.open(traj_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    for(int j=0; j< Num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
    }
}

void Chromatin::export_for_resume(char* buffer, int MD_step){
    ofstream write_resume_file;
    string resume_file_name="Results/Resume_";
    resume_file_name+=buffer;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        //Node_force=0
    }
}

void Chromatin::generate_report(char* buffer)
{
    string Report_file_name;
    Report_file_name= "Results/Membrane_Report_";
    Report_file_name+=buffer;
    Report_file_name+=".txt";
    
    ofstream Report;
    Report.open(Report_file_name.c_str());
    Report<< std:: fixed;
    Report<<"Node Mass"<< setw(20)<<Node_Mass<<endl;
//    Report<<"Radius"<< setw(20)<<Radius<<endl;
    Report<<"Minimum node pair length"<< setw(20)<<Min_node_pair_length<<endl;
    Report<<"Maximum node pair length"<< setw(20)<<Max_node_pair_length<<endl;
    Report<<"Average node pair length"<< setw(20)<<Average_node_pair_length<<endl;
    Report<<"# of Nodes "<< setw(20)<<return_num_of_nodes()<<endl;
//    Report<<"# of Triangles "<< setw(20)<<return_num_of_triangle()<<endl;
    Report<<"Spring model"<< setw(20)<<spring_model<<endl;
    Report<<"Spring coefficient"<< setw(20)<<Spring_coefficient<<endl;
//    Report<<"Bending coefficient"<< setw(20)<<Bending_coefficient<<endl;
    Report<<"Damping coefficient"<< setw(20)<<Damping_coefficient<<endl;
    if (spring_model==1)
    {
        Report<<"Membrane Spring Model:"<< setw(20)<<"FENE"<<endl;
        
    }
    if (spring_model==2)
    {
        Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;
        
    }
    
}
