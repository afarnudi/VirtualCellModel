//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"

using std::string;
using std::endl;

void Chromatin::write_parameters(int MD_Step){
    string traj_file_name;
    
    traj_file_name="Results/CM_"+GenConst::trajectory_file_name+"Chromatin_"+std::to_string(chrom_index)+"_"+file_time+".txt";
    std::ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), std::ios::app);
    Trajectory << std:: fixed;
    
    for (int i=0; i<Num_of_Nodes-2; i++) {
        for (int j=i+2; Num_of_Nodes; j++) {
            Trajectory<<Contact_Matrix[i][j]<<"\t";
        }
        Trajectory<<"\n";
    }
}

void Chromatin::packing_traj (void)
{
    string energy_file_name;
    string traj_file_name;
    
    traj_file_name="Results/Relaxation/Packing_"+GenConst::trajectory_file_name+"Chromatin_"+std::to_string(chrom_index)+"_"+file_time+".xyz";
    
    std::ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), std::ios::app);
    Trajectory << std:: fixed;
    Trajectory <<Num_of_Nodes<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< Num_of_Nodes; j++) // saving trajectory
    {
        Trajectory << "chem" <<std::setprecision(5)<< std::setw(20)<<Node_Position[j][0]<<std::setw(20)<<Node_Position[j][1]<< std::setw(20)<<Node_Position[j][2]<<endl;
    }
    
}

void Chromatin::write_traj (string traj_name, string label){
    std::ofstream Trajectory;
    Trajectory.open(traj_name.c_str(), std::ios::app);
    Trajectory << std:: fixed;
    string label_A =label+"_A";
    string label_B =label+"_B";
    for(int j=0; j< Num_of_Nodes;j++) // saving trajectory
    {
        if (AB_index[j]!=0) {
            Trajectory << label_A <<std::setprecision(5)<< std::setw(20)<<Node_Position[j][0]<< std::setw(20)<<Node_Position[j][1]<< std::setw(20)<<Node_Position[j][2]<<endl;
        } else {
            Trajectory << label_B <<std::setprecision(5)<< std::setw(20)<<Node_Position[j][0]<< std::setw(20)<<Node_Position[j][1]<< std::setw(20)<<Node_Position[j][2]<<endl;
        }
        
    }
}

void Chromatin::export_for_resume(int MD_step){
    std::ofstream write_resume_file;
    string resume_file_name="Results/Resumes/Resume_Chromatin_"+std::to_string(chrom_index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        write_resume_file<<AB_index[i]<<"\n";
        //Node_force=0
    }
}

void Chromatin::generate_report(void)
{
    string Report_file_name;
    Report_file_name= "Results/Reports/Report_Chromatin_"+std::to_string(chrom_index)+"_";
    Report_file_name+=file_time;
    Report_file_name+=".txt";
    
    std::ofstream Report;
    Report.open(Report_file_name.c_str());
    Report<< std:: fixed;
    Report<<"Node Mass"<< std::setw(20)<<Node_Mass<<endl;
//    Report<<"Radius"<< setw(20)<<Radius<<endl;
    Report<<"Minimum node pair length"<< std::setw(20)<<Min_node_pair_length<<endl;
    Report<<"Maximum node pair length"<< std::setw(20)<<Max_node_pair_length<<endl;
    Report<<"Average node pair length"<< std::setw(20)<<Average_node_pair_length<<endl;
    Report<<"# of Nodes "<< std::setw(20)<<return_num_of_nodes()<<endl;
//    Report<<"# of Triangles "<< setw(20)<<return_num_of_triangle()<<endl;
    Report<<"Spring model"<< std::setw(20)<<spring_model<<endl;
    Report<<"Spring coefficient"<< std::setw(20)<<Spring_coefficient<<endl;
//    Report<<"Bending coefficient"<< setw(20)<<Bending_coefficient<<endl;
    Report<<"Damping coefficient"<< std::setw(20)<<Damping_coefficient<<endl;
    if (spring_model==1)
    {
        Report<<"Membrane Spring Model:"<< std::setw(20)<<"FENE"<<endl;
        
    }
    if (spring_model==2)
    {
        Report<<"Membrane Spring Model:"<< std::setw(20)<<"Houkian"<<endl;
        
    }
    
}
