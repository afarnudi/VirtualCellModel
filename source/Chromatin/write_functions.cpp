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
using std::cout;

void Chromatin::write_parameters(int MD_Step){
    string traj_file_name;
    
    traj_file_name="Results/CM_"+GenConst::trajectory_file_name+"Chromatin_"+std::to_string(index)+"_"+file_time+".txt";
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
    
    traj_file_name="Results/Relaxation/Packing_"+GenConst::trajectory_file_name+"Chromatin_"+std::to_string(index)+"_"+file_time+".xyz";
    
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
    string resume_file_name="Results/Resumes/Resume_Chromatin_"+std::to_string(index)+"_";
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
    std::string Report_file_name;
    Report_file_name= "Results/Reports/Report_Chromatin_" + std::to_string(index) + "_";
    Report_file_name+=file_time;
    Report_file_name+=".txt";
    
    std::ofstream Report;
    Report.open(Report_file_name.c_str());
    
    if (!Report.is_open()){
        cout<<"Couldn't generate report. Please check that the required directories exist. The required directory should be ./bin/Results/Reports or ./Results/Reports"<<endl;
        Report.close();
        cout<<"Generating report in an alternative directory:\nbin/Results/Reports/"<<endl;
        Report_file_name= "bin/Results/Reports/Report_Chromatin_" + std::to_string(index) + "_";
        Report_file_name+=file_time;
        Report_file_name+=".txt";
        Report.open(Report_file_name.c_str());
    }
    
    Report<< std:: fixed;
    Report<<"MD Params in the general configuration file:\n---------------\n";
    Report<<"Simulation_Time_In_Ps\t"<<GenConst::Simulation_Time_In_Ps<<endl;
    Report<<"Step_Size_In_Fs\t"<<GenConst::Step_Size_In_Fs<<endl;
    Report<<"Report_Interval_In_Fs\t"<<GenConst::Report_Interval_In_Fs<<endl;
    Report<<"MD_num_of_steps\t"<<GenConst::MD_num_of_steps<<endl;
    Report<<"MD_traj_save_step\t"<<GenConst::MD_traj_save_step<<endl;
    Report<<"K\t"<<GenConst::K<<endl;
    Report<<"MD_T\t"<<GenConst::MD_T<<endl;
    Report<<"MD_thrmo_step\t"<<GenConst::MD_thrmo_step<<endl;
    Report<<"Bussi_tau\t"<<GenConst::Bussi_tau<<endl;
    Report<<"MC_step\t"<<GenConst::MC_step<<endl;
    Report<<"Mem_fluidity\t"<<GenConst::Mem_fluidity<<endl;
    Report<<"Lbox\t"<<GenConst::Lbox<<endl;
    Report<<"Periodic_condtion_status\t"<<GenConst::Periodic_condtion_status<<endl;
    Report<<"trajectory_file_name\t"<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"Chromatin Params:\n---------------\n";
    Report<<"Num_of_Nodes\t"<<Num_of_Nodes<<endl;
    Report<<"Node_Mass\t"<<Node_Mass<<endl;
    Report<<"Node_radius\t"<<Node_radius<<endl;
    Report<<"Spring model\t"<<spring_model<<endl;
    Report<<"Spring coefficient\t"<<Spring_coefficient<<endl;
    Report<<"Damping coefficient\t"<<Damping_coefficient<<endl;
    if (spring_model==1)
    {
        Report<<"Membrane Spring Model:\tFENE"<<endl;
        
    }
    if (spring_model==2)
    {
        Report<<"Membrane Spring Model:\tHoukian"<<endl;
        
    }
    Report<<"Shift_in_X_direction\t"<<Shift_in_X_direction<<endl;
    Report<<"Shift_in_Y_direction\t"<<Shift_in_Y_direction<<endl;
    Report<<"Shift_in_Z_direction\t"<<Shift_in_Z_direction<<endl;
    
}
