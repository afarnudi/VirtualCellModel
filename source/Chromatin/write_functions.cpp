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

void Chromatin::export_for_resume(int MD_step){
    std::ofstream write_resume_file;
    string resume_file_name="Results/Resumes/Resume_Chromatin_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    
    write_resume_file<<num_of_node_types<<endl;
    for (int i=0; i<num_of_node_types; i++) {
        write_resume_file<<epsilon_LJ[i]<<"\t"<<sigma_LJ[i]<<endl;
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        write_resume_file<<ABC_index[i]<<"\n";
        //Node_force=0
    }
}

void Chromatin::export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count){
    std::ofstream write_resume_file;
    string resume_file_name="Results/Resumes/Resume_Chromatin_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    
    write_resume_file<<num_of_node_types<<endl;
    for (int i=0; i<num_of_node_types; i++) {
        write_resume_file<<epsilon_LJ[i]<<"\t"<<sigma_LJ[i]<<endl;
    }
    
    for (int i=atom_count; i<atom_count+Num_of_Nodes; i++) {
        Node_Position[i-atom_count][0] = atoms[i].posInNm[0];
        Node_Position[i-atom_count][1] = atoms[i].posInNm[1];
        Node_Position[i-atom_count][2] = atoms[i].posInNm[2];
        
        Node_Velocity[i-atom_count][0] = atoms[i].velocityInNmperPs[0];
        Node_Velocity[i-atom_count][1] = atoms[i].velocityInNmperPs[1];
        Node_Velocity[i-atom_count][2] = atoms[i].velocityInNmperPs[2];
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        write_resume_file<<ABC_index[i]<<"\n";
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
    Report<<"MC_step\t"<<GenConst::MC_step<<endl;
    Report<<"Mem_fluidity\t"<<GenConst::Mem_fluidity<<endl;
    Report<<"Lbox\t"<<GenConst::Lbox<<endl;
    Report<<"Periodic_box\t"<<GenConst::Periodic_box<<endl;
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
//    Report<<"Shift_in_X_direction\t"<<Shift_in_X_direction<<endl;
//    Report<<"Shift_in_Y_direction\t"<<Shift_in_Y_direction<<endl;
//    Report<<"Shift_in_Z_direction\t"<<Shift_in_Z_direction<<endl;
    
    Report<<"num of node types = "<<num_of_node_types<<endl;
    Report<<"epsilon_LJ\tsigma_LJ\n--------------------\n";
    for (int i=0; i<num_of_node_types; i++) {
        Report<<std::to_string(i)<<" "<<epsilon_LJ[i]<<"\t"<<sigma_LJ[i]<<endl;
    }
    
}

void Chromatin::export_coordinates(void){
    std::string traj_name="Results/"+GenConst::trajectory_file_name+file_time+"_chromatin_"+std::to_string(index)+".txt";
    
    std::ofstream exportcoords;
    exportcoords.open(traj_name.c_str());
    
    if (!exportcoords.is_open()){
        string errorMessage = TWARN;
        errorMessage+="I can't write the chromatin initial coordinates to storage. \nPath: "+traj_name;
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
    for (int i=0; i<Num_of_Nodes; i++) {
        
        exportcoords<<std::fixed << std::setprecision(8)<<Node_Position[i][0]<<" "<<Node_Position[i][1]<<" "<<Node_Position[i][2]<<" "<<Node_Velocity[i][0]<<" "<<Node_Velocity[i][2]<<" "<<Node_Velocity[i][3]<<"\n";
    }
    exportcoords.close();
}

