//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"
#include "General_constants.h"

using std::endl;
using std::cout;


void Actin::generate_report()
{
    std::string Report_file_name;
    Report_file_name= "Results/Reports/Report_Actin_"+std::to_string(index)+"_";
    Report_file_name+=file_time;
    Report_file_name+=".txt";
    
    std::ofstream Report;
    Report.open(Report_file_name.c_str());
    
    if (!Report.is_open()){
        cout<<"Couldn't generate report. Please check that the required directories exist. The required directory should be ./bin/Results/Reports or ./Results/Reports"<<endl;
        Report.close();
        cout<<"Generating report in an alternative directory:\nbin/Results/Reports/"<<endl;
        Report_file_name= "bin/Results/Reports/Report_Actin_"+std::to_string(index)+"_";
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
    Report<<"Num_of_Membranes\t"<<GenConst::Num_of_Membranes<<endl;
    Report<<"trajectory_file_name\t"<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"Actin Params:\n---------------\n";
    Report<<"Mesh_file_name\t"<<Mesh_file_name<<endl;
    Report<<"Node_Mass\t"<<Node_Mass<<endl;
    Report<<"Node_radius\t"<<Node_radius<<endl;
    Report<<"spring_model\t"<<spring_model<<endl;
    if (spring_model==0)
    {
        Report<<"Actin Spring Model:\tHookian"<<endl;
    }
    if (spring_model==2)
    {
        Report<<"Actin Spring Model:\tKelvin"<<endl;
    }
    Report<<"Spring_coefficient\t"<<Spring_coefficient<<endl;
    Report<<"Kelvin Damping Coefficient\t"<<Kelvin_Damping_Coefficient<<endl;
    Report<<"Shift_in_X_direction\t"<<Shift_in_X_direction<<endl;
    Report<<"Shift_in_Y_direction\t"<<Shift_in_Y_direction<<endl;
    Report<<"Shift_in_Z_direction\t"<<Shift_in_Z_direction<<endl;
    Report<<"rescale_factor\t"<<rescale_factor<<endl;
    Report<<"# of Nodes\t"<<Num_of_Nodes<<endl;
    Report<<"# of Node Pairs\t"<<Num_of_Node_Pairs<<endl;
    Report<<"Contractile_force\t"<<Contractile_force<<endl;
    Report<<"Contractile_k1\t"<<Contractile_k1<<endl;
    Report<<"Contractile_k2\t"<<Contractile_k2<<endl;
    Report<<"Contractile_rmin factor\t"<<Contractile_rmin<<endl;
    Report<<"Contractile_rmax factor\t"<<Contractile_rmax<<endl;
    Report<<"abp_force\t"<<abp_force<<endl;
    Report<<"abp_k1\t"<<abp_k1<<endl;
    Report<<"abp_k2\t"<<abp_k2<<endl;
    Report<<"abp_rmin factor\t"<<abp_rmin<<endl;
    Report<<"abp_rmax factor\t"<<abp_rmax<<endl;
    Report<<"MT_force\t"<<MT_force<<endl;
    Report<<"MT_k1\t"<<MT_k1<<endl;
    Report<<"MT_k2\t"<<MT_k2<<endl;
    Report<<"MT_rmin factor\t"<<MT_rmin<<endl;
    Report<<"MT_rmax factor\t"<<MT_rmax<<endl;
    
}

void Actin::export_for_resume(int MD_step){
    std::ofstream write_resume_file;
    std::string resume_file_name="Results/Resumes/Resume_Actin_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Node_Pairs<<endl;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        write_resume_file<<Node_Bond_list[i][0]<<"\t"<<Node_Bond_list[i][1]<<"\n";
    }
    
    write_resume_file<<Max_node_pair_length<<"\t"<<Min_node_pair_length<<"\t"<<Average_node_pair_length<<endl;
}

void Actin::export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count){
    std::ofstream write_resume_file;
    std::string resume_file_name="Results/Resumes/Resume_Actin_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    
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
    }
    
    write_resume_file<<Num_of_Node_Pairs<<endl;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        write_resume_file<<Node_Bond_list[i][0]<<"\t"<<Node_Bond_list[i][1]<<"\n";
    }
    
    write_resume_file<<Max_node_pair_length<<"\t"<<Min_node_pair_length<<"\t"<<Average_node_pair_length<<endl;
}
