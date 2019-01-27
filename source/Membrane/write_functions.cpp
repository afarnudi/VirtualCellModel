//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"
#include "General_constants.h"

void Membrane::write_traj (string traj_name, string label){
    ofstream Trajectory;
    Trajectory.open(traj_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    for(int j=0; j< Num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
    }
}

void Membrane::relaxation_traj (void)
{
    string energy_file_name;
    string traj_file_name;
    
    traj_file_name="Results/Relaxation/Relaxation_"+GenConst::trajectory_file_name+"Membrane_"+to_string(mem_index)+"_"+file_time+".xyz";
    //trajectory:
    
    ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    Trajectory <<Num_of_Nodes<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< Num_of_Nodes; j++) // saving trajectory
    {
        Trajectory << "mem" <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
    }
    
}


void Membrane::export_for_resume(int MD_step){
    ofstream write_resume_file;
    string resume_file_name="Results/Resumes/Resume_Membrane_"+to_string(mem_index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        //Node_force=0
    }
    
    write_resume_file<<Num_of_Triangles<<endl;
    for (int i=0; i<Num_of_Triangles; i++) {
        write_resume_file<<Triangle_list[i][0]<<"\t"<<Triangle_list[i][1]<<"\t"<<Triangle_list[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Node_Pairs<<endl;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        write_resume_file<<Node_Bond_list[i][0]<<"\t"<<Node_Bond_list[i][1]<<"\n";
    }
    //In the import function we should call the neighbour list constructor
    write_resume_file<<Num_of_Triangle_Pairs<<endl;
    for (int i=0; i<Num_of_Triangle_Pairs; i++) {
        write_resume_file<<Triangle_pair_list[i][0]<<"\t"<<Triangle_pair_list[i][1]<<"\n";
        write_resume_file<<Triangle_Pair_Nodes[i][0]<<"\t"<<Triangle_Pair_Nodes[i][1]<<"\t"<<Triangle_Pair_Nodes[i][2]<<"\t"<<Triangle_Pair_Nodes[i][3]<<"\t"<<"\n";
    }
    write_resume_file<<Max_node_pair_length<<"\t"<<Min_node_pair_length<<"\t"<<Average_node_pair_length<<endl;
    //run_check_ for max and min node distances
}


void Membrane::generate_report()
{
    string Report_file_name;
    Report_file_name= "Results/Reports/Report_Membrane_"+to_string(mem_index)+"_";
    Report_file_name+=file_time;
    Report_file_name+=".txt";
    
    ofstream Report;
    Report.open(Report_file_name.c_str());
    Report<< std:: fixed;
    Report<<"General MD Params:\n---------------\n";
    Report<<"MD_num_of_steps"<<setw(20)<<GenConst::MD_num_of_steps<<endl;
    Report<<"MD_traj_save_step"<<setw(20)<<GenConst::MD_traj_save_step<<endl;
    Report<<"MD_Time_Step"<<setw(20)<<GenConst::MD_Time_Step<<endl;
    Report<<"MD_T"<<setw(20)<<GenConst::MD_T<<endl;
    Report<<"MD_thrmo_step"<<setw(20)<<GenConst::MD_thrmo_step<<endl;
    Report<<"MC_step"<<setw(20)<<GenConst::MC_step<<endl;
    Report<<"Mem_fluidity"<<setw(20)<<GenConst::Mem_fluidity<<endl;
    Report<<"Lbox"<<setw(20)<<GenConst::Lbox<<endl;
    Report<<"Periodic_condtion_status"<<setw(20)<<GenConst::Periodic_condtion_status<<endl;
    Report<<"Num_of_Membranes"<<setw(20)<<GenConst::Num_of_Membranes<<endl;
    Report<<"trajectory_file_name"<<setw(20)<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"Membrane Params:\n---------------\n";
    Report<<"Node Mass"<< setw(20)<<Node_Mass<<endl;
    Report<<"Radius"<< setw(20)<<Radius<<endl;
    Report<<"spring_model"<< setw(20)<<spring_model<<endl;
    if (spring_model==1)
    {
        Report<<"Membrane Spring Model:"<< setw(20)<<"FENE"<<endl;
    }
    if (spring_model==2)
    {
        Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;
    }
    Report<<"Spring coefficient"<< setw(20)<<Spring_coefficient<<endl;
    Report<<"Bending coefficient"<< setw(20)<<Bending_coefficient<<endl;
    Report<<"Damping coefficient"<< setw(20)<<Damping_coefficient<<endl;
    Report<<"K_surfaceConstant_local"<< setw(20)<<K_surfaceConstant_local<<endl;
    Report<<"Shift_in_X_direction"<< setw(20)<<Shift_in_X_direction<<endl;
    Report<<"Shift_in_Y_direction"<< setw(20)<<Shift_in_Y_direction<<endl;
    Report<<"Shift_in_Z_direction"<< setw(20)<<Shift_in_Z_direction<<endl;
    Report<<"Downward_speed"<< setw(20)<<Downward_speed<<endl;
    Report<<"X_in_mem"<< setw(20)<<X_in_mem<<endl;
    Report<<"Y_in_mem"<< setw(20)<<Y_in_mem<<endl;
    Report<<"Z_in_mem"<< setw(20)<<Z_in_mem<<endl;
    
    Report<<"Minimum node pair length"<< setw(20)<<Min_node_pair_length<<endl;
    Report<<"Maximum node pair length"<< setw(20)<<Max_node_pair_length<<endl;
    Report<<"Average node pair length"<< setw(20)<<Average_node_pair_length<<endl;
    Report<<"# of Nodes "<< setw(20)<<return_num_of_nodes()<<endl;
    Report<<"# of Triangles "<< setw(20)<<return_num_of_triangle()<<endl;
    
    
    
    
    
}

void Membrane::write_parameters(int MD_Step){
    //    string energy_file_name;
    string traj_file_name;
    omega_calculator();
    double to_T=2.0/(3.0*Num_of_Nodes-3);
    double a[3]={Omega[0],Omega[1],Omega[2]};
    double Omega_len=vector_length(a);
    
    traj_file_name="Results/Param_"+GenConst::trajectory_file_name+"Membrane_"+to_string(mem_index)+"_"+file_time+".txt";
    ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    if (GenConst::File_header==false) {
        Trajectory<<"MD Step\t"<<"Total Kinetic Energy\t"<<"Temperature\t"
        <<"Error\t"<<"Error_com\t"
        <<"omega_x\t"<<"omega_y\t"<<"omega_z\t"<<"omega\t"
        <<"k_Omega\t"<<"delta_k_omega"<<endl;
        GenConst::File_header=true;
    }
    Trajectory<<MD_Step<<"\t"<<Total_Kinetic_Energy<<"\t"<<to_T*Total_Kinetic_Energy
    <<"\t"<<error<<"\t"<<error_com
    <<"\t"<<Omega[0]<<"\t"<<Omega[1]<<"\t"<<Omega[2]<<"\t"<<Omega_len
    <<"\t"<<k_angular<<"\t"<<delta_k_angular<<endl;
}
