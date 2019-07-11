//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"
#include "General_constants.h"

void ECM::write_traj (std::string traj_name, std::string label){
    std::ofstream Trajectory;
    Trajectory.open(traj_name.c_str(), std::ios::app);
    Trajectory << std:: fixed;
    for(int j=0; j< Num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<std::setprecision(5)<< std::setw(20)<<Node_Position[j][0]<< std::setw(20)<<Node_Position[j][1]<< std::setw(20)<<Node_Position[j][2]<<endl;
    }
}





void ECM::generate_report()
{
    std::string Report_file_name;
    Report_file_name= "Results/Reports/Report_ECM_"+std::to_string(index)+"_";
    Report_file_name+=file_time;
    Report_file_name+=".txt";
    
    std::ofstream Report;
    Report.open(Report_file_name.c_str());
    Report<< std:: fixed;
    Report<<"General MD Params:\n---------------\n";
    Report<<"MD_num_of_steps"<<std::setw(20)<<GenConst::MD_num_of_steps<<endl;
    Report<<"MD_traj_save_step"<<std::setw(20)<<GenConst::MD_traj_save_step<<endl;
    Report<<"Step_Size_In_Fs"<<std::setw(20)<<GenConst::Step_Size_In_Fs<<endl;
    Report<<"MD_T"<<std::setw(20)<<GenConst::MD_T<<endl;
    Report<<"MD_thrmo_step"<<std::setw(20)<<GenConst::MD_thrmo_step<<endl;
    Report<<"Bussi_tau"<<std::setw(20)<<GenConst::Bussi_tau<<endl;
    Report<<"MC_step"<<std::setw(20)<<GenConst::MC_step<<endl;
    Report<<"Mem_fluidity"<<std::setw(20)<<GenConst::Mem_fluidity<<endl;
    Report<<"Lbox"<<std::setw(20)<<GenConst::Lbox<<endl;
    Report<<"Periodic_condtion_status"<<std::setw(20)<<GenConst::Periodic_condtion_status<<endl;
    Report<<"Num_of_ECMs"<<std::setw(20)<<GenConst::Num_of_ECMs<<endl;
    Report<<"trajectory_file_name"<<std::setw(20)<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"ECM Params:\n---------------\n";
    Report<<"# of Nodes"<< std::setw(20) <<Num_of_Nodes<<endl;
    Report<<"Node Mass"<< std::setw(20)<<Node_Mass<<endl;
    Report<<"Node Radius"<< std::setw(20)<<Node_radius<<endl;
    Report<<"Number of Node Pairs"<< std::setw(20)<<Num_of_Node_Pairs<<endl;
    Report<<"spring_model"<< std::setw(20)<<spring_model<<endl;
    if (spring_model==0)
    {
        Report<<"ECM Spring Model:"<< std::setw(20)<<"Hookian"<<endl;
    }
    if (spring_model==2)
    {
        Report<<"ECM Spring Model:"<< std::setw(20)<<"Kelvin"<<endl;
    }
    Report<<"Spring coefficient"<< std::setw(20)<<Spring_coefficient<<endl;
    Report<<"Kelvin Damping Coefficient"<< std::setw(20)<<Kelvin_Damping_Coefficient<<endl;
    Report<<"Shift_in_X_direction"<< std::setw(20)<<Shift_in_X_direction<<endl;
    Report<<"Shift_in_Y_direction"<< std::setw(20)<<Shift_in_Y_direction<<endl;
    Report<<"Shift_in_Z_direction"<< std::setw(20)<<Shift_in_Z_direction<<endl;
    Report<<"Downward_speed"<< std::setw(20)<<Downward_speed<<endl;
    
}
